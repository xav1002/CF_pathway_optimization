"""a basic stoichiometric model with thermodynamics."""

# The MIT License (MIT)
#
# Copyright (c) 2013 Weizmann Institute of Science
# Copyright (c) 2018-2020 Institute for Molecular Systems Biology,
# ETH Zurich
# Copyright (c) 2018-2020 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import logging
import sys
import warnings
from typing import Dict, Iterable, Optional, Tuple, Union

import numpy as np
import pandas as pd
from equilibrator_cache import Compound, Reaction
from scipy import stats
from slugify import slugify

from equilibrator_api import Q_, ComponentContribution, R, default_T
from equilibrator_api.phased_reaction import PhasedReaction

from . import Bounds, SBtabDocument, SBtabError, SBtabTable, open_sbtabdoc


class StoichiometricModel(object):
    """A basic stoichiometric model with thermodynamics.

    Designed as a base model for 'Pathway' which also includes flux directions
    and magnitudes.
    """

    MINIMAL_STDEV = 1e-3  # TODO: move to params dict

    def __init__(
        self,
        S: pd.DataFrame,
        compound_dict: Dict[str, Compound],
        reaction_dict: Dict[str, Reaction],
        comp_contrib: Optional[ComponentContribution] = None,
        standard_dg_primes: Optional[Q_] = None,
        dg_sigma: Optional[Q_] = None,
        bounds: Optional[Bounds] = None,
        config_dict: Optional[Dict[str, str]] = None,
    ) -> None:
        """Initialize a StoichiometricModel object.

        Parameters
        ----------
        S : DataFrame
            stoichiometric matrix, where rows are the compound IDs and columns
            are the reaction IDs
        compound_dict : Dict[str, Compound]
            a dictionary of Compound objects, where the keys are their IDs
        reaction_dict : Dict[str, Reaction]
            a dictionary of Reaction objects, where the keys are their IDs
        comp_contrib : ComponentContribution
            a ComponentContribution object
        standard_dg_primes : Quantity, optional
            reaction energies (in kJ/mol)
        dg_sigma : Quantity, optional
            square root of the uncertainty covariance matrix (in kJ/mol)
        bounds : Bounds, optional
            bounds on metabolite concentrations (by default uses the
            "data/cofactors.csv" file in `equilibrator-api`)
        config_dict : dict, optional
            configuration parameters for Pathway analysis
        """
        self.comp_contrib = comp_contrib or ComponentContribution()
        self.water = self.comp_contrib.ccache.water
        self.is_proton = self.comp_contrib.ccache.is_proton
        self.is_water = self.comp_contrib.ccache.is_water
        self.reaction_dict = reaction_dict
        self.compound_dict = compound_dict

        if bounds is None:
            self._bounds = Bounds.get_default_bounds(self.comp_contrib)
        else:
            self._bounds = bounds.copy()

        self.config_dict = config_dict or dict()
        self.configure()

        self.S = S

        # we keep the index to H2O handy, since we often need to remove
        # water from the list of compounds
        compounds = list(self.compounds)

        if self.water in compounds:
            self.idx_water = compounds.index(self.water)
            # make sure the bounds on Water are (1, 1)
            self._bounds.set_bounds(
                compounds[self.idx_water], Q_(1.0, "M"), Q_(1.0, "M")
            )
        else:
            self.idx_water = -1

        # TODO: properly deal with H+ (like water, they should have no
        # effect on the results)

        self.Nc, self.Nr = self.S.shape

        if standard_dg_primes is None:
            assert dg_sigma is None, (
                "If standard_dg_primes are not " "provided, dg_sigma must also be None"
            )
            self.update_standard_dgs()
        else:
            assert standard_dg_primes.check("[energy]/[substance]")
            assert standard_dg_primes.shape == (self.Nr,)
            if dg_sigma is not None:
                # the format for the uncertainty should be a representation of
                # the square root of the covariance. note that we need to have
                # a rank-reduced matrix, so that the number of columns can
                # also be smaller than Nr
                assert dg_sigma.check("[energy]/[substance]")
                assert dg_sigma.shape[0] == self.Nr

            # dGr should be orthogonal to nullspace of S
            # If not, dGr is not contained in image(S) and then there
            # is no consistent set of dGfs that generates dGr and the
            # first law of thermo is violated by the model.
            S_T = self.S.T.values
            S_inv = np.linalg.pinv(S_T)
            null_proj = np.eye(self.S.shape[1]) - S_T @ S_inv
            projected = null_proj @ standard_dg_primes.T
            assert (projected < Q_("1e-8 kJ/mol")).all(), (
                "Supplied reaction standard deltaG values are inconsistent "
                "with the stoichiometric matrix."
            )
            self.standard_dg_primes = standard_dg_primes
            self.dg_sigma = dg_sigma

    def configure(self) -> None:
        """Configure the Component Contribution aqueous conditions."""
        if "p_h" in self.config_dict:
            self.comp_contrib.p_h = Q_(self.config_dict["p_h"])

        if "p_mg" in self.config_dict:
            self.comp_contrib.p_mg = Q_(self.config_dict["p_mg"])

        if "ionic_strength" in self.config_dict:
            self.comp_contrib.ionic_strength = Q_(self.config_dict["ionic_strength"])

        if "temperature" in self.config_dict:
            self.comp_contrib.temperature = Q_(self.config_dict["temperature"])

        # dg_confidence is the confidence we want for the bounds on the
        # standard_dg variables (which is described by a multivariate normal
        # distribution). 0.95 there is a 5% chance that the correct
        # solution is outside the bounds. If you don't want to allow any
        # change to the mean standard_dg estimates, set dg_confidence to 0.
        if "dg_confidence" in self.config_dict:
            self.dg_confidence = float(self.config_dict["dg_confidence"])
        else:
            self.dg_confidence = 0.95

        # ln_conc_confidence is the confidence we want for the bounds on each
        # individual log concnetration variables. 0.95 means there is a 5%
        # chance that the correct solution is outside the bounds (for each
        # single metabolite indipendently). If you don't want to allow any
        # change to the mean standard_dg estimates, set ln_conc_confidence to 0.
        if "ln_conc_confidence" in self.config_dict:
            self.ln_conc_confidence = float(self.config_dict["ln_conc_confidence"])
        else:
            self.ln_conc_confidence = 0.95

    @property
    def compound_ids(self) -> Iterable[str]:
        """Get the list of compound IDs."""
        return self.S.index

    @property
    def compounds(self) -> Iterable[Compound]:
        """Get the list of Compound objects."""
        return map(self.compound_dict.get, self.compound_ids)

    @property
    def compound_df(self) -> pd.DataFrame:
        """Get a DataFrame with all the compound data.

        The columns are:
            compound_id
            lower_bound
            upper_bound

        """
        lbs, ubs = self.bounds
        return pd.DataFrame(
            data=list(
                zip(
                    self.compound_ids,
                    lbs,
                    ubs,
                )
            ),
            columns=[
                "compound_id",
                "lower_bound",
                "upper_bound",
            ],
        )

    @property
    def reaction_ids(self) -> Iterable[str]:
        """Get the list of reaction IDs."""
        return self.S.columns

    @property
    def reactions(self) -> Iterable[Reaction]:
        """Get the list of Reaction objects."""
        return map(self.reaction_dict.get, self.reaction_ids)

    @property
    def reaction_formulas(self) -> Iterable[str]:
        """Iterate through all the reaction formulas.

        :return: the reaction formulas
        """
        reaction_formulas = []
        for rid in self.reaction_ids:
            r_series = self.S[rid]
            sparse = {
                cid: r_series[cid] for cid in self.compound_ids if r_series[cid] != 0
            }
            reaction_formulas.append(StoichiometricModel._sparse_to_formula(sparse))
        return reaction_formulas

    @property
    def reaction_df(self) -> pd.DataFrame:
        """Get a DataFrame with all the reaction data.

        The columns are:
            reaction_id
            reaction_formula
            standard_dg_prime

        """
        return pd.DataFrame(
            data=list(
                zip(
                    self.reaction_ids,
                    self.reaction_formulas,
                    self.standard_dg_primes,
                )
            ),
            columns=[
                "reaction_id",
                "reaction_formula",
                "standard_dg_prime",
            ],
        )

    def update_standard_dgs(self) -> None:
        """Calculate the standard ΔG' values and uncertainties.

        Use the Component Contribution method.
        """
        (
            self.standard_dg_primes,
            self.dg_sigma,
        ) = self.comp_contrib.standard_dg_prime_multi(
            self.reactions,
            uncertainty_representation="fullrank",
        )

    def set_bounds(
        self, cid: str, lb: Optional[Q_] = None, ub: Optional[Q_] = None
    ) -> None:
        """Set the lower and upper bound of a compound.

        Parameters
        ----------
        compound_id: str
            the compound ID
        lb: Quantity, optional
            the new concentration lower bound (ignored if the value is None)
        ub: Quantity, optional
            the new concentration upper bound (ignored if the value is None)
        """

        old_ld, old_ub = self.get_bounds(cid)
        lb = lb or old_ld
        ub = ub or old_ub

        if not lb.check("[concentration]"):
            raise ValueError("Bounds must be in units of concentration")
        if not ub.check("[concentration]"):
            raise ValueError("Bounds must be in units of concentration")

        self._bounds.set_bounds(self.compound_dict[cid], lb, ub)

    def get_bounds(self, cid: str) -> Tuple[Q_, Q_]:
        """Get the lower and upper bound of a compound.

        Parameters
        ----------
        compound_id: str
            the compound ID

        Returns
        -------
        lb: Quantity, optional
            the new concentration lower bound (ignored if the value is None)
        ub: Quantity, optional
            the new concentration upper bound (ignored if the value is None)
        """

        if cid not in self.compound_dict:
            raise KeyError(f"There is no compound with the ID {cid}")

        lb = self._bounds.get_lower_bound(self.compound_dict[cid])
        ub = self._bounds.get_upper_bound(self.compound_dict[cid])
        return lb, ub

    @property
    def bounds(self) -> Tuple[Iterable[Q_], Iterable[Q_]]:
        """Get the concentration bounds.

        The order of compounds is according to the stoichiometric matrix index.

        Returns
        -------
        tuple of (lower bounds, upper bounds)
        """
        return self._bounds.get_bounds(self.compounds)

    @property
    def bound_df(self) -> pd.DataFrame:
        """Get a DataFrame with all the bounds data."""
        return self._bounds.to_data_frame()

    @property
    def ln_conc_lb(self) -> np.array:
        """Get the log lower bounds on the concentrations.

        The order of compounds is according to the stoichiometric matrix index.

        Returns
        -------
        a NumPy array of the log lower bounds
        """
        return np.array(
            list(map(float, self._bounds.get_ln_lower_bounds(self.compounds)))
        )

    @property
    def ln_conc_ub(self) -> np.ndarray:
        """Get the log upper bounds on the concentrations.

        The order of compounds is according to the stoichiometric matrix index.

        Returns
        -------
        a NumPy array of the log upper bounds
        """
        return np.array(
            list(map(float, self._bounds.get_ln_upper_bounds(self.compounds)))
        )

    @property
    def ln_conc_mu(self) -> np.array:
        """Get mean of log concentration distribution based on the bounds.

        The order of compounds is according to the stoichiometric matrix index.

        Returns
        -------
        a NumPy array with the mean of the log concentrations
        """
        return (self.ln_conc_ub + self.ln_conc_lb) / 2.0

    @property
    def ln_conc_sigma(self) -> np.array:
        """Get stdev of log concentration distribution based on the bounds.

        Returns
        -------
        a NumPy array with the stdev of the log concentrations
        """
        ci = (self.ln_conc_ub - self.ln_conc_lb) / 2.0
        return ci / stats.chi2.ppf(self.ln_conc_confidence, 1) ** (0.5)

    @staticmethod
    def _sparse_to_formula(sparse: Dict[str, float]) -> str:
        """Write a reaction as a chemical formula."""

        def write_compound_and_coeff(cid: str, coeff: float) -> str:
            if np.abs(coeff - 1.0) < sys.float_info.epsilon:
                return cid
            else:
                return "%g %s" % (coeff, cid)

        l_tokens = []
        r_tokens = []
        for cid, coeff in sorted(sparse.items()):
            if coeff < 0:
                l_tokens.append(write_compound_and_coeff(cid, -coeff))
            elif coeff > 0:
                r_tokens.append(write_compound_and_coeff(cid, coeff))

        l_str = " + ".join(l_tokens)
        r_str = " + ".join(r_tokens)
        return f"{l_str} = {r_str}"

    @staticmethod
    def read_thermodynamics(
        thermo_sbtab: SBtabTable, config_dict: Dict[str, str]
    ) -> Dict[str, Q_]:
        """Read the 'thermodynamics' table from an SBtab.

        Parameters
        ----------
        thermo_sbtab : SBtabTable
            A SBtabTable containing the thermodynamic data
        config_dict : dict
            A dictionary containing the configuration arguments

        Returns
        -------
        A dictionary mapping reaction IDs to standard ΔG' values.
        """
        try:
            std_conc = thermo_sbtab.get_attribute("StandardConcentration")
            assert Q_(std_conc) == Q_(
                "1 M"
            ), "We only support a standard concentration of 1 M."
        except SBtabError:
            pass

        if "temperature" in config_dict:
            temperature = Q_(config_dict["temperature"])
            assert temperature.check("[temperature]")
        else:
            temperature = default_T

        thermo_df = thermo_sbtab.to_data_frame()
        cols = ["QuantityType", "Value", "Compound", "Reaction", "Unit"]

        rid2dg0 = dict()
        for i, row in thermo_df.iterrows():
            typ, val, cid, rid, unit = [row[c] for c in cols]
            val = Q_(float(val), unit)

            if typ == "equilibrium constant":
                if not val.check(""):
                    raise ValueError(
                        f"Error in parameter table row {i}: "
                        "equilibrium constants must have no units.\n" + row
                    )
                rid2dg0[rid] = -np.log(val.m_as("")) * R * temperature
            elif typ == "reaction gibbs energy":
                if not val.check("[energy]/[substance]"):
                    raise ValueError(
                        f"Error in parameter table row {i}: "
                        "Gibbs energies must be in kJ/mol (or equivalent).\n" + row
                    )
                val.ito("kJ/mol")
                rid2dg0[rid] = val
            else:
                raise ValueError(
                    f"Error in parameter table row {i}: "
                    "unrecognized Rate Constant Type ({typ}).\n" + row
                )
        return rid2dg0

    @staticmethod
    def _read_network_sbtab(
        filename: Union[str, SBtabDocument],
        comp_contrib: Optional[ComponentContribution] = None,
        freetext: bool = True,
    ):
        """Initialize a Pathway object using a 'network'-only SBtab.

        Parameters
        ----------
        filename : str, SBtabDocument
             a filename containing an SBtabDocument (or the SBtabDocument
             object itself) defining the network (topology) only
        comp_contrib : ComponentContribution, optional
             a ComponentContribution object needed for parsing and searching
             the reactions.
             also used to set the aqueous parameters (pH, I, etc.)
        freetext : bool, optional
             a flag indicating whether the reactions are given as free-text
             (i.e. common names for compounds) or by standard database
             accessions. (Default value: `True`)

        Returns
        -------
        sbtabdoc
        S
        reaction_dict
        compound_dict
        config_dict
        comp_contrib
        """

        def formula_to_obj(rxn_formula: str) -> Reaction:
            """Convert the reaction formulas to Reaction objects."""
            logging.debug("formula = %f", rxn_formula)

            if freetext:
                rxn = comp_contrib.search_reaction(rxn_formula)
            else:
                rxn = comp_contrib.parse_reaction_formula(rxn_formula)
            if not rxn.is_balanced():
                raise Exception(f"This reaction is not balanced: {rxn_formula}")
            return rxn

        comp_contrib = comp_contrib or ComponentContribution()

        sbtabdoc = open_sbtabdoc(filename)

        reaction_df = sbtabdoc.get_sbtab_by_id("Reaction").to_data_frame()
        reaction_df["ReactionObject"] = reaction_df.ReactionFormula.apply(
            formula_to_obj
        )
        if "ID" in reaction_df.columns:
            reaction_df.index = reaction_df.ID
        else:
            reaction_df.index = [f"R{i:05d}" for i in range(reaction_df.shape[0])]
        for row in reaction_df.itertuples():
            row.ReactionObject.rid = row.Index

        config_dict = {
            "p_h": str(comp_contrib.p_h),
            "p_mg": str(comp_contrib.p_mg),
            "ionic_strength": str(comp_contrib.ionic_strength),
            "temperature": str(comp_contrib.temperature),
        }
        reaction_dict = reaction_df.ReactionObject.to_dict()

        S = comp_contrib.create_stoichiometric_matrix_from_reaction_objects(
            reaction_df.ReactionObject
        )
        S.columns = reaction_df.index

        compound_dict = {}
        compound_ids = []
        for cpd in S.index:
            cid = slugify(cpd.get_common_name(), separator="_", lowercase=False)
            compound_ids.append(cid)
            compound_dict[cid] = cpd
        S.index = compound_ids
        return (
            sbtabdoc,
            S,
            reaction_dict,
            compound_dict,
            config_dict,
            comp_contrib,
        )

    @staticmethod
    def _read_model_sbtab(
        filename: Union[str, SBtabDocument],
        comp_contrib: Optional[ComponentContribution] = None,
    ):
        """Parse and SBtabDocument and return a StoichiometricModel.

        Parameters
        ----------
        filename : str or SBtabDocument
            a filename containing an SBtabDocument (or the SBtabDocument
            object itself) defining the pathway
        comp_contrib : ComponentContribution, optional
            a ComponentContribution object needed for parsing and searching
            the reactions. also used to set the aqueous parameters (pH, I, etc.)

        Returns
        -------
        stoich_model: StoichiometricModel
            A StoichiometricModel object based on the configuration SBtab
        S
        compound_dict
        reaction_didct
        standard_dg_primes
        dg_sigma
        bounds
        config_dict
        """
        sbtabdoc = open_sbtabdoc(filename)

        # Read the configuration table if it exists
        config_sbtab = sbtabdoc.get_sbtab_by_id("Configuration")
        if config_sbtab:
            config_df = config_sbtab.to_data_frame()
            assert (
                "Option" in config_df.columns
            ), "Configuration table must have an Option column"
            assert (
                "Value" in config_df.columns
            ), "Configuration table must have an Value column"
            config_dict = config_df.set_index("Option").Value.to_dict()
        else:
            config_dict = dict()

        table_ids = ["Compound", "ConcentrationConstraint", "Reaction"]
        dfs = []

        for table_id in table_ids:
            sbtab = sbtabdoc.get_sbtab_by_id(table_id)
            if sbtab is None:
                tables = ", ".join(map(lambda s: s.table_id, sbtabdoc.sbtabs))
                raise ValueError(
                    f"The SBtabDocument must have a table "
                    f"with the following ID: {table_id}, "
                    f"however, only these tables were "
                    f"found: {tables}"
                )
            dfs.append(sbtab.to_data_frame())

        compound_df, bounds_df, reaction_df = dfs

        # remove trailing whitespaces from the ID columns in the compound and
        # reaction tables, and also in the 'Compound' column of the bounds table
        # which points to compound IDs
        compound_df["ID"] = compound_df["ID"].str.strip()
        reaction_df["ID"] = reaction_df["ID"].str.strip()
        bounds_df["Compound"] = bounds_df["Compound"].str.strip()

        # Read the Compound table
        # -----------------------
        # use equilibrator-cache to build a dictionary of compound IDs to
        # Compound objects
        compound_df.set_index("ID", inplace=True)
        reaction_df.set_index("ID", inplace=True)

        if "Identifiers" not in compound_df.columns:
            # Fallback legacy option, KEGG IDs could be provided without the
            # namespace when the column title explicitly indicates that they
            # were from KEGG
            if "Identifiers:kegg.compound" in compound_df.columns:
                compound_df["Identifiers"] = (
                    "kegg:" + compound_df["Identifiers:kegg.compound"]
                )
            else:
                raise KeyError(
                    "There is no column of Identifiers in the Compound table"
                )
        compound_df["Identifiers"] = compound_df["Identifiers"].str.strip()

        comp_contrib = comp_contrib or ComponentContribution()

        compound_objects = []
        for i, row in enumerate(compound_df.itertuples()):
            cpd = None
            if row.Identifiers:
                cpd = comp_contrib.get_compound(row.Identifiers)
            if cpd is None:
                # generate a novel Compound object (with a unique negative ID)
                # just to use it later in this specific pathway
                cpd = Compound(id=-i)
            ph_cpd = PhasedReaction.to_phased_compound(cpd)
            if "RedoxPotential" in compound_df.columns and row.RedoxPotential:
                ph_cpd.abundance = Q_(row.RedoxPotential)
            compound_objects.append(ph_cpd)

        compound_df["CompoundObject"] = compound_objects

        if pd.isnull(compound_df.CompoundObject).any():
            accessions_not_found = compound_df.loc[
                pd.isnull(compound_df.CompoundObject), "Identifiers"
            ]
            error_msg = str(accessions_not_found.to_dict())
            raise KeyError(
                f"Some compounds not found in equilibrator-cache: " f"{error_msg}"
            )

        # Read the ConcentrationConstraints table
        # ---------------------------------------
        # convert compound IDs in the bounds table to Compound objects
        bounds_df = bounds_df.join(compound_df[["CompoundObject"]], on="Compound")
        bounds_df.set_index("CompoundObject", inplace=True)

        try:
            bounds_sbtab = sbtabdoc.get_sbtab_by_id("ConcentrationConstraint")
            bounds_unit = bounds_sbtab.get_attribute("Unit")
            lbs = bounds_df["Min"].apply(lambda x: Q_(float(x), bounds_unit))
            ubs = bounds_df["Max"].apply(lambda x: Q_(float(x), bounds_unit))
        except SBtabError:
            # if the unit is not defined in the header, we assume it is given
            # next to each bound individually
            lbs = bounds_df["Min"].apply(Q_)
            ubs = bounds_df["Max"].apply(Q_)

        bounds = Bounds(lbs.to_dict(), ubs.to_dict())
        bounds.check_bounds()

        # Read the Reaction table
        # -----------------------
        S = comp_contrib.create_stoichiometric_matrix_from_reaction_formulas(
            reaction_df.ReactionFormula
        )
        S.columns = reaction_df.index

        missing_ids = set(S.index).difference(compound_df.index)

        assert len(missing_ids) == 0, (
            "Some compounds used in the `Reaction` table are not present in "
            "the `Compound` table: " + ", ".join(missing_ids)
        )

        # TODO: if H+ is one of the compounds in S, drop it
        # TODO: if water is not one of the compounds in S, add another row

        reaction_objects = []
        for rid in reaction_df.index:
            sparse = dict()
            sparse_with_phases = dict()
            for cid in S.index:
                if S.at[cid, rid] != 0:
                    cpd = compound_df.at[cid, "CompoundObject"]
                    sparse[cpd.compound] = S.at[cid, rid]
                    sparse_with_phases[cpd] = S.at[cid, rid]
            rxn = PhasedReaction(
                sparse=sparse,
                arrow="=",
                rid=rid,
                sparse_with_phases=sparse_with_phases,
            )
            if not rxn.is_balanced(ignore_atoms=("H",)):
                warnings.warn(f"Reaction {rxn.rid} is not balanced", stacklevel=3)
            reaction_objects.append(rxn)
        reaction_df["ReactionObject"] = reaction_objects

        thermo_sbtab = sbtabdoc.get_sbtab_by_id("Thermodynamics")
        if thermo_sbtab:
            rid2dg0 = StoichiometricModel.read_thermodynamics(thermo_sbtab, config_dict)

            # make sure all reactions have a standard Gibbs energy
            missing_ids = set(S.columns).difference(rid2dg0.keys())
            assert len(missing_ids) == 0, (
                "Some reactions are missing a Keq or ΔG'0 "
                "in the `Thermodynamics` table: " + ", ".join(missing_ids)
            )
            standard_dg_primes = np.array(
                [rid2dg0[rid].m_as("kJ/mol") for rid in reaction_df.index]
            ) * Q_("1 kJ/mol")
            dg_sigma = None
        else:
            standard_dg_primes = None
            dg_sigma = None

        return (
            sbtabdoc,
            S,
            compound_df.CompoundObject.to_dict(),
            reaction_df.ReactionObject.to_dict(),
            standard_dg_primes,
            dg_sigma,
            bounds,
            config_dict,
        )

    @classmethod
    def from_network_sbtab(
        cls,
        filename: Union[str, SBtabDocument],
        comp_contrib: Optional[ComponentContribution] = None,
        freetext: bool = True,
        bounds: Optional[Bounds] = None,
    ) -> "StoichiometricModel":
        """Initialize a Pathway object using a 'network'-only SBtab.

        Parameters
        ----------
        filename : str, SBtabDocument
            a filename containing an SBtabDocument (or the SBtabDocument
            object itself) defining the network (topology) only
        comp_contrib : ComponentContribution, optional
            a ComponentContribution object needed for parsing and searching
            the reactions. also used to set the aqueous parameters (pH, I, etc.)
        freetext : bool, optional
            a flag indicating whether the reactions are given as free-text (i.e.
            common names for compounds) or by standard database accessions
            (Default value: `True`)
        bounds : Bounds, optional
            bounds on metabolite concentrations (by default uses the
            "data/cofactors.csv" file in `equilibrator-api`)

        Returns
        -------
            a Pathway object
        """
        (
            sbtabdoc,
            S,
            reaction_dict,
            compound_dict,
            config_dict,
            comp_contrib,
        ) = StoichiometricModel._read_network_sbtab(filename, comp_contrib, freetext)

        stoich_model = StoichiometricModel(
            S=S,
            compound_dict=compound_dict,
            reaction_dict=reaction_dict,
            comp_contrib=comp_contrib,
            bounds=bounds,
            config_dict=config_dict,
        )
        return stoich_model

    @classmethod
    def from_sbtab(
        cls,
        filename: Union[str, SBtabDocument],
        comp_contrib: Optional[ComponentContribution] = None,
    ) -> "StoichiometricModel":
        """Parse and SBtabDocument and return a StoichiometricModel.

        Parameters
        ----------
        filename : str or SBtabDocument
            a filename containing an SBtabDocument (or the SBtabDocument
            object itself) defining the pathway
        comp_contrib : ComponentContribution, optional
            a ComponentContribution object needed for parsing and searching
            the reactions. also used to set the aqueous parameters (pH, I, etc.)

        Returns
        -------
        stoich_model: StoichiometricModel
            A StoichiometricModel object based on the configuration SBtab

        """
        (
            sbtabdoc,
            S,
            compound_dict,
            reaction_dict,
            standard_dg_primes,
            dg_sigma,
            bounds,
            config_dict,
        ) = cls._read_model_sbtab(filename, comp_contrib)

        return StoichiometricModel(
            S=S,
            compound_dict=compound_dict,
            reaction_dict=reaction_dict,
            comp_contrib=comp_contrib,
            standard_dg_primes=standard_dg_primes,
            dg_sigma=None,
            bounds=bounds,
            config_dict=config_dict,
        )

    def to_sbtab(self) -> SBtabDocument:
        """Export the model to an SBtabDocument."""
        sbtabdoc = SBtabDocument("pathway", filename="pathway.tsv")

        config_dict = dict(self.config_dict)
        config_df = pd.DataFrame(
            data=[(k, v, "") for k, v in config_dict.items()],
            columns=["!Option", "!Value", "!Comment"],
        )
        config_sbtab = SBtabTable.from_data_frame(
            df=config_df, table_id="Configuration", table_type="Config"
        )

        reaction_df = pd.DataFrame(
            data=list(zip(self.reaction_ids, self.reaction_formulas)),
            columns=["!ID", "!ReactionFormula"],
        )
        reaction_sbtab = SBtabTable.from_data_frame(
            df=reaction_df, table_id="Reaction", table_type="Reaction"
        )

        compound_df = pd.DataFrame(
            data=[
                (
                    cid,
                    self.compound_dict[cid].get_common_name(),
                    self.compound_dict[cid].get_accession(),
                )
                for cid in self.compound_ids
            ],
            columns=["!ID", "!Name", "!Identifiers"],
        )
        compound_sbtab = SBtabTable.from_data_frame(
            df=compound_df, table_id="Compound", table_type="Compound"
        )

        thermo_df = pd.DataFrame(
            data=[
                (
                    "reaction gibbs energy",
                    rxn.rid,
                    "",
                    f"{dg.m_as('kJ/mol'):.2f}",
                    "kJ/mol",
                )
                for rxn, dg in zip(self.reactions, self.standard_dg_primes)
            ],
            columns=[
                "!QuantityType",
                "!Reaction",
                "!Compound",
                "!Value",
                "!Unit",
            ],
        )
        thermo_sbtab = SBtabTable.from_data_frame(
            df=thermo_df, table_id="Thermodynamics", table_type="Quantity"
        )
        thermo_sbtab.change_attribute("StandardConcentration", "M")

        lbs, ubs = self.bounds
        lbs = map(lambda x: x.m_as("mM"), lbs)
        ubs = map(lambda x: x.m_as("mM"), ubs)
        conc_df = pd.DataFrame(
            data=[
                (
                    "concentration",
                    cid,
                    f"{lb:.3g}",
                    f"{ub:.3g}",
                )
                for cid, lb, ub in zip(self.compound_ids, lbs, ubs)
            ],
            columns=["!QuantityType", "!Compound", "!Min", "!Max"],
        )
        conc_sbtab = SBtabTable.from_data_frame(
            df=conc_df,
            table_id="ConcentrationConstraint",
            table_type="Quantity",
        )
        conc_sbtab.change_attribute("Unit", "mM")

        sbtabdoc.add_sbtab(config_sbtab)
        sbtabdoc.add_sbtab(reaction_sbtab)
        sbtabdoc.add_sbtab(compound_sbtab)
        sbtabdoc.add_sbtab(thermo_sbtab)
        sbtabdoc.add_sbtab(conc_sbtab)
        return sbtabdoc

    def write_sbtab(self, filename: str) -> None:
        """Write the pathway to an SBtab file."""
        self.to_sbtab().write(filename)
