"""A biochemical compound module."""

# The MIT License (MIT)
#
# Copyright (c) 2018 Institute for Molecular Systems Biology, ETH Zurich
# Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability,
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


import re
from functools import reduce
from typing import Iterable, List, Optional, Union

import numpy as np
from scipy.special import logsumexp
from sqlalchemy import Column, Float, Integer, PickleType, String
from sqlalchemy.orm import Mapped, relationship

from ..exceptions import MissingDissociationConstantsException
from ..thermodynamic_constants import (
    LOG10,
    Q_,
    R,
    debye_hueckel_derivative,
    default_pMg,
    ureg,
)
from . import Base
from .mixins import TimeStampMixin


class Compound(TimeStampMixin, Base):
    """
    Model a chemical compound in the context of component contribution.

    Attributes
    ----------
    id : int
        The primary key in the table.
    inchi_key : str
        InChIKey is a hash of the full InChI with a constant length.
    inchi : str
        InChI descriptor of the molecule.
    smiles : str
        SMILES descriptor of the molecule, taken from MetaNetX but not used.
    mass : float
        Molecualr mass of the molecule.
    atom_bag : dict
        The chemical formula, where keys are the atoms and values are the
        stoichiometric coefficient.
    dissociation_constants : list
        A list of float, which are the pKa values of this molecule.
    group_vector : list
        A list of groups and their counts
    magnesium_dissociation_constants : list
        The compound's MagnesiumDissociationConstants in a one-to-many
        relationship.
    microspecies : list
        The compound's microspecies in a one-to-many relationship
    identifiers : list
        The compound's identifiers in a one-to-many relationship.

    """

    __tablename__ = "compounds"

    # SQLAlchemy column descriptors.
    id: int = Column(Integer, primary_key=True, autoincrement=True)
    inchi_key: Optional[str] = Column(String(), default=None, nullable=True, index=True)
    inchi: Optional[str] = Column(String(), default=None, nullable=True, index=True)
    smiles: Optional[str] = Column(String(), default=None, nullable=True, index=True)
    mass: Optional[float] = Column(Float, default=None, nullable=True)
    atom_bag: Optional[dict] = Column(PickleType, nullable=True)
    dissociation_constants: Optional[list] = Column(PickleType, nullable=True)
    group_vector: Optional[list] = Column(PickleType, nullable=True)
    magnesium_dissociation_constants: Mapped[List["MagnesiumDissociationConstant"]] = (
        relationship(
            "MagnesiumDissociationConstant",
            cascade="all, delete-orphan",
            lazy="select",
        )
    )
    microspecies: Mapped[List["CompoundMicrospecies"]] = relationship(
        "CompoundMicrospecies", cascade="all, delete-orphan", lazy="select"
    )
    identifiers: Mapped[List["CompoundIdentifier"]] = relationship(
        "CompoundIdentifier", cascade="all, delete-orphan", lazy="select"
    )

    def __repr__(self) -> str:
        """Return a string representation of this object."""
        return f"{type(self).__name__}(id={self.id}, inchi_key={self.inchi_key})"

    def __eq__(self, other: "Compound") -> bool:
        """Compare two compound objects, return True if they are equal."""
        return self.id == other.id

    def __lt__(self, other: "Compound") -> bool:
        """Compare two compound objects.

        Return True if the first has a smaller ID.
        """
        return self.id < other.id

    def __hash__(self) -> int:
        """Return a hash for this compound (same as its internal ID)."""
        return self.id

    @property
    def formula(self) -> Optional[str]:
        """Return the chemical formula."""
        if self.atom_bag is None:
            return None

        return "".join(
            [
                element if count == 1 else f"{element}{count}"
                for element, count in sorted(self.atom_bag.items())
                if (count > 0 and element != "e-")
            ]
        )

    def can_be_transformed(self) -> bool:
        """Check if this compound can be transformed.

        In other words, it has been analyzed by ChemAxon and its microspecies
        are populated

        Returns
        -------
        True if the compound can be transformed
        """
        return len(self.microspecies) != 0 and None not in self.microspecies

    def _get_ms_ddg_over_rt(
        self,
        p_h: Q_,
        ionic_strength: Q_,
        temperature: Q_,
        p_mg: Q_ = default_pMg,
    ) -> Iterable:
        if not self.can_be_transformed():
            raise MissingDissociationConstantsException(
                f"{self} has not yet been analyzed by ChemAxon."
            )

        pH = p_h.m_as("")
        pMg = p_mg.m_as("")
        ionic_strength_M = ionic_strength.m_as("M")
        T_in_K = temperature.m_as("K")
        return map(
            lambda ms: -ms.transform(
                pH=pH, pMg=pMg, ionic_strength_M=ionic_strength_M, T_in_K=T_in_K
            ),
            self.microspecies,
        )

    @ureg.check(None, "", "[concentration]", "[temperature]", "")
    def transform(
        self,
        p_h: Q_,
        ionic_strength: Q_,
        temperature: Q_,
        p_mg: Q_ = default_pMg,
    ) -> Q_:
        r"""Calculate the Legendre transform for a compound.

        Parameters
        ----------
        p_h : Quantity
            The pH value, i.e., the logarithmic scale for the molar
            concentration of hydrogen ions :math:`-log10([H+])`
        ionic_strength : Quantity
            Set the ionic strength
        temperature : Quantity
            Set the temperature
        p_mg : Quantity, optional
            The logarithmic molar concentration of magnesium ions
            :math:`-log10([Mg2+])`, (Default value = default_pMg)

        Returns
        -------
        Quantity
            The transformed relative :math:`\Delta G` of this compound.

        """

        # Since the ddg_over_rt values of different microspecies can span
        # a large range, summing their exponents directly will cause floating
        # point overflow issue. Therefore, we use logaddexp instead.
        return (
            -R
            * temperature
            * reduce(
                np.logaddexp,
                self._get_ms_ddg_over_rt(
                    p_h=p_h,
                    p_mg=p_mg,
                    ionic_strength=ionic_strength,
                    temperature=temperature,
                ),
            )
        )

    @ureg.check(None, "", "[concentration]", "[temperature]", "")
    def sensitivity_to_p_h(
        self,
        p_h: Q_,
        ionic_strength: Q_,
        temperature: Q_,
        p_mg: Q_ = default_pMg,
    ) -> Q_:
        r"""Calculate the derivative of the formation ΔG' w.r.t. pH.

        By derivating the Legendre transform as a function of pH, we get
        the slope which is a linear approximation of the pH response.

        Parameters
        ----------
        p_h : Quantity
            The pH value, i.e., the logarithmic scale for the molar
            concentration of hydrogen ions :math:`-log10([H+])`
        ionic_strength : Quantity
            Set the ionic strength
        temperature : Quantity
            Set the temperature
        p_mg : Quantity, optional
            The logarithmic molar concentration of magnesium ions
            :math:`-log10([Mg2+])`, (Default value = default_pMg)

        Returns
        -------
        Quantity
            The derivative of  :math:`\Delta G_f` with respect to pH.

        """

        # The Legendre transform term is:
        # T = -ln( sum[ exp(-g_i) ] ) = -lse[-g_i]
        # where we use lse[] as a shorthand for log-sum-exp,
        # and g_i as the transform term for microspecies i (in RT).
        # Note that the values returned by ms_ddg_over_rt_list() are already
        # including the minus sign.
        #
        # Therefore, the derivative of the entire transform w.r.t ph is:
        # dT/dph = -sum[ (-dg_i/dph) * exp(-g_i) ] / sum[ exp(-g_i) ]
        #        = sum[ exp( -g_i + ln(dg_i/dph) ) ] / sum[ exp(-g_i) ]
        #        = exp( lse[ -g_i + ln(dg_i/dph) ] - lse[ -g_i ] )

        # g_i has only one term that depends on pH which is:
        # N_H_i * ln(10) * pH
        # where N_H_i is the number of protons in microspecies i.
        # Therefore: dg_i/dph = N_H_i * ln(10)
        #
        # So the final expression for the sensitivity is:
        # dT/dph = -exp( lse[ -g_i + ln(N_H_i * ln(10)) ] - lse[ -g_i ] )
        #
        # the expression inside the ln() is always positive, so we don't need
        # to consider a case where the sign of the logsumexp is negative

        ms_ddg_over_rt_list = np.array(
            list(
                self._get_ms_ddg_over_rt(
                    p_h=p_h,
                    p_mg=p_mg,
                    ionic_strength=ionic_strength,
                    temperature=temperature,
                )
            )
        )
        weights = np.array([ms.number_protons for ms in self.microspecies]) * LOG10
        unweighted_sum = logsumexp(ms_ddg_over_rt_list)
        weighted_sum = logsumexp(ms_ddg_over_rt_list, b=weights, return_sign=False)
        return R * temperature * np.exp(weighted_sum - unweighted_sum)

    @ureg.check(None, "", "[concentration]", "[temperature]", "")
    def sensitivity_to_I(
        self,
        p_h: Q_,
        ionic_strength: Q_,
        temperature: Q_,
        p_mg: Q_ = default_pMg,
    ) -> Q_:
        r"""Calculate the derivative of the formation ΔG' w.r.t. pH.

        By derivating the Legendre transform as a function of pH, we get
        the slope which is a linear approximation of the pH response.

        Parameters
        ----------
        p_h : Quantity
            The pH value, i.e., the logarithmic scale for the molar
            concentration of hydrogen ions :math:`-log10([H+])`
        ionic_strength : Quantity
            Set the ionic strength
        temperature : Quantity
            Set the temperature
        p_mg : Quantity, optional
            The logarithmic molar concentration of magnesium ions
            :math:`-log10([Mg2+])`, (Default value = default_pMg)

        Returns
        -------
        Quantity
            The derivative of  :math:`\Delta G_f` with respect to pH.

        """

        # See the derivation in sensitivity_to_p_h(). So we can similary
        # say that:
        # dT/dI = exp( lse[ -g_i + ln(dg_i/dI) ] - lse[ -g_i ] )

        # g_i has only one term that depends on I which is:
        # DH(I) * (z_i^2 - N_H_i - 4 * N_Mg_i)
        # where DH() is the Debye-Hueckel() term and
        # z_i, N_H_i, and N_Mg_i are the charge,
        # number of protons, and number of Mg2+ ions in microspecies i.
        # Therefore: dg_i/dI = dDH/dI * (z_i^2 - N_H_i - 4 * N_Mg_i)
        #
        # So the final expression for the sensitivity is:
        # dT/dph = -exp(
        #               lse[ -g_i + ln(dDH/dI*(z_i^2 - N_H_i - 4 * N_Mg_i)) ] -
        #               lse[ -g_i ]
        #           )
        #
        # the expression inside the ln() can sometimes be negative, so we do need
        # to consider a case where the sign of the logsumexp is negative. However,
        # this is solved by using the `return_sign=True` option.

        ms_ddg_over_rt_list = np.array(
            list(
                self._get_ms_ddg_over_rt(
                    p_h=p_h,
                    p_mg=p_mg,
                    ionic_strength=ionic_strength,
                    temperature=temperature,
                )
            )
        )
        weights = np.array(
            [
                ms.charge**2 - ms.number_protons - 4 * ms.number_magnesiums
                for ms in self.microspecies
            ]
        ) * debye_hueckel_derivative(ionic_strength.m_as("M"), temperature.m_as("K"))
        unweighted_sum = logsumexp(ms_ddg_over_rt_list)
        weighted_sum, sign = logsumexp(ms_ddg_over_rt_list, b=weights, return_sign=True)
        return -sign * np.exp(weighted_sum - unweighted_sum) * Q_("1 kJ/mol/molar")

    def get_common_name(self) -> Union[str, None]:
        """Return the most common name for this compound.

        Right now, we simply choose the first identifier with the smallest
        ID in the list, which is not ideal.
        """
        best_priority = np.inf
        synonym = None
        for identifier in self.identifiers:
            if identifier.registry.name == "Synonyms":
                priority = identifier.id  # TODO: find a better way to choose
                if priority < best_priority:
                    synonym = identifier.accession.split("|")[0]
                    best_priority = priority
        return synonym

    ORDER_OF_REGISTRIES = (
        "MIR:00000567",  # MetaNetX chemical
        "MIR:00000578",  # KEGG
        "MIR:00000556",  # BiGG
        "MIR:00000002",  # ChEBI
        "MIR:00000552",  # SEED
    )

    @staticmethod
    def _identifier_sorting_key(identifier) -> tuple:
        try:
            priority = Compound.ORDER_OF_REGISTRIES.index(
                identifier.registry.identifier
            )
        except ValueError:
            return np.inf, np.inf

        # try to find a the numeric value in the accession
        numbers = re.findall(r"(\d+)", identifier.accession)
        try:
            return priority, int(numbers[0])
        except IndexError:
            # since we lack any other options, just use the length of the
            # accession string
            return priority, len(identifier.accession)

    def get_accession(self) -> Union[str, None]:
        """Get an accession code for this compound."""
        try:
            best_id = min(self.identifiers, key=Compound._identifier_sorting_key)
            return best_id.registry.namespace + ":" + best_id.accession
        except ValueError:
            return None
