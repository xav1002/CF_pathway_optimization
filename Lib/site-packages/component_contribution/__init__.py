# The MIT License (MIT)
#
# Copyright (c) 2023 The Weizmann Institute of Science.
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

from equilibrator_cache import ZenodoSettings

# Chemical group definitions for the component-contribution method
GROUP_DEFINITIONS_SETTINGS = ZenodoSettings(
    doi="10.5281/zenodo.4010930",
    filename="group_summary.csv",
    md5="c3279061cfba3b2370350744e84bbad1",
    url="https://zenodo.org/api/",
)

# Training data for the component contribution method
TRAINING_DATA_TECR_SETTINGS = ZenodoSettings(
    doi="10.5281/zenodo.3978440",
    filename="TECRDB.csv",
    md5="b39e950384004d6aa6dbca60a94c72c1",
    url="https://zenodo.org/api/",
)

TRAINING_DATA_REDOX_SETTINGS = ZenodoSettings(
    doi="10.5281/zenodo.3978440",
    filename="redox.csv",
    md5="c78c03f262f430c16760cc6032800de2",
    url="https://zenodo.org/api/",
)

TRAINING_DATA_FORMATION_SETTINGS = ZenodoSettings(
    doi="10.5281/zenodo.3978440",
    filename="formation_energies_transformed.csv",
    md5="a8b9f3c29fc0fb5f6588d470a5330927",
    url="https://zenodo.org/api/",
)

# Model parameters trained using the Component Contribution method
DEFAULT_CC_PARAMS_SETTINGS = ZenodoSettings(
    doi="10.5281/zenodo.4013789",
    filename="cc_params.npz",
    md5="d20c72a47fe934228a6be9f4e49f96e4",
    url="https://zenodo.org/api/",
)

# Model parameters trained using the Component Contribution method (legacy version)
LEGACY_CC_PARAMS_SETTINGS = ZenodoSettings(
    doi="10.5281/zenodo.4037939",
    filename="cc_params.npz",
    md5="0a449ddc4903ded610c4a88b9aac3514",
    url="https://zenodo.org/api/",
)

from .parameters import CCModelParameters
from .preprocessor import Preprocessor
