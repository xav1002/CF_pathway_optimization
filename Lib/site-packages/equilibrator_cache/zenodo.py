"""Handles downloading and caching of files from Zenodo."""

# The MIT License (MIT)
#
# Copyright (c) 2013 The Weizmann Institute of Science.
# Copyright (c) 2018 Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
# Copyright (c) 2018 Institute for Molecular Systems Biology,
# ETH Zurich, Switzerland.
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
import pathlib
from typing import NamedTuple

import appdirs
import pooch


logger = logging.getLogger(__name__)


class ZenodoSettings(NamedTuple):
    """
    Bundle the configuration for interacting with Zenodo.org.

    Attributes
    ----------
    doi : str
        The DOI of the equilibrator cache entry.
    filename : str
        The filename of the SQLite database.
    md5 : str
        The MD5 checksum of the SQLite database file.
    url : str
        The base URL of the API.

    """

    doi: str
    filename: str
    md5: str
    url: str


DEFAULT_COMPOUND_CACHE_SETTINGS = ZenodoSettings(
    doi="10.5281/zenodo.4128543",
    filename="compounds.sqlite",
    md5="9b66b85a886926d09755a66a3b452b6b",
    url="https://zenodo.org/api/",
)


def get_cached_filepath(settings: ZenodoSettings) -> pathlib.Path:
    """Get data from a file stored in Zenodo (or from cache, if available).

    Parameters
    ----------
    settings : ZenodoSettings
        Configuration for the interaction with Zenodo.org.

    Returns
    -------
    pathlib.Path
        The path to the locally cached file.

    """
    cache_directory = pathlib.Path(appdirs.user_cache_dir(appname="equilibrator"))
    cache_directory.mkdir(parents=True, exist_ok=True)

    cache_fname = pooch.retrieve(
        path=cache_directory,
        fname=settings.filename,
        url="doi:" + settings.doi + "/" + settings.filename,
        known_hash="md5:" + settings.md5,
        progressbar=True,
    )
    return pathlib.Path(cache_fname)
