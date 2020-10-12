#!/usr/bin/env python3

# -*- coding: utf-8 -*-

# Copyright (C) 2020  Doguhan Sariturk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__author__ = "Doguhan Sariturk"
__version__ = "0.1.0"
__email__ = "dogu.sariturk@gmail.com"
__status__ = "Development"
__maintainer__ = "Doguhan Sariturk"
__license__ = "GPL"

from itertools import combinations_with_replacement

import numpy as np

from PyMDL.Parsers.LAMMPS import Dump


class RDF:
    """A class for calculating partial radial distribution functions from LAMMPS Dump files.

    Note
    -----
    TODO: RDF algorithm from mpmorph, include LICENSE

    Parameters
    ----------
    dump : Dump
        The object name of the dump file.

    elements : list
        A list of element names correspond to LAMMPS atom types in the dump file.

    cutoff : float, optional
        Maximum distance for which pair RDFs will be calculated. Defaults to 8.0 Angstroms.

    n_bin : int, optional
        Number of RDF bins

    smooth : bool, optional
        Whether recursively smooth RDFs using Savitzky-Golay filter. Defaults to True.

    Attributes
    ----------
    pairs : list
        A list of 2-length tuples of element pairs.

    natoms: int
        Number of atoms in the system.

    nsnapshot: int
        The number of snapshots the dump file has.

    Returns
    -------
    RDF : dict
        A dictionary with element pairs as keys and RDFs as values.


    """
    def __init__(self, dump, elements, cutoff=8.0, n_bin=100, smooth=True):
        super().__init__()

        self.elements = list(elements)
        self.cutoff = cutoff
        self.n_bin = n_bin
        self.pairs = [pair for pair in combinations_with_replacement(self.elements, 2)]
        self.smooth = smooth

        self.natoms = dump.natoms
        self.nsnapshot = dump.nsnapshot
        data = dump.data

        self.atom_types = data[:, 1]
        self.coords = data[:, 2:]

        self.distance_matrix = self._distance_matrix(self.coords)

        self.RDF = {}

    def _distance_matrix(self, coords):
        """"Returns the distance matrix between all sites in the structure.

        Note
        -----
        TODO: Distance Matrix from pymatgen, include LICENSE
        Parameters
        ----------

        Returns
        -------

        """
        z = (coords[:, None, :] - coords[None, :, :]) ** 2
        return np.sum(z, axis=-1) ** 0.5

    def calculate(self):
        pass
