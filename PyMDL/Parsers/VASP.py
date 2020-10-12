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


class OSZICAR:
    """A class for a VASP OSZICAR file.

    Information on each electronic and ionic SCF step.

    Parameters
    ----------
    filename : str
        The location of the OSZICAR file.
    """
    def __init__(self, filename):
        super().__init__()
        self._filename = filename

    def parse(self):
        pass


class OUTCAR:
    """A class for a VASP OUTCAR file.

    Main output file.

    Parameters
    ----------
    filename : str
        The location of the OUTCAR file.
    """
    def __init__(self, filename):
        super().__init__()
        self._filename = filename

    def parse(self):
        pass


class PCDAT:
    """A class for a VASP PCDAT file.

    Contains the pair correlation function.

    Parameters
    ----------
    filename : str
        The location of the PCDAT file.
    """
    def __init__(self, filename):
        super().__init__()
        self._filename = filename

    def parse(self):
        pass


class XDATCAR:
    """A class for a VASP XDATCAR file.

    Contains ionic configuration for each output step of molecular dynamics simulations.

    Parameters
    ----------
    filename : str
        The location of the XDATCAR file.
    """
    def __init__(self, filename):
        super().__init__()
        self._filename = filename

    def parse(self):
        pass


class CONTCAR:
    """A class for a VASP CONTCAR file.

    Is the updated POSCAR file after each calculation, whether ionic movement was performed or not.

    Parameters
    ----------
    filename : str
        The location of the XDATCAR file.
    """
    def __init__(self, filename):
        super().__init__()
        self._filename = filename

    def parse(self):
        pass