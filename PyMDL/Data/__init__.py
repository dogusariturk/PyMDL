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

from PyMDL.Data.Elements import _Element, _element_data
from PyMDL.Data.Mixing import Mixing

__author__ = 'Doguhan Sariturk'
__version__ = "0.1.0"
__email__ = 'dogu.sariturk@gmail.com'
__status__ = 'Development'
__maintainer__ = 'Doguhan Sariturk'
__license__ = "GPL"


def Element(name=None):
    """TODO

    Parameters
    ----------
    name : [type], optional
        [description], by default None

    Returns
    -------
    [type]
        [description]
    """
    if name is None:
        print("Usage: Element('ELEMENT_NAME')")
    elif name in _element_data:
        return _Element(name, _element_data[name]['melting_point'], _element_data[name]['atomic_number'],
                        _element_data[name]['atomic_volume'], _element_data[name]['atomic_weight'],
                        _element_data[name]['atomic_radius'], _element_data[name]['nvalence'])
    else:
        print('Enter a valid element name.')
