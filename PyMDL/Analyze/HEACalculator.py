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

import math
import itertools

import numpy as np
from thermo.elements import nested_formula_parser

from PyMDL.Data import Element, Mixing

__author__ = "Doguhan Sariturk"
__version__ = "1.2.0"
__email__ = "dogu.sariturk@gmail.com"
__status__ = "Development"


class HEACalculator:
    """General class for the high entropy alloys.

    Parameters
    ----------
    formula : str
        Alloy formula.

    Attributes
    ----------

    mixing_enthalpy: float
        Mixing enthalpy of the alloy. [1]_

    density: float
        Approximate density of the alloy.

    VEC: int
        Valence electron concentration of the alloy.

    melting_temperature: float
        Approximate melting temperature of the alloy.

    delta: float
        Atomic size difference of the alloy. [2]_

    mixing_entropy: float
        Mixing entropy of the alloy.

    omega: float
        Omega parameter. [3]_

    isSolidSolution: bool
        Return True if alloy forms solid solution.

    microstructure: str
        Expected crystal structure of the alloy. [4]_

    References
    ----------
    .. [1] Zhang, Y.; Zuo, T.T.; Tang, Z.; Gao, M.C.; Dahmen, K.A.; Liaw, P.K.; Lu, Z.P. Prog. Mater. Sci. 2014, 61.
    .. [2] S.S.Fang, X. S. Xiao, L. Xia, W. H. Li, Y. D. Dong,J. Non-Cryst. Solids 2003, 321, 120.
    .. [3] Yang, X.; Zhang, Y. Mater. Chem. Phys. 2012, 132, 233–238.
    .. [4] Guo, S.; Ng, C.; Lu, J.; Liu, C.T. J. Appl. Phys. 2011, 109, 103505.
    """

    _GAS_CONSTANT = 8.314462618

    def __init__(self, formula):
        self.formula = formula
        self._alloy = nested_formula_parser(self.formula)
        self._atomic_percentage = [num / sum(self._alloy.values())
                                   for num
                                   in self._alloy.values()]

        self.mixing_enthalpy = self.density = self.VEC = self.melting_temperature = self.delta = None
        self.mixing_entropy = self.omega = self.isSolidSolution = self.microstructure = None

    def calculate(self):
        """This method calculates phenomenological parameters based on thermodynamics and physics
        in order to predict the formation of solid solutions in High Entropy Alloys (HEAs).

        Raises
        ------
        TypeError
            If one of the entries are not in the database.
        """
        try:
            pair_list = list(itertools.combinations(self._alloy, 2))
            pair_enthalpy = [Mixing(pair)
                             for pair in pair_list]
            percentage = [(self._alloy[each[0]] / sum(self._alloy.values()))
                          * (self._alloy[each[1]] / sum(self._alloy.values()))
                          for each in pair_list]
            self.mixing_enthalpy = 4 * sum([percent * enthalpy
                                            for percent, enthalpy in zip(percentage, pair_enthalpy)])

            total_weight = sum(Element(elm).atomic_weight * af
                               for elm, af in self._alloy.items())
            total_volume = sum(Element(elm).atomic_volume * af
                               for elm, af in self._alloy.items())
            self.density = total_weight / total_volume

            nvalence_list = [Element(elm).nvalence
                             for elm in self._alloy.keys()]
            self.VEC = sum([percentage * valence
                            for percentage, valence in zip(self._atomic_percentage, nvalence_list)])

            melting_temperature_list = [Element(elm).melting_point
                                        for elm in self._alloy.keys()]
            self.melting_temperature = math.ceil(sum([percentage * melting_temp
                                                      for percentage, melting_temp
                                                      in zip(self._atomic_percentage, melting_temperature_list)]))

            atomic_radius_list = [Element(elm).atomic_radius
                                  for elm in self._alloy.keys()]
            average_atomic_radius = sum([percentage * radius
                                         for percentage, radius in zip(self._atomic_percentage, atomic_radius_list)])
            _delta = sum([percentage * (1 - (radius / average_atomic_radius)) ** 2
                          for percentage, radius in zip(self._atomic_percentage, atomic_radius_list)])
            self.delta = math.sqrt(_delta) * 100

            self.mixing_entropy = -1 * self._GAS_CONSTANT * sum(self._atomic_percentage * np.log(self._atomic_percentage))

            self.omega = (self.melting_temperature * self.mixing_entropy) / (abs(self.mixing_enthalpy) * 1000)

            self.isSolidSolution = True if self.omega >= 1.1 and 0 < self.delta < 6.6 and 5 > self.mixing_enthalpy > -15 else False

            if 2.5 <= self.VEC <= 3.5:
                self.microstructure = "HCP"
            elif self.VEC >= 8.0:
                self.microstructure = "FCC"
            elif self.VEC <= 6.87:
                self.microstructure = "BCC"
            else:
                self.microstructure = "BCC+FCC"

        except TypeError:
            print("TypeError: Formula contains elements which are not in the database!")

    def __repr__(self):
        try:
            return f"\n{self.formula}\n\n"\
               f"\tDensity: \t\t{self.density:.2f} g/cm^3\n"\
               f"\tDelta: \t\t\t{self.delta:.2f}\n"\
               f"\tMixing Enthalpy: \t{self.mixing_enthalpy:.2f} kJ/mol\n"\
               f"\tVEC: \t\t\t{self.VEC}\n"\
               f"\tMixing Entropy: \t{self.mixing_entropy:.2f} J/K.mol\n"\
               f"\tMelting Temperature: \t{self.melting_temperature} K\n"\
               f"\tOmega: \t\t\t{self.omega:.2f}\n"\
               f"\tCrystal Structure: \t{self.microstructure}\n"\
               f"\tIs Solid Solution: \t{self.isSolidSolution}\n"
        except AttributeError:
            return "AttributeError: Call calculate() method first."
