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

import numpy as np
import pandas as pd

from PyMDL.Parsers.LAMMPS import RDF

__author__ = "Doguhan Sariturk"
__version__ = "0.1.0"
__email__ = "dogu.sariturk@gmail.com"
__status__ = "Development"
__maintainer__ = "Doguhan Sariturk"
__license__ = "GPL"


class StructureFactor:
    """" TODO

    TODO

    Parameters
    ----------
    rdf : pd.DataFrame or PyMDL.Parsers.LAMMPS.RDF object
        Radial Distribution Function. First column is the position.
        Number of columns depends on whether total or partial RDFs are provided.
        e.g., 2-columns (position, rdf) for total RDfs
              5-columns (position, rdf_11, rdf_12, rdf_21, rdf_22) for binary partial RDFs.

    number_density: float
        Number density in atoms/A^3

    Attributes
    ----------

    Returns
    -------

    """
    # TODO : Include a modifier function

    def __init__(self, rdf, number_density):
        super().__init__()

        if isinstance(rdf, RDF):
            self.rdf = rdf.rdf
        elif isinstance(rdf, pd.DataFrame):
            self.rdf = rdf

        self.number_density = number_density

        self.r = self.rdf.iloc[:, 0]

        self.q = self.r
        self.sq = self.c1 = self.c2 = None

    def _structureFactor(self):
        """ Calculate the total structure factor S(Q) from the pair correlation function g(r)
        using a naive Fourier Transformation.

        Returns
        -------
        s_q: pd.DataFrame
            A DataFrame of two or five columns. Number of columns depends on whether total or partial RDFs are provided.
            First column is the wavenumber, the remainders are Structure Factors.
        """
        dr = self.r.iloc[1] - self.r.iloc[0]
        multp = (4 * np.pi * self.number_density)

        if len(self.rdf.columns) == 2:
            self.gr = self.rdf.iloc[:, 1]
            sq = [np.sum(self.r * (self.gr - 1) * dr * np.sin(q * self.r)) / q
                  for q in self.r]
            sq = np.array(sq) * multp
            return 1 + sq

        elif len(self.rdf.columns) == 5:
            self.gr_1 = self.rdf.iloc[:, 1]
            self.gr_2 = self.rdf.iloc[:, 2]
            self.gr_3 = self.rdf.iloc[:, 3]
            self.gr_4 = self.rdf.iloc[:, 4]
            sq_1 = [np.sum(self.r * (self.gr_1 - 1) * dr * np.sin(q * self.r)) / q
                    for q in self.r]
            sq_1 = np.array(sq_1) * multp
            sq_2 = [np.sum(self.r * (self.gr_2 - 1) * dr * np.sin(q * self.r)) / q
                    for q in self.r]
            sq_2 = np.array(sq_2) * multp
            sq_3 = [np.sum(self.r * (self.gr_3 - 1) * dr * np.sin(q * self.r)) / q
                    for q in self.r]
            sq_3 = np.array(sq_3) * multp
            sq_4 = [np.sum(self.r * (self.gr_4 - 1) * dr * np.sin(q * self.r)) / q
                    for q in self.r]
            sq_4 = np.array(sq_4) * multp
            return 1 + np.vstack((sq_1, sq_2, sq_3, sq_4))

    @staticmethod
    def _scatteringFactor(row, el):
        """" TODO

        TODO

        Note
        ----

        Parameters
        ----------

        Attributes
        ----------

        Returns
        -------

        Raises
        ------

        References
        ----------

        """
        # A1 B1 A2 B2 A3 B3 A4 B4 C

        CONST = {"Al": [6.4202, 3.0387, 1.9002, 0.7426, 1.5936, 31.5472, 1.9646, 85.0886, 1.1151],
                 "Sm": [24.0042, 2.47274, 19.4258, 0.196451, 13.4396, 14.3996, 2.89604, 128.007, 2.20963]}

        f = CONST[el][0] * np.exp(-1 * CONST[el][1] * row["r"] ** 2) + CONST[el][8]
        + CONST[el][2] * np.exp(-1 * CONST[el][3] * row["r"] ** 2) + CONST[el][8]
        + CONST[el][4] * np.exp(-1 * CONST[el][5] * row["r"] ** 2) + CONST[el][8]
        + CONST[el][6] * np.exp(-1 * CONST[el][7] * row["r"] ** 2) + CONST[el][8]
        return f

    @staticmethod
    def _weightFactor_11(row, c1, c2):
        """" TODO

        TODO

        Note
        ----

        Parameters
        ----------

        Attributes
        ----------

        Returns
        -------

        Raises
        ------

        References
        ----------

        """
        return (c1 ** 2 * row["Al"] ** 2) / (c1 * row["Al"] + c2 * row["Sm"]) ** 2

    @staticmethod
    def _weightFactor_12(row, c1, c2):
        """" TODO

        TODO

        Note
        ----

        Parameters
        ----------

        Attributes
        ----------

        Returns
        -------

        Raises
        ------

        References
        ----------

        """
        return (2 * c1 * c2 * row["Al"] * row["Sm"]) / (c1 * row["Al"] + c2 * row["Sm"]) ** 2

    @staticmethod
    def _weightFactor_22(row, c1, c2):
        """" TODO

        TODO

        Note
        ----

        Parameters
        ----------

        Attributes
        ----------

        Returns
        -------

        Raises
        ------

        References
        ----------

        """
        return (c2 ** 2 * row["Sm"] ** 2) / (c1 * row["Al"] + c2 * row["Sm"]) ** 2

    def _faberZimanFormalism(self, sq, c1, c2):
        """Calculate the total structure factor S(Q) from the partial pair correlation functions g_ab(r)
        using the Faber-Ziman Formalism.

        ONLY FOR Al Sm binary alloys, for now.

        Using:
        :math:: S(q) = \\omega_{11}S_{11}(q) + \\omega_{12}S_{12}(q) + \\omega_{22}S_{22}(q)

        where:
        .. math::
            S_{ij}(q) = 1 + 4\\pi\\rho \\int_0^{\\infty} \\left[ g_{ij}(r) - 1 \\right] \\frac{\\sin(qr)}{r}rdr

        and
        .. math::
            \\omega_{11} &= \\frac{c_1^2f_1^2(q)}{\\left[ c_1f_1(q) + c_2f_2(q)\\right]^2} \\
            \\omega_{12} &= \\frac{2c_1c_2f_1(q)f_2(q)}{\\left[ c_1f_1(q) + c_2f_2(q)\\right]^2} \\
            \\omega_{22} &= \\frac{c_2^2f_2^2(q)}{\\left[ c_1f_1(q) + c_2f_2(q)\\right]^2}

        where
        .. math::
            f\\left(\\sin \\theta/\\lambda\\right) = \\sum_{i = 1}^4 a_i \\exp\\left(-b_i\\left(\\sin \\theta/\\lambda\\right)^2\\right) + c


        Parameters
        ----------
        sq : pd.DataFrame
            Partial structure factors of species obtained by Fourier Transform
            of the corresponding partial pair correlation functions.
        c1: float
            Composition (molar fraction) of specie 1
        c2: float
            Composition (molar fraction) of specie 2

        References
        ----------

        """
        r = pd.DataFrame(self.r)

        r["Al"] = r.apply(self._scatteringFactor, el="Al", axis=1)  # Get scattering factors for Al
        r["Sm"] = r.apply(self._scatteringFactor, el="Sm", axis=1)  # Get scattering factors for Sm

        r["w11"] = r.apply(self._weightFactor_11, c1=c1, c2=c2, axis=1)  # Get weight factor w_11
        r["w12"] = r.apply(self._weightFactor_12, c1=c1, c2=c2, axis=1)  # Get weight factor w_12
        r["w22"] = r.apply(self._weightFactor_22, c1=c1, c2=c2, axis=1)  # Get weight factor w_22

        r["Total"] = r["w11"] * sq[0] + r["w12"] * sq[1] + r["w22"] * sq[3]

        self.sq = r["Total"].to_numpy()

    def calculate_total_structure_factor(self, method="direct", c1=None, c2=None):
        """" TODO

        TODO

        Note
        ----

        Parameters
        ----------
        method : str
            Structure Factor calculation method, by default 'Direct'
            Available methods are 'Direct' and 'Faber-Ziman'
        c1 : float
            TODO
        c2 : float
            TODO
        Attributes
        ----------

        Returns
        -------

        Raises
        ------

        References
        ----------

        """

        if method.lower() == "direct":
            assert len(self.rdf.columns) == 2, "Total RDFs should be provided for total structure factor calculation " \
                                               "via 'Direct' method. "
            self.sq = self._structureFactor()
        elif method.lower() == "faber-ziman":
            assert len(self.rdf.columns) == 5, "Partial RDFs should be provided for total structure factor " \
                                               "calculation via 'Faber-Ziman' method. "
            if all(isinstance(i, (float, int))
                   for i in [c1, c2]):
                self.c1 = c1
                self.c2 = c2
            else:
                print(f"Error: Please provide compositions.")

            self.sq = self._structureFactor()
            self._faberZimanFormalism(self.sq, self.c1, self.c2)
        else:
            print(f"NotImplementedError: Given method, {method}, is not implemented.")

    def calculate_partial_structure_factor(self):
        """" TODO

        TODO

        Note
        ----

        Parameters
        ----------

        Attributes
        ----------

        Returns
        -------

        Raises
        ------

        References
        ----------

        """
        assert len(self.rdf.columns) == 5, "Partial RDFs should be provided for partial structure factor calculation."
        return self._structureFactor()

    def plot(self, ax=None):
        """" TODO

        TODO

        Note
        ----

        Parameters
        ----------

        Attributes
        ----------

        Returns
        -------

        Raises
        ------

        References
        ----------

        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        plt.style.use("grayscale")

        plt_params = {
                "text.usetex": True,
                "font.family": "serif",
                "axes.labelsize": 10,
                "font.size": 10,
                "axes.linewidth": 1,
                "legend.fontsize": 8,
                "xtick.labelsize": 8,
                "ytick.labelsize": 8,
        }
        mpl.rcParams.update(plt_params)

        if not ax:
            fig, ax = plt.subplots(1, 1, figsize=(6, 4))

        if self.sq.shape == (self.rdf.shape[0],):
            ax.plot(self.q, self.sq)
        elif self.sq.shape == (self.rdf.shape[0], 4):
            ax.plot(self.q, self.sq["sq_1"], label="sq-1")
            ax.plot(self.q, self.sq["sq_2"], label="sq-2")
            ax.plot(self.q, self.sq["sq_3"], label="sq-3")
            ax.plot(self.q, self.sq["sq_4"], label="sq-4")

            plt.legend()

        ax.set_xlabel(r"q(Å$^{-1}$)", labelpad=10)
        ax.set_ylabel("S(q)", labelpad=10)

        plt.tight_layout()
        plt.show()

        return plt
