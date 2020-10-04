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

__author__ = 'Doguhan Sariturk'
__version__ = "0.1.0"
__email__ = 'dogu.sariturk@gmail.com'
__status__ = 'Development'
__maintainer__ = 'Doguhan Sariturk'
__license__ = "GPL"


# TODO : refactor exceptions


class SnapshotDump:
    """Class for representing a snapshot from LAMMPS Dump.


    Parameters
    ----------
    natoms: int
        Number of atoms in the system.

    timestep: int(s)
        The timesteps at which the snapshots created.

    box: numpy.ndarray(s), [[xlo, xhi],[ylo, yhi], [zlo, zhi]]
        Simulation box boundaries.

    data: numpy.ndarray(s), [[id, x, y, z], ]
        IDs and coordinates of all atoms in the system.
    """

    def __init__(self, box, data, natoms, timestep):
        self.box = box
        self.data = data
        self.natoms = natoms
        self.timestep = timestep

    def __repr__(self):
        return f"Dump.get_snapshot({self.timestep})"


class Dump:
    """A class for a LAMMPS dump file in 'atom' style.


    Parameters
    ----------
    filename : str
        The location of the dump file.

    Attributes
    ----------
    natoms: int
        Number of atoms in the system.

    nsnapshot: int
        The number of snapshots the dump file has.

    timestep: int(s)
        The timesteps at which the snapshots created.

    box: numpy.ndarray(s), [[xlo, xhi, ylo, yhi, zlo, zhi], ]
        Simulation box boundaries.

    data: numpy.ndarray(s), [[id, x, y, z], ]
        IDs and coordinates of all atoms in the system.

    Returns
    -------
    None
    """

    def __init__(self, filename):
        # TODO : Implement element names.
        # TODO : combine _get_natoms() and _get_nsnapshot() into _parse_header()
        self._filename = filename

        self._timestep = self.box = self.data = None

        self._get_natoms()
        self._get_nsnapshot()

    def _get_natoms(self):
        """Retrieves the number of atoms in the system.

        Returns
        -------
        None

        Raises
        ------
        FileNotFoundError
            If 'filename' does not exist.
        """
        try:
            with open(self._filename, 'rt') as f:
                f.readline()  # 'ITEM: TIMESTEP\n'
                f.readline()
                f.readline()  # 'ITEM: NUMBER OF ATOMS\n'
                self.natoms = int(f.readline())
        except FileNotFoundError as e:
            print(f"FileNotFoundError: \n \t {e.strerror}: '{e.filename}'")

    def _get_nsnapshot(self):
        """Retrieves the number of snapshots the dump file has.

        Returns
        -------
        None

        Raises
        ------
        FileNotFoundError
            If 'filename' does not exist.
        """
        try:
            with open(self._filename, 'rt') as f:
                d = f.read()
                self.nsnapshot = d.count('ITEM: TIMESTEP\n')
        except FileNotFoundError as e:
            print(f"FileNotFoundError: \n \t {e.strerror}: '{e.filename}'")

    def parse_one(self):
        """Parses the first snapshot of a LAMMPS dump file in 'atom' style.

        Returns
        -------
        SnapshotDump

        Raises
        ------
        FileNotFoundError
            If 'filename' does not exist.
        """
        try:
            with open(self._filename, 'rt') as f:
                f.readline()  # 'ITEM: TIMESTEP\n'
                timestep = int(f.readline())
                f.readline()  # 'ITEM: NUMBER OF ATOMS\n'
                natoms = int(f.readline())
                f.readline()
                xlo, xhi = [float(s) for s in f.readline().split()]
                ylo, yhi = [float(s) for s in f.readline().split()]
                zlo, zhi = [float(s) for s in f.readline().split()]

                box = np.array([xlo, xhi, ylo, yhi, zlo, zhi])

                f.readline()  # 'ITEM: ATOMS id type xs ys zs\n'

                data = np.zeros((natoms, 5), np.float)

                for i in range(natoms):
                    data[i] = f.readline().split()

            return SnapshotDump(box, data, natoms, timestep)

        except FileNotFoundError as e:
            print(f"FileNotFoundError: \n \t {e.strerror}: '{e.filename}'")

    def parse_all(self):
        """Parses all snapshots of a LAMMPS dump file in 'atom' style.

        Returns
        -------
        None

        Raises
        ------
        FileNotFoundError
            If 'filename' does not exist.
        """
        try:
            with open(self._filename, 'rt') as f:

                self._timestep = np.array([], dtype='int')
                __box = list()
                __data = list()

                for line in f:
                    if line.startswith('ITEM: TIMESTEP'):
                        __timestep = int(next(f).split()[0])
                        self._timestep = np.append(self._timestep, __timestep)

                    if line.startswith('ITEM: BOX'):
                        xlo, xhi = [float(s) for s in next(f).split()]
                        ylo, yhi = [float(s) for s in next(f).split()]
                        zlo, zhi = [float(s) for s in next(f).split()]

                        __box.append([xlo, xhi, ylo, yhi, zlo, zhi])

                    if line.startswith('ITEM: ATOMS'):
                        __data_snapshot = list()
                        for _ in range(int(self.natoms)):
                            item = list(next(f).split())
                            __data_snapshot.append(item)
                        __data.append(__data_snapshot)

                self.box = np.array(__box)
                self.data = np.array(__data, dtype='float64')

            _string = ' '.join(map(str, self._timestep))
            print(f"Success:\n\tParsed {self.nsnapshot} snapshots, at timesteps {_string}.")
            print("\tNow, use get_snapshot() method to get a SnapshotDump object at that timestep.")

        except FileNotFoundError as e:
            print(f"FileNotFoundError: \n \t {e.strerror}: '{e.filename}'")

    def get_snapshot(self, requested_timestep):
        """Gets the snapshot at a given timestep, if it exists.

        Returns
        -------
        SnapshotDump object

        Raises
        ------
        AttributeError
            If no parsing was done.

        ValueError
            If requested timestep does not exist.
        """
        try:
            if requested_timestep not in self._timestep:
                print("ValueError: \n \t Requested timestep does not exist.")
            else:
                _index = int(np.where(self._timestep == requested_timestep)[0])
                return SnapshotDump(self.box[_index], self.data[_index], self.natoms, self._timestep[_index])
        except AttributeError:
            print(
                    "Error: \n \t Use parse_all() before calling get_snapshot(requested_timestep).")


class SnapshotLog:
    """Class for representing a single LAMMPS Log as a snapshot.


    Parameters
    ----------
    header: str
        List of vector names.

    nrun: int
        Simulation run time.

    data: numpy.ndarray(s), [[id, x, y, z], ]
        IDs and coordinates of all atoms in the system.

    Returns
    -------
    None
    """

    def __init__(self, box, data, natoms, timestep):
        self.box = box
        self.data = data
        self.natoms = natoms
        self.timestep = timestep

    def __repr__(self):
        # TODO : change repr string
        return "This method returns a SnapshotLog class and should be instantiated."


class Log:
    """A class for a LAMMPS log file.

    Parameters
    ----------
    filename : str
        The location of the dump file.

    Returns
    -------
    None
    """

    def __init__(self, filename):
        self._filename = filename
        self._nrun = []

    def _parse_header(self):
        """Parses the LAMMPS log file header.

        Returns
        -------
        None

        Raises
        ------
        FileNotFoundError
            If 'filename' does not exist.
        """
        try:
            with open(self._filename, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    if line.startswith('run'):
                        self._nrun.append(line.split()[1])
        except FileNotFoundError as e:
            print(f"FileNotFoundError: \n \t {e.strerror}: '{e.filename}'")


class SnapshotRDF:
    """Class for representing a single LAMMPS RDF as a snapshot.


    Parameters
    ----------
    rdf : Pandas.DataFrame
        TODO
    n_column : int
        TODO

    Returns
    -------
    None
    """

    def __init__(self, rdf, n_column):
        self.rdf = rdf
        self.n_column = n_column

    def plot(self):
        """" TODO

        Returns
        -------
        None

        """
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        plt.style.use('grayscale')

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

        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
        if self.n_column == 4:
            ax.plot(self.rdf['r'], self.rdf['gr'])
        elif self.n_column == 10:
            ax.plot(self.rdf['r'], self.rdf['gr_1'], label='gr-1')
            ax.plot(self.rdf['r'], self.rdf['gr_2'], label='gr-2')
            ax.plot(self.rdf['r'], self.rdf['gr_3'], label='gr-3')
            ax.plot(self.rdf['r'], self.rdf['gr_4'], label='gr-4')

        plt.xlabel("r(Å)", labelpad=10)
        plt.ylabel("g(r)", labelpad=10)

        plt.legend()

        plt.show()

        return ax


class RDF:
    """A class for a LAMMPS RDF file to obtain total RDF.

    LAMMPS RDF file should be created using fix ave/time.

    Parameters
    ----------
    filename : str
        The location of the RDF file.

    Attributes
    ----------
    n_column : int
        TODO
    n_bin : int
        TODO
    num_rdf : int
        TODO

    Returns
    -------
    None
    """

    def __init__(self, filename):
        self._filename = filename

        self._parse_header()

    def _parse_header(self):
        """Parses the LAMMPS rdf file header.

        Returns
        -------
        None

        Raises
        ------
        FileNotFoundError
            If 'filename' does not exist.
        """
        try:
            self._file = pd.read_csv(self._filename, skiprows=[0, 1], sep=' ')

            self.n_column = int(len(self._file.axes[1]) - 1)
            self.n_bin = int(self._file.iloc[0][1])
            self.num_rdf = int((len(self._file)) / (self.n_bin + 1))
            print(f"Success:\n\tRDF file contains {self.num_rdf} snapshots.")
        except FileNotFoundError as e:
            print(f"FileNotFoundError: \n \t {e.strerror}: '{e.filename}'")

    def _smooth_rdf(self, rdf):
        """TODO

        Parameters
        ----------
        rdf : pandas.DataFrame
            TODO

        Returns
        -------
        pandas.DataFrame
            TODO
        """
        from scipy.interpolate import make_interp_spline
        new_rdf = pd.DataFrame()
        r_new = np.linspace(rdf['r'].min(), rdf['r'].max(), 500)
        if self.n_column == 4:
            spl = make_interp_spline(rdf['r'], rdf['gr'], k=3)
            gr_smooth = spl(r_new)
            new_rdf['r'] = r_new
            new_rdf['gr'] = gr_smooth
        elif self.n_column == 10:
            spl_1 = make_interp_spline(rdf['r'], rdf['gr_1'], k=3)
            spl_2 = make_interp_spline(rdf['r'], rdf['gr_2'], k=3)
            spl_3 = make_interp_spline(rdf['r'], rdf['gr_3'], k=3)
            spl_4 = make_interp_spline(rdf['r'], rdf['gr_4'], k=3)
            gr_1_smooth = spl_1(r_new)
            gr_2_smooth = spl_2(r_new)
            gr_3_smooth = spl_3(r_new)
            gr_4_smooth = spl_4(r_new)
            new_rdf = pd.DataFrame()
            new_rdf['r'] = r_new
            new_rdf['gr_1'] = gr_1_smooth
            new_rdf['gr_2'] = gr_2_smooth
            new_rdf['gr_3'] = gr_3_smooth
            new_rdf['gr_4'] = gr_4_smooth

        return new_rdf

    def get_snapshot(self, required_snapshot, smooth=False):
        """TODO

        Parameters
        ----------
        required_snapshot : pandas.DataFrame
            TODO

        smooth : bool
            TODO, by default False

        Returns
        -------
        pandas.DataFrame
            TODO
        """
        # TODO : check if required_snapshot >= 0
        start = required_snapshot
        stop = (required_snapshot + 1)
        return self.average_snapshots(start=start, stop=stop, smooth=smooth)

    def average_snapshots(self, start=0, stop=1, step=1, smooth=False):
        """TODO

        Parameters
        ----------
        start : int
            TODO, by default 0

        stop : int
            TODO , by default 1

        step : int
            TODO , by default 1

        smooth : bool
            TODO, by default False

        Returns
        -------
        pandas.DataFrame
            TODO
        """

        # TODO : check if start >= 0
        # TODO: check if stop <= self.num_rdf

        rdf = pd.DataFrame()
        for i in range(start, stop, step):
            df = self._file.iloc[(self.n_bin * i) + i + 1: (self.n_bin * i) + i + (self.n_bin + 1), :]
            if i == start:
                rdf['r'] = df.iloc[:, 1].values
                if self.n_column == 4:
                    rdf['gr'] = df.iloc[:, 2].values
                elif self.n_column == 10:
                    rdf['gr_1'] = df.iloc[:, 2].values
                    rdf['gr_2'] = df.iloc[:, 4].values
                    rdf['gr_3'] = df.iloc[:, 6].values
                    rdf['gr_4'] = df.iloc[:, 8].values
            else:
                if self.n_column == 4:
                    rdf['gr'] = rdf['gr'] + df.iloc[:, 2].values
                elif self.n_column == 10:
                    rdf['gr_1'] = rdf['gr_1'] + df.iloc[:, 2].values
                    rdf['gr_2'] = rdf['gr_2'] + df.iloc[:, 4].values
                    rdf['gr_3'] = rdf['gr_3'] + df.iloc[:, 6].values
                    rdf['gr_4'] = rdf['gr_4'] + df.iloc[:, 8].values

        if self.n_column == 4:
            rdf['gr'] = rdf['gr'] / ((stop - start) / step)
        elif self.n_column == 10:
            rdf['gr_1'] = rdf['gr_1'] / ((stop - start) / step)
            rdf['gr_2'] = rdf['gr_2'] / ((stop - start) / step)
            rdf['gr_3'] = rdf['gr_3'] / ((stop - start) / step)
            rdf['gr_4'] = rdf['gr_4'] / ((stop - start) / step)

        if not smooth:
            return SnapshotRDF(rdf, self.n_column)
        elif smooth:
            new_rdf = self._smooth_rdf(rdf)
            return SnapshotRDF(new_rdf, self.n_column)


class Monitor:
    """A class for a LAMMPS Monitor file.

    LAMMPS Monitor file should be created using fix print.

    Parameters
    ----------
    filename : str
        The location of the Monitor file.

    string: str
        Text string used to print the Monitor file.
    """

    def __init__(self, filename, string):
        super().__init__()

        self.filename = filename
        self.string = string
