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

import numpy as _np

__author__ = 'Doguhan Sariturk'
__version__ = "0.1.0"
__email__ = 'dogu.sariturk@gmail.com'
__status__ = 'Development'
__maintainer__ = 'Doguhan Sariturk'
__license__ = "GPL"


class _SnapshotDump:
    def __init__(self, box, data, natoms, timestep):
        """Class for representing a single LAMMPS Dump as a snapshot.


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
        self.box = box
        self.data = data
        self.natoms = natoms
        self.timestep = timestep

    def __repr__(self):
        return "This method returns a SnapshotDump class and should be instantiated."


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
    """

    def __init__(self, filename):
        self._filename = filename
        self._get_natoms()
        self._get_nsnapshot()

    def _get_natoms(self):
        """Retrieves the number of atoms in the system.

        Returns
        -------
        None

        Raises
        ------
        IOError
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
        IOError
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
        SnapshotDump object

        Raises
        ------
        IOError
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

                box = _np.array([xlo, xhi, ylo, yhi, zlo, zhi])

                f.readline()  # 'ITEM: ATOMS id type xs ys zs\n'

                data = _np.zeros((natoms, 5), _np.float)

                for i in range(natoms):
                    data[i] = f.readline().split()

            return _SnapshotDump(box, data, natoms, timestep)

        except FileNotFoundError as e:
            print(f"FileNotFoundError: \n \t {e.strerror}: '{e.filename}'")

    def parse_all(self):
        """Parses all snapshots of a LAMMPS dump file in 'atom' style.

        Returns
        -------
        None

        Raises
        ------
        IOError
            If 'filename' does not exist.
        """
        try:
            with open(self._filename, 'rt') as f:

                self._timestep = _np.array([], dtype='int')
                __box = list()
                __data = list()

                for line in f:
                    if line.startswith('ITEM: TIMESTEP'):
                        __timestep = int(next(f).split()[0])
                        self._timestep = _np.append(self._timestep, __timestep)

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

                self.box = _np.array(__box)
                self.data = _np.array(__data)

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
                _index = int(_np.where(self._timestep == requested_timestep)[0])
                return _SnapshotDump(self.box[_index], self.data[_index], self.natoms, self._timestep[_index])
        except AttributeError:
            print(
                "Error: \n \t Use parse_all() before calling get_snapshot(requested_timestep).")


class _SnapshotLog:
    def __init__(self, box, data, natoms, timestep):
        """Class for representing a single LAMMPS Log as a snapshot.


        Parameters
        ----------
        header: str
            List of vector names.

        nrun: int
            Simulation run time.

        data: numpy.ndarray(s), [[id, x, y, z], ]
            IDs and coordinates of all atoms in the system.
        """
        self.box = box
        self.data = data
        self.natoms = natoms
        self.timestep = timestep

    def __repr__(self):
        return "This method returns a SnapshotDump class and should be instantiated."


class Log:
    """A class for a LAMMPS log file.


    Parameters
    ----------
    filename : str
        The location of the dump file.
    """

    def __init__(self, filename):
        self._filename = filename
        self._nrun = []

    def _parse_header(self):
        """Parses the LAMMPS log file.

        Returns
        -------
        None


        Raises
        ------
        IOError
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


class RDF:
    """A class for a LAMMPS RDF file.

    LAMMPS RDF file should be created using fix ave/time.

    Parameters
    ----------
    filename : str
        The location of the RDF file.
    """
    def __init__(self):
        super().__init__()


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
    def __init__(self):
        super().__init__()


if __name__ == "__main__":
    print("")
    print("LAMMPS Parser")
    print("\nUsage:")
    print("\t - Import as 'from PyMDL.Parsers import Dump' from the base dir and initialize an object with Dump('DUMP_FILE')")
    print("\t   Then, parse the DUMP_FILE via obj.parse_one() or obj.parse_all().\n")
    print("\t - If obj.parse_all() was used, get the data at desired timestep via obj.get_snapshot(requested_timestep).")
