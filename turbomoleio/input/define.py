# -*- coding: utf-8 -*-
# The turbomoleio package, a python interface to Turbomole
# for preparing inputs, parsing outputs and other related tools.
#
# Copyright (C) 2018-2021 BASF SE, Matgenix SRL.
#
# This file is part of turbomoleio.
#
# Turbomoleio is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Turbomoleio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with turbomoleio (see ~turbomoleio/COPYING). If not,
# see <https://www.gnu.org/licenses/>.

"""
Module containing the required classes to run define based on a set of input parameters.
"""

import os
import shutil
import re
import signal
import pexpect
from pexpect import popen_spawn
from copy import deepcopy
from monty.os import cd
from turbomoleio.core.control import sdg, cdg, mdgo, Control


class DefineError(Exception):
    """Base exception for define errors"""

    def __init__(self, message=None):
        """
        Args:
            message (str): the message describing the error
        """
        super().__init__(message)
        self.message = message


class DefineExpectError(DefineError):
    """
    Exception raised when an expect fail, going in timeout or EOF
    """

    def __init__(self, message, pattern):
        """

        Args:
            message (str): the message describing the error
            pattern (list): the pattern passed to define that failed to be found
        """
        super().__init__(message=message)
        self.pattern = pattern


class DefineParameterError(DefineError):
    """
    Exception raised when parameters are not defined properly, e.g. incompatible among them or
    with the other files provided, like coord or a previous control.
    """
    pass


class DefineIredError(DefineError):
    """
    Exception raised for problems related to the generation of the internal coordinates.
    Even in the case of incomplete set in internal coordinates and ired is True.
    """
    pass


class UnsupportedDefineStateError(DefineError):
    """
    Exception raised when according to the current implementation a specific behavior is
    expected from define, but a proper way of handling with the case has not been implemented.
    """
    pass


class DefineRunner:
    """
    Class that runs the define executable using pexpect to communicate with the code.
    """

    def __init__(self, parameters, timeout=60, log_filepath="define.log", workdir=None,
                 executable=None, use_popen=True):
        """
        Args:
            parameters (dict): dictionary with the parameters defining the execution. The dictionary will be
                copied to avoid side effects during the execution.
            timeout (int): integer with the maximum time waited by "expect" before raising a timeout exception.
            log_filepath (str): string with a path to the file used to log the input and output of define.
                If None no log will be created.
            workdir (str): directory where the define will be executed. If None runs in the current directory.
            executable (str): path to the executable that will be used to run define. If None will first
                search in an environment variable and then default to "define", that should be in the $PATH
            use_popen (bool): if True the pexpect.popen_spawn.PopenSpawn implementation of pexpect based
                on popen is used. Otherwise the standard spawn (based on pty).
        """

        self.parameters = deepcopy(parameters)
        self.timeout = timeout
        self.log_filepath = log_filepath
        workdir = workdir or os.getcwd()
        self.workdir = os.path.abspath(workdir)
        self.executable = executable
        self.history = []
        self.sent_commands = []
        self.define = None
        self.use_popen = use_popen
        # set the internal value with the spawn command that should be used based when
        # starting pexpect
        if self.use_popen:
            self._spawn = pexpect.popen_spawn.PopenSpawn
        else:
            self._spawn = pexpect.spawn

    def _expect(self, pattern, timeout=-1, action=None):
        """
        Runs self.define.expect with the options provided.
        Handles the timeout and eof exceptions, logs and raises them as specific define exceptions.

        Args:
            pattern: the pattern given to pexpect.expect(). See pexpect documentations for more details.
            timeout (int): number of second that should be waited before raising the TIMEOUT exception.
                The default (-1) keeps the timeout value set in spawn.
            action (str): a string describing the action that is being performed. Will be used for logging and
                tracing of errors.

        Returns:
            int: The index of the options provided, as returned by pexpect.expect().
        """

        log_msg = "Expect: {}".format(pattern)
        if action:
            log_msg += ". Action: {}".format(action)
        self.history.append(log_msg)

        try:
            return self.define.expect(pattern, timeout=timeout)
        except (pexpect.exceptions.TIMEOUT, pexpect.exceptions.EOF) as e:
            msg = "Expect could not find the pattern: {}.".format(pattern)

            if self._is_alive:
                msg += " The job was alive."
            else:
                msg += " The job was not alive."
            raise DefineExpectError(message=msg, pattern=pattern)

    def _sendline(self, line, action=""):
        """
        Calls self.define.sendline.
        Logs the action and the list of commands sent.

        Args:
            line (str): a string with the line that should be sent to define.
            action (str): a string describing the action that is being performed. Will be used for logging and
                tracing of errors.

        Returns:
            None
        """

        log_msg = "Sendline: {}".format(line)
        if action:
            log_msg += ". Action: {}".format(action)
        self.history.append(log_msg)
        self.sent_commands.append(line)

        self.define.sendline(line)

    @property
    def _is_alive(self):
        """
        Utility function to determine if the process is alive. Calls the correct method
        depending on the pexpect implementation used (popen or pty)

        Returns:
            bool: True is the process is alive
        """
        if self.use_popen:
            return self.define.proc.poll() is None
        else:
            return self.define.isalive()

    def _close(self):
        """
        Closes the connection with the child application and/or kills the child process
        depending on the pexpect implementation used (popen or pty)

        Returns:
            None
        """
        if self.use_popen:
            if self._is_alive:
                self.define.kill(sig=signal.SIGKILL)
        else:
            self.define.close(force=True)

    def _get_bin_path(self):
        """
        Gives the path to the define executable based on the parameters.

        Returns:
            str: The path to the define executable
        """
        bin_path = self.executable

        if not bin_path:
            bin_path = "define"

        return bin_path

    def run_full(self):
        """
        Runs define going through all the menus using the parameters. If The jobs reaches the end
        without exceptions being raised it will check if the job ended normally or abnormally.

        Returns:
            bool: True if define ended normally otherwise abnormally.
        """
        self._set_metric()

        ended_normally = False

        with cd(self.workdir):
            with open(self.log_filepath, "wb") as logfile:
                self.define = self._spawn(self._get_bin_path(), timeout=self.timeout, logfile=logfile)
                try:
                    self._initialize_control(update=False)

                    # coordinate menu
                    self._geometry_menu(new_coords=True)

                    # leave the previous menu and move to the next one
                    self._switch_to_atomic_attribute_menu()

                    # atomic attribute menu - only basis set is used here
                    self._define_basis_sets()

                    # leave the previous menu and move to the next one
                    self._switch_to_molecular_orbital_definition_menu()

                    # molecular orbital definition menu
                    if self.parameters.get("usemo", None) is not None:
                        self._set_mo_from_control_file()
                    else:
                        self._extended_hueckel_theory_menu()

                    # excited state menu
                    ex_keys = {"ex_method", "ex_multi", "ex_irrep_states", "ex_all_states", "ex_frequency", "ex_exopt"}
                    method = self.parameters["method"].lower()
                    if any(self.parameters.get(k, None) for k in ex_keys) and method != "adc(2)":
                        self._excited_state_menu()

                    # general menu
                    self._general_menu()

                    self._quit_general_menu()

                    case = self._expect(["define ended normally", "define ended abnormally"],
                                        action="check end of define")

                    if case == 0:
                        ended_normally = True

                    # wait for the EOF with a timeout
                    self._expect(pexpect.EOF, timeout=20, action="wait for EOF")

                    # this is a blocking call. It may be removed
                    self.define.wait()

                finally:
                    # Always close define
                    self._close()

            if self.parameters.get("copymo", None) is not None:
                self._copy_mo_files()

            self._postprocess()

        return ended_normally

    def run_update_internal_coords(self):
        """
        Runs the define in a folder where a control file is already present and uses the parameters to
        update the internal coordinates of the system. Quits afterwards terminating define "abnormally".

        Returns:
            bool: True if internal coordinates were successfully updated.
        """

        coords_updated = False
        self._set_metric()

        with cd(self.workdir):
            with open(self.log_filepath, "wb") as logfile:
                self.define = self._spawn(self._get_bin_path(), timeout=self.timeout, logfile=logfile)
                try:
                    self._initialize_control(update=True)

                    # coordinate menu
                    self._geometry_menu(new_coords=False, update=True)

                    self._sendline("*", action="leave geometry menu")

                    case = self._expect(["GEOMETRY DATA WILL BE WRITTEN TO FILE",
                                         "IF YOU DO NOT WANT TO USE INTERNAL COORDINATES ENTER"])

                    if case == 0:
                        coords_updated = True
                    else:
                        raise DefineIredError("An incomplete set of internal coordinate is present with ired True")

                    self._sendline("qq", "quit after geometry modification")

                    case = self._expect(["define ended abnormally"], action="check end of define")

                    # wait for the EOF with a timeout
                    self._expect(pexpect.EOF, timeout=20, action="wait for EOF")

                    # this is a blocking call. It may be removed
                    self.define.wait()

                finally:
                    # Always close define
                    self._close()

        return coords_updated

    def run_generate_mo_files(self):
        """
        Runs define in a folder where a control file is already present and only regenerates
        the molecular orbital files. If the usemo option is defined the related molecular orbital
        will be taken from there and updated, otherwise they will be generated from scratch with
        the extended hueckel guess. (N.B. the control should not be the one currently being edited
        by define)

        Also allows to change some information in the geometry. If at least one of the following
        options is enabled it will enter the geometry menu: sym, desy, ired. Notice that if the
        geometry menu should be skipped these options should be disabled.

        Quits terminating define "abnormally".

        The function can only return True or fail with an exception. The return value is given
        to provide an interface similar to the other running methods.

        Returns:
            bool: always returns True if the function completes.
        """

        with cd(self.workdir):
            with open(self.log_filepath, "wb") as logfile:
                self.define = self._spawn(self._get_bin_path(), timeout=self.timeout, logfile=logfile)
                try:
                    self._initialize_control(update=True)

                    # False if all the options are False, None or empty strings
                    update_geom = any(self.parameters.get(opt) for opt in ["sym", "desy", "ired"])

                    self._geometry_menu(new_coords=False, update=update_geom)

                    self._skip_atomic_attribute_menu(update_geom)

                    case = self._expect([".*MOLECULAR ORBITAL DEFINITION.*",
                                         "MOLECULAR ORBITAL DATA.*DO YOU WANT TO CHANGE THESE DATA"],
                                        action="switch to molecular orbital definition menu")

                    if case == 1:
                        self._sendline("y", action="start modification of mo data")

                    # molecular orbital definition menu
                    if self.parameters.get("usemo", None) is not None:
                        self._set_mo_from_control_file()
                    else:
                        self._extended_hueckel_theory_menu()

                    # quit here. At this point some questions may still be asked (like if some groups
                    # should be deleted or similar). They do not seem to affect the generation of the
                    # mos, that has already been completed, and will add more complication without
                    # any particular benefit.
                    self._sendline("qq", "quit after geometry modification")

                    case = self._expect(["define ended abnormally"], action="check end of define")

                    # wait for the EOF with a timeout
                    self._expect(pexpect.EOF, timeout=20, action="wait for EOF")

                    # this is a blocking call. It may be removed
                    self.define.wait()

                finally:
                    # Always close define
                    self._close()

        return True

    def _set_metric(self):
        if self.parameters.get("metric", None) is not None :
            if not os.path.isfile("control"):
                # create the control file with redund_inp/metric entry
                c = Control.from_metric(self.parameters["metric"])
                c.to_file("control")
            else:
                mdgo("$redund_inp", {"metric": "  metric {}".format(self.parameters["metric"])})

    def _initialize_control(self, update=False):
        """
        Initialize the control file in define checking if a control file is already existing in the folder
        and setting the title. If update is True a previous control file should be present.

        Args:
            update (bool): determines if the current run is supposed to update a fully defined control.

        Returns:
            None
        """

        case = self._expect(["DATA WILL BE TAKEN FROM control BY DEFAULT",
                             "DATA WILL BE WRITTEN TO THE NEW FILE control"],
                            action="start inizialize control")

        # appending data to fragmentary control file
        if case == 1:
            if update:
                raise DefineError("No control file in the folder {}. Required for update".format(os.getcwd()))

            self._sendline("", action="don't use other control files")

        case = self._expect(["INPUT TITLE"], action="initialize control title")

        if self.parameters.get("title", None) is not None:
            self._sendline(self.parameters["title"], action="set title")
        else:
            self._sendline("", action="no or default title")

    def _geometry_menu(self, new_coords=True, update=False):
        """
        Sets the options for the molecule geometry menu.

        Args:
            new_coords (bool): if True it expects that the coordinates should be updated in a previously
                existing control file.
            update (bool): whether already present coordinates should be updated or not (i.e. update
                internal coordinates). Ignored if new_coord is True.

        Returns:
            None

        Raises:
            DefineError if incompatibility between current folder state and new_coords.
        """

        # the first case should be when update is False, other cases should only happen
        # when in update, since the coord file is read immediately.
        # The message with SCHOENFLIES SYMBOL happens when a symbol incompatible with the
        # structure is present in control.
        case = self._expect(["SPECIFICATION OF MOLECULAR GEOMETRY",
                             "DO YOU WANT TO CHANGE THE GEOMETRY DATA",
                             "CARTESIAN COORDINATES AND VALUES OF INTERNAL COORDINATES DO  N O T   AGREE",
                             "THESE COORDINATES AND SCHOENFLIES SYMBOL ARE IGNORED"],
                            action="start geometry menu")

        if (not new_coords and case == 0) or (new_coords and case in (1, 2, 3)):
            which_calc = "full define" if new_coords else "update"
            present = "not " if update else ""
            raise DefineError("To run a {} previous data should {}be present".format(which_calc, present))

        if new_coords:
            self._add_atomic_coordinates_from_file()
        else:
            self._existing_coordinates(case, update)

        # if the section should be skipped this will not be used.
        if new_coords or update:
            # if copymo is set, use the symmetry of control file
            if self.parameters.get("copymo", None) is not None:
                self.parameters["sym"] = sdg("$symmetry", os.path.join(self.parameters["copymo"], "control")).strip()

            if self.parameters.get("sym", None) is not None:
                self._define_symmetry()
            elif self.parameters.get("desy", False):
                self._determine_symmetry()

            if self.parameters.get("ired", False):
                self._generate_internal_coordinates()

    def _general_menu(self):
        """
        Runs all the options for the general menu.

        Returns:
            None
        """

        method = self.parameters["method"].lower()

        if method == "dft":
            self._set_dft_options(True)

        elif method == "hf":
            self._set_dft_options(False)

        elif method in ("mp2", "ccsdt", "ccsd(t)", "adc(2)"):
            if method == "ccsdt":
                method = "ccsd(t)"

            # switch off ri for rimp2.
            # ri-flag means ri for ridft in this context
            self.parameters["ri"] = False

            energy_only = self.parameters.get("mp2energy", False)

            # perform mp2 settings
            self._set_mp2_options(energy_only, method)

        else:
            raise DefineParameterError("Unknown method {}".format(method))

        self._set_scf_options()
        self._set_ri_state()

    def _add_atomic_coordinates_from_file(self):
        """
        Adds atomic coordinates from a file.

        Returns:
            None
        """

        filename = self.parameters.get("coord_file", "coord")

        self._sendline("a {}".format(filename), action="set coordinates from file")

        case = self._expect(["CARTESIAN COORDINATES AND VALUES OF INTERNAL COORDINATES DO  N O T   AGREE",
                             "SPECIFICATION OF MOLECULAR GEOMETRY"],
                            action="atomic coordinates from file")

        if case == 0:
            # take cartesian coordinates as default
            # if cartesian and internal do not match
            self._sendline("c", action="choose cartesian coordinates")
            case = self._expect(["SPECIFICATION OF MOLECULAR GEOMETRY"], action="end atomic coordinates from file")

    def _existing_coordinates(self, case, update):
        """
        Handles the coordinates menu for a previously existing control file.

        Args:
            case (int): The first case received from leaving the title menu. Should not be 0.
            update (bool): whether already present coordinates should be updated or not
                (i.e. update internal coordinates).

        Returns:
            None
        """
        choice = "y" if update else "n"

        if case == 1:
            self._sendline(choice, action="change the geometry data")
        elif case == 2:
            self._sendline("c", action="keep cartesian coordinates")
            case = self._expect(["DO YOU WANT TO CHANGE THE GEOMETRY DATA"], action="update coordinates")
            self._sendline(choice, action="change the geometry data")
        elif case == 3:
            raise UnsupportedDefineStateError("Coordinates cannot be simply updated.")
        else:
            raise DefineError("Incompatible choice with update coordinates")

        if update:
            case = self._expect(["SPECIFICATION OF MOLECULAR GEOMETRY"], action="end existing coordinates")

    def _define_symmetry(self):
        """
        Sets the symmetry with the "sy" key when "sym" is among the parameters.
        Adds the tolerance if sym_eps is present.

        Returns:
            None
        """
        self._sendline("sy {} {}".format(self.parameters["sym"], self.parameters.get("sym_eps", "")),
                       action="specify symmetry")
        case = self._expect(["INVALID SCHOENFLIES SYMBOL", "SPECIFICATION OF MOLECULAR GEOMETRY.*"],
                            action="end define symmetry")

        if case == 0:
            raise DefineParameterError("Wrong symmetry symbol: {}".format(self.parameters["sym"]))

    def _determine_symmetry(self):
        """
        Lets define determine the symmetries when the "desy" is present and True.
        Adds the tolerance if "desy_eps" is present.

        Returns:
            None
        """

        action = "determine symmetry"
        if self.parameters.get("desy_eps", None) is None:
            self._sendline("desy", action=action)
        else:
            self._sendline("desy {}".format(self.parameters["desy_eps"]), action=action)
        self._expect(["SPECIFICATION OF MOLECULAR GEOMETRY.*"], action="end determine symmetry")

    def _generate_internal_coordinates(self):
        """
        Generates the redundant internal coordinates when "ired" is among the keys.

        Returns:
            None
        """

        self._sendline("ired", action="generate redundant internal coordinates")

        case = self._expect(["SPECIFICATION OF MOLECULAR GEOMETRY.*", "CAUTION, within .*? steps NO convergence",
                             r"THE B\*m\*Bt-MATRIX IS SINGULAR"],
                            action="end generate internal coordinates")
        if case in (1, 2):
            raise DefineIredError("Some error occurred during the ired procedure.")

    def _switch_to_atomic_attribute_menu(self):
        """
        Switches to atomic attribute menu after setting the coordinates.

        Returns:
            None
        """

        self._sendline("*", action="switch to atomic attribute menu")

        case = self._expect(["IF YOU DO NOT WANT TO USE INTERNAL COORDINATES ENTER*",
                             "ATOMIC ATTRIBUTE DEFINITION MENU.*"], action="switch to atomic attribute menu")

        if case == 0:
            if self.parameters.get("ired", False):
                raise DefineIredError("An incomplete set of internal coordinate is present with ired True")
            self._sendline("no", action="discard internal coords when leaving geometry menu")
            case = self._expect(["ATOMIC ATTRIBUTE DEFINITION MENU.*"], action="end switch to atomic attribute menu")

    def _skip_atomic_attribute_menu(self, from_geometry_menu=False):
        """
        Switches to atomic attribute menu after setting the coordinates.

        Args:
            from_geometry_menu (bool): if True the code start at the moment of closing the
                geometry menu. Otherwise the the geometry menu has been entirely skipped
                before.

        Returns:
            None
        """

        if from_geometry_menu:
            self._sendline("*", action="switch to atomic attribute menu")

        case = self._expect(["ATOMIC ATTRIBUTE DATA.*DO YOU WANT TO CHANGE THESE DATA",
                             "ATOMIC ATTRIBUTE DEFINITION MENU.*",
                             "IF YOU DO NOT WANT TO USE INTERNAL COORDINATES ENTER*"],
                            action="skip to atomic attribute menu")

        # if needed keep cartesian coordinate and expect again the other options
        if case == 2:
            self._sendline("no", action="discard internal coords when leaving geometry menu")
            case = self._expect(["ATOMIC ATTRIBUTE DATA.*DO YOU WANT TO CHANGE THESE DATA",
                                 "ATOMIC ATTRIBUTE DEFINITION MENU.*"], action="end skip to atomic attribute menu")

        if case == 0:
            self._sendline("n", action="do not enter in the atomic attribute menu")
        elif case == 1:
            # here define entered in the definition of the menu, meaning that the
            # atomic properties are not defined. This case is not supported here.
            raise UnsupportedDefineStateError("atomic menu should be already filled, but apparently it is not")

    def _set_basis(self, atom_type, basis):
        """
        Sets the atomic basis set.

        Args:
            atom_type (str): defines the type type of atom(s) addressed. Can be "all", a list "1,2,4-6"
                or the atomic specie "\"c\"" (quotes are required).
            basis (stR): the type of basis.

        Returns:
            None
        """

        self._sendline("b", action="go to basis menu")

        case = self._expect(["ENTER A SET OF ATOMS TO WHICH YOU WANT TO ASSIGN BASIS SETS.*"],
                            action="set basis")

        self._sendline("{} {}".format(atom_type, basis))

        case = self._expect(["ATOMIC ATTRIBUTE DEFINITION MENU.*", "THERE ARE NO DATA SETS CATALOGUED IN FILE"],
                            action="set basis")

        if case == 1:
            raise DefineParameterError("Define did not recognize basis {} for {}".format(basis, atom_type))

    def _define_basis_sets(self):
        """
        Defines the basis set using the "basis" and "basis_atom" keywords in the parameters.

        Returns:
            None
        """

        if self.parameters.get("basis", None) is not None:
            self._set_basis("all", self.parameters["basis"])

        if self.parameters.get("basis_atom", None) is not None:
            for atom_type, basis in self.parameters["basis_atom"].items():
                # convert to string in case an integer slips through since it can also represent an index
                atom_type = str(atom_type).strip()
                # if it is a symbol and does not already contain quotations add them
                if re.fullmatch("[A-Za-z]+", atom_type):
                    atom_type = "\"{}\"".format(atom_type)

                self._set_basis(atom_type, basis)

    def _switch_to_molecular_orbital_definition_menu(self):
        """
        Switches to the molecular orbital definition menu.

        Returns:
            None
        """

        self._sendline("*", action="switch to molecular orbital definition menu")

        case = self._expect([".*MOLECULAR ORBITAL DEFINITION.*"], action="switch to molecular orbital definition menu")

    def _go_to_general_menu(self):
        """
        Returns to the general menu.

        Returns:
            None
        """

        # go back to general menu
        self._sendline("", action="go to general menu")
        case = self._expect(["GENERAL MENU.*"], action="go to general menu")

    def _set_mo_from_control_file(self):
        """
        Sets the molecular orbitals from a previous control file when the path is present in
        the "usemo" key is in the options.
        It checks the actual presence of the file.

        TODO needs to be checked if it works with alpha/beta as well and if it should.

        Returns:
            None
        """

        control_path = self.parameters["usemo"]
        if os.path.isdir(control_path):
            control_path = os.path.join(control_path, "control")

        if not os.path.isfile(control_path):
            raise DefineParameterError("No control file at the path: {}".format(control_path))

        try:
            # Apparently define doesn't like long paths. Try with a relative one if it is too long.
            # The limit seems to be 77 letters. Keep lower value to stay on the safe side.
            if len(control_path) > 70:
                rel_path = os.path.relpath(control_path)

                # if still too long create a symbolic link to the folder
                if len(rel_path) > 70:
                    os.symlink(os.path.dirname(control_path), "tmp_mo")
                    control_path = os.path.join("tmp_mo", os.path.basename(control_path))
                else:
                    control_path = rel_path

            # set path for usemo
            self._sendline("use {}".format(control_path), action="set mo from control file")

            # leave the atomic attribute menu
            case = self._expect(["GENERAL MENU.*", "DO YOU REALLY WANT TO WRITE OUT NATURAL ORBITALS.*",
                                 "LEFT OVER FROM PREVIOUS CALCULATIONS", "DO YOU REALLY WANT TO USE",
                                 "TO CONTINUE, ENTER <return>"],
                                action="back to general menu after use mo")

            if case == 4:
                self._sendline("", action="continue")
                case = self._expect(["GENERAL MENU.*", "DO YOU REALLY WANT TO WRITE OUT NATURAL ORBITALS.*",
                                     "LEFT OVER FROM PREVIOUS CALCULATIONS", "DO YOU REALLY WANT TO USE"],
                                    action="back to general menu after use mo")

            if case in [0, 2, 3]:
                # 0 is if we are already in the general menu, everything is ok.
                # Other cases cover messages when regenerating the mo files.
                pass

            elif case == 1:
                self._go_to_general_menu()

        finally:
            if os.path.isdir("tmp_mo"):
                os.unlink("tmp_mo")

    def _extended_hueckel_theory_menu(self):
        """
        Handles the whole menu of the extended Hueckel theory.

        Returns:
            None
        """

        self._sendline("eht", action="start menu for extended Hueckel theory")

        case = self._expect(["FOUND ATOMIC ORBITALS.*",
                             "DO YOU WANT THE DEFAULT PARAMETERS FOR THE EXTENDED HUECKEL CALCULATION.*",
                             "DO YOU WANT THESE*"],
                            action="initial choice of menu for extended Hueckel theory")

        if case == 0:
            # send default (=y; take orbitals that have been found)
            self._sendline("y", action="default values for Hueckel")
            case = self._expect(["DO YOU WANT THE DEFAULT PARAMETERS FOR THE EXTENDED HUECKEL CALCULATION.*"],
                                action="default values for Hueckel")

        elif case == 1:
            # extended hueckel will be treated in next block
            pass

        elif case == 2:
            # send default (=y)
            self._sendline("y", action="default values for Hueckel")
            case = self._expect(["DO YOU WANT THE DEFAULT PARAMETERS FOR THE EXTENDED HUECKEL CALCULATION.*"],
                                action="default values for Hueckel")

        # send default (=y; use default parameters for extended hueckel calculation)
        self._sendline("y", action="default values for Hueckel")
        case = self._expect(["ENTER THE MOLECULAR CHARGE*"], action="finish setting default values for Hueckel")

        # set molecular charge
        if self.parameters.get("charge", None) is not None:
            # use given charge
            self._sendline(str(self.parameters["charge"]), action="set molecular charge for Hueckel")
        else:
            # use default (=0)
            self._sendline("0", action="set 0 charge for Hueckel")

        case = self._expect(["DO YOU ACCEPT THIS OCCUPATION.*","OCCUPATION NUMBER ASSIGNMENT MENU*"],
                            action="charge for Hueckel")

        unpaired_electrons = self.parameters.get("unpaired_electrons", None)
        if case == 0:
            if unpaired_electrons is not None:
                self._sendline("no", action="don't keep default electrons occupation")
                self._expect(["OCCUPATION NUMBER ASSIGNMENT MENU*"], action="unpaired electrons")
            else:
                # send default (=y; of course we accept this occupation...)
                self._sendline("y", action="don't keep default electrons occupation")

        elif case == 1:
            if unpaired_electrons is None:
                raise DefineParameterError("unpaired_electrons keywork is required if no default is present")

        if unpaired_electrons is not None:
            self._sendline("u {}".format(unpaired_electrons), action="set unpaired electrons")
            case = self._expect(["OCCUPATION NUMBER ASSIGNMENT MENU*"], action="set unpaired electrons")

            self._sendline("*", action="leave the unpaired electron menu")

        # leave the atomic attribute menu
        case = self._expect(["GENERAL MENU.*", "DO YOU REALLY WANT TO WRITE OUT NATURAL ORBITALS.*",
                             "LEFT OVER FROM PREVIOUS CALCULATIONS", "DO YOU REALLY WANT TO USE"],
                            action="back to general menu after huckel")

        if case in [0, 2, 3]:
            # 0 is if we are already in the general menu, everything is ok.
            # Other cases cover messages when regenerating the mo files.
            pass

        elif case == 1:
            self._go_to_general_menu()

    def _excited_state_menu(self):
        """
        Handles the whole excited state menu.

        Returns:
            None
        """

        # check if values are set, set default values if not
        ex_method = self.parameters.get("ex_method", "rpa")
        ex_multi = self.parameters.get("ex_multi", "singlet")
        ex_irrep_states = self.parameters.get("ex_irrep_states", None)
        ex_all_states = self.parameters.get("ex_all_states", None if ex_irrep_states else 10)
        if self.parameters.get("ex_frequency", None) is None:
            # d-linie in nm for natrium
            ex_frequency = "589"
            ex_frequency_unit = "nm"
        else:
            ex_frequency = self.parameters["ex_frequency"]
            ex_frequency_unit = self.parameters.get("ex_frequency_unit", "nm")

        if self.parameters.get("ex_exopt", None) is not None:
            cdg("$exopt", "    {}".format(self.parameters.get("ex_exopt")))

        self._sendline("ex", action="move to response calculations menu")
        case = self._expect(["MAIN MENU FOR RESPONSE CALCULATIONS"], action="start excited states")

        # search the match with the list of available response calculation options
        case = self._expect(["----.*ENTER <OPTION> TO SWITCH ON/OFF OPTION"],
                            action="identify available response calculations")

        response_opt_avail = []
        for l in self.define.match.group().splitlines():
            l = str(l, "utf-8").strip()
            if not l or "OPTION" in l or "----" in l:
                continue
            split = [x.strip() for x in l.split('|')]
            response_opt_avail.append(split[0])

        # we assume that in the list of available options there will be either only "urpa" or
        # both "rpas" and "rpat". Here we check that this is actually the case.
        # Potentially this should not happen, but check to avoid setting the wrong options.
        # the same for "ucis", "ciss" "cist"
        for m in ("cis", "rpa"):
            opt1 = "u{}".format(m) in response_opt_avail
            opt2 = "{}s".format(m) in response_opt_avail and "{}t".format(m) in response_opt_avail
            # this is a NOT XOR as we want exactly one of the cases to be true
            if opt1 == opt2:
                raise UnsupportedDefineStateError("The current list response options is not compatible "
                                                  "with the choices for cis and rpa: {}".format(response_opt_avail))

        if ex_method in ("cis", "rpa"):
            if "u{}".format(ex_method) in response_opt_avail:
                option = "u{}".format(ex_method)
            else:
                option = "{}{}".format(ex_method, ex_multi[0])

        elif ex_method in ("dynpol", "polly"):
            option = str(ex_method)
        else:
            raise DefineParameterError("ex_method not among the supported values: {}".format(ex_method))

        # this check should be redundant, but an error can still occur
        if option not in response_opt_avail:
            msg = "Option {} is not present among those available: {}".format(option, response_opt_avail)
            raise DefineParameterError(msg)

        self._sendline(str(option), action="set excited state option")
        # search the match with the list of available response calculation options
        case = self._expect(["----.*ENTER <OPTION> TO SWITCH ON/OFF OPTION"],
                            action="check that the option has been activated")

        for l in self.define.match.group().splitlines():
            l = str(l, "utf-8").strip()
            split = [x.strip() for x in l.split('|')]
            if split[0] == option and split[1] == "on":
                break
        else:
            raise DefineError("Option {} was not activated".format(option))

        self._sendline("*", action="leave response calculation main menu")

        case = self._expect(["STATE SELECTION MENU",
                             "SELECT FREQUENCY UNIT*",
                             "GENERAL OPTIONS FOR RESPONSE CALCULATIONS*"],
                            action="after response calculation main menu")
        if case == 0:

            # search the match with the list of available irreps
            case = self._expect(["----.*SELECT IRREP AND NUMBER OF STATES"], action="identify available irreps")
            irreps_avail = {}
            for l in self.define.match.group().splitlines():
                l = str(l, "utf-8").strip()
                if not l or "IRREP" in l or "----" in l:
                    continue
                split = l.split("|")
                # The line can also contain "    |   |   | xz yz", being the continuation of the line above
                # Skip all the lines whose second split cannot be converted to int
                try:
                    num_states = int(split[1].strip())
                except ValueError as e:
                    continue

                irreps_avail[split[0].strip()] = num_states

            # only use number of states that are <= than those available
            selected_irreps = {}
            if ex_all_states is not None:
                for i, n in irreps_avail.items():
                    selected_irreps[i] = ex_all_states if ex_all_states <= n else n

            if ex_irrep_states:
                for i, n in ex_irrep_states.items():
                    if i in irreps_avail:
                        selected_irreps[i] = n if n <= irreps_avail[i] else irreps_avail[i]

            if not selected_irreps:
                raise DefineParameterError("No excited state could be set given the input parameters")

            for i, n in selected_irreps.items():
                self._sendline("{} {}".format(i, n), action="set excited states number for specific irrep")
                self._expect(["STATE SELECTION MENU*"], action="back to state selection menu")

            self._sendline("*", action="leave excited states menu")

            case = self._expect(["GENERAL OPTIONS FOR RESPONSE CALCULATIONS*"],
                                action="back to response calulation menu")

        elif case == 1:
            self._sendline(str(ex_frequency_unit), action="set excitation units")
            case = self._expect(["FREQUENCY INPUT MENU*"], action="set excitations")

            self._sendline("a")
            self._expect(["INPUT FREQUENCY NO.*"], action="set excitations")
            self._sendline(str(ex_frequency), action="set excitation frequency")
            self._expect(["INPUT FREQUENCY NO.*"], action="set excitations")
            self._sendline("*", action="leave frequency menu")
            self._expect(["FREQUENCY INPUT MENU*"], action="set excitations")
            self._sendline("*", action="leave excited states menu")

            case = self._expect(["GENERAL OPTIONS FOR RESPONSE CALCULATIONS*"],
                                action="back to response calulation menu")

        elif case == 2:
            # already at the end of the other cases
            pass

        self._sendline("*", action="leave excited states menu")
        self._expect(["SET SCF DENSITY CONVERGENCE THRESHOLD*"], action="excited states: scf")
        self._sendline("y", action="default for density convergence threshold")

    def _set_ri_state(self):
        """
        Sets the parameters for ri calculation.

        Returns:
            None
        """

        if self.parameters.get("rijk", False):
            self._sendline("rijk", action="choose rijk option")
            self._expect(["ENTER RI-OPTION TO BE MODIFIED.*"], action="ri menu")

            self._sendline("on", action="activate rijk option")
            self._expect(["ENTER RI-OPTION TO BE MODIFIED.*"], action="ri menu")

            self._go_to_general_menu()

        elif self.parameters.get("ri", False):
            self._sendline("ri", action="choose ri option")
            self._expect(["ENTER RI-OPTION TO BE MODIFIED.*"], action="ri menu")

            self._sendline("on", action="activate ri option")
            self._expect(["ENTER RI-OPTION TO BE MODIFIED.*"], action="ri menu")

            self._go_to_general_menu()

            if self.parameters.get("marij", False):
                self._sendline("marij", action="choose marij option")
                self._expect(["Enter the number to change a value or <return> to accept all.*"],
                             action = "ri menu")

                self._sendline("", "accept all")
                self._expect(["GENERAL MENU.*"], action="back to general menu")

    def _set_dft_options(self, use_dft=True):
        """
        Sets the dft options.

        Returns:
            None
        """

        self._sendline("dft", action="move to dft menu")
        case = self._expect(["ENTER DFT-OPTION TO BE MODIFIED.*"], action="start dft options")

        if use_dft:
            self._sendline("on", action="enable dft")
            case = self._expect(["ENTER DFT-OPTION TO BE MODIFIED.*"], action="back to dft options")

            functional = self.parameters.get("functional", "b-p")
            self._sendline("func {}".format(functional), action="set functional")
            case = self._expect(["ENTER DFT-OPTION TO BE MODIFIED.*", "SPECIFIED FUNCTIONAL not SUPPORTED"],
                                action="back to dft options")
            if case == 1:
                raise DefineParameterError("XC functional not supported by define: {}".format(functional))

            if self.parameters.get("gridsize", None) is not None:
                self._sendline("grid {}".format(self.parameters["gridsize"]), action="set gridsize")
                self._expect(["ENTER DFT-OPTION TO BE MODIFIED.*"], action="back to dft options")
        else:
            self._sendline("off", action="disable dft")
            case = self._expect(["ENTER DFT-OPTION TO BE MODIFIED.*"], action="back to dft options")

        self._go_to_general_menu()

    def _set_mp2_options(self, energy_only=True, method="mp2"):
        """
        Sets mp2 options.

        Args:
            energy_only (bool): if true the job is energy only.
            method (str): the method used for the mp2 calculation.

        Returns:
            None
        """

        # change to mp2 menu
        self._sendline("cc", action="change to mp2 menu")
        self._expect(["INPUT MENU FOR CALCULATIONS WITH ricc2*"], action="start ricc2 options")

        # go to freeze menu and back again to let
        # define create the correct default settings for "freeze"
        self._sendline("freeze", action="go to freeze menu")
        self._expect(["frozen core / frozen virtual assignment menu*"], action="at free menu")
        self._sendline("*", action="quit freeze menu")
        self._expect(["INPUT MENU FOR CALCULATIONS WITH ricc2*"], action="back to ricc2 options")

        # go to cbas menu and back again to let
        # define create the references to the auxiliary basis sets
        self._sendline("cbas", action="go to cbas menu")
        self._expect(["AUXILIARY BASIS SET DEFINITION MENU*"], action="at cbas menu")
        self._sendline("*", action="quit cbas menu")
        self._expect(["INPUT MENU FOR CALCULATIONS WITH ricc2*"], action="back to ricc2 options")

        if self.parameters.get("maxcor", None) is not None:
            # set maximum core memory (maxcor)
            self._sendline("memory", action="input memory")
            self._sendline(str(self.parameters["maxcor"]), action="memory value")
            self._expect(["INPUT MENU FOR CALCULATIONS WITH ricc2*"], action="back to ricc2 options")

        if method == "adc(2)":
            if self.parameters.get("ex_mp2", None) is not None:
                self._sendline("exci", action="go to excited states menu")
                for irrep, l in self.parameters["ex_mp2"].items():
                    irrep_cmd = "irrep={} multiplicity={} nexc={}".format(irrep, l[0], l[1])
                    self._sendline(irrep_cmd,
                                   action="set excited states")
                    case = self._expect(["unkown irrep", "incomprehensible", "error", "irrep"],
                                        action="set excited states")

                    if case in (0, 1, 2):
                        raise DefineParameterError('Error when providing the irrep "{}" in ex_mp2'.format(irrep_cmd))

                self._sendline("spectrum  states=all", action="set spectrum")
                self._sendline("exprop  states=all unrelaxed", action="set exprop")
                self._sendline("*", action="quit excited states menu")
                self._expect(["INPUT MENU FOR CALCULATIONS WITH ricc2*"], action="back to ricc2 options")

        if self.parameters.get("use_f12", False):
            self._sendline("f12", action="enter f12 menu")
            case = self._expect(["INPUT MENU FOR MP2-F12 CALCULATIONS"], action="f12 menu")
            self._sendline("*", action="quit f12 menu")
            case = self._expect(["INPUT MENU FOR CALCULATIONS WITH ricc2*"], action="back to ricc2 options")

            self._sendline("cabs", action="go to cabs menu")
            case = self._expect(["AUXILIARY BASIS SET DEFINITION MENU"], action="at cabs menu")
            self._sendline("*", action="quit cabs menu")
            case = self._expect(["INPUT MENU FOR CALCULATIONS WITH ricc2*"], action="back to ricc2 options")

            # Note: it seems that the line ccsdapprox ccsd(f12*) is already there with just use_f12
            # this option is kept for safety.
            if self.parameters.get("use_f12*", False) and method == "ccsd(t)":
                mdgo("rir12", {"ccsdapprox": "ccsdapprox  ccsd(f12*)"})

        # go to submenu "ricc2"
        self._sendline("ricc2", action="go to ricc2 menu")
        if energy_only:
            # set mp2energy
            self._sendline(method, action="set method")
        else:
            # set geopt mp2energy
            self._sendline("geoopt {}".format(method), action="set method")

        if self.parameters.get("maxiter", None) is not None:
            self._sendline("maxiter {}".format(self.parameters["maxiter"]), action="set maximum iterations")

        # go back to mp2 menu
        self._sendline("*", action="go back to mp2 menu")
        case = self._expect(["INPUT MENU FOR CALCULATIONS WITH ricc2*"], action="back to ricc2 options")

        self._sendline("*", action="go to general menu")
        case = self._expect(["GENERAL MENU.*"], action="go to general menu")

    def _switch_to_scf_menu(self):
        """
        Goes back to scf menu.

        Returns:
            None
        """

        self._sendline("scf", action="go back to scf menu")
        case = self._expect(["ENTER SCF-OPTION TO BE MODIFIED.*"], action="at scf menu")

    def _set_scf_options(self):
        """
        Sets all the options for scf.

        Returns:
            None
        """

        self._switch_to_scf_menu()

        if self.parameters.get("scfconv", None) is not None:
            self._sendline("conv", action="go to conv menu")
            self._expect(["ENTER DESIRED ACCURACY OF SCF-ENERGY.*"], action="at scf-energy menu")

            self._sendline(str(self.parameters["scfconv"]), action="set scfconv")
            self._expect(["ENTER SCF-OPTION TO BE MODIFIED.*"], action="back to scf menu")

        if self.parameters.get("scfiterlimit", None) is not None:
            self._sendline("iter", action="go to iter")
            case = self._expect(["ENTER NEW VALUE FOR MAXIMUM NUMBER OF SCF-ITERATIONS.*"],
                                action="scf iterations")

            self._sendline(str(self.parameters["scfiterlimit"]), action="insert maximum number of operations")
            case = self._expect(["ENTER SCF-OPTION TO BE MODIFIED.*"], action="back to scf menu")

        self._go_to_general_menu()

    def _quit_general_menu(self):
        """
        Quits the general menu, exiting from define.

        Returns:
            None
        """

        self._sendline("*", action="quit general menu")

    def _copy_mo_files(self):
        """
        Copies the mos, alpha and beta files from the "copymo" path given in the parameters "copymo".

        Returns:
            None
        """

        path = None

        # if flag copymo - copy files (mos/alpha/beta) to current folder
        if self.parameters["copymo"]:

            path = self.parameters["copymo"]
            source_mos = os.path.join(path, "mos")
            source_alpha = os.path.join(path, "alpha")
            source_beta = os.path.join(path, "beta")

            if not os.path.isfile(source_mos) and \
                    (not os.path.isfile(source_alpha) and not os.path.isfile(source_beta)):
                raise DefineParameterError("The files that should be copied do not exist")

            dest_mos = os.path.join(self.workdir, "mos")
            dest_alpha = os.path.join(self.workdir, "alpha")
            dest_beta = os.path.join(self.workdir, "beta")

            # copy mos file from copymo dir to calculation dir if exist in both folders
            if os.path.isfile(source_mos) and os.path.isfile(dest_mos):
                shutil.copy(source_mos, self.workdir)
            # copy alpha and beta file from copymo dir to calculation dir if exist in both folders
            if os.path.isfile(source_alpha) and os.path.isfile(dest_alpha) and os.path.isfile(source_beta) \
                    and os.path.isfile(dest_beta):
                shutil.copy(source_alpha, self.workdir)
                shutil.copy(source_beta, self.workdir)

            # manipulate mos file in copymo dir and save it to alpha and beta
            # if mos file in copymo dir an alpha/beta file in calculation dir
            if os.path.isfile(source_mos) and os.path.isfile(dest_alpha) and os.path.isfile(dest_beta):
                with open(source_mos, "r") as f:
                    new_alpha = f.readlines()
                new_beta = list(new_alpha)

                for i in range(len(new_alpha)):
                    if "scfmo" in new_alpha[i]:
                        new_alpha[i] = new_alpha[i].replace("scfmo", "uhfmo_alpha")
                        new_beta[i] = new_beta[i].replace("scfmo", "uhfmo_beta")
                        break

                with open(dest_alpha, "w") as f:
                    f.writelines(new_alpha)

                with open(dest_beta, "w") as f:
                    f.writelines(new_beta)

    def _add_cosmo(self):
        """
        Adds $cosmo datagroup based on the parameters available.

        Returns:
            None
        """

        cosmo_keys = ["epsilon", "nppa", "nspa", "disex", "rsolv", "routf", "cavity"]

        d = {k: v for k, v in self.parameters.items() if k in cosmo_keys}

        c = Control.from_file("control")
        c.add_cosmo(**d)
        c.to_file("control")

    def _postprocess(self):
        """
        Adds keywords to the control file that should be set outside define:
            * cosmo
            * disp

        Returns:
            None
        """

        if self.parameters.get("use_cosmo", False):
            self._add_cosmo()

        c = Control.from_file("control")
        c.set_disp(self.parameters.get("disp", None))
        c.to_file("control")

    def dump_command_file(self, filepath="define_commands"):
        """
        Utility function that dumps a file with the commands given to define.
        Can be used as an input for define to test its behavior:
        define < define_commands

        Args:
            filepath: path to the file that should be written.

        Returns:
            None
        """
        with open(filepath, "w") as f:
            # add one file empty space to be sure of sending the last command
            f.write("\n".join(self.sent_commands+[""]))
