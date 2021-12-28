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

import os
import cerberus
from monty.serialization import loadfn
from turbomoleio.input.define import DefineError


def get_define_template(name):
    """
    Returns the dict generate from the default templates available.

    Args:
        name (str): name of the template present in the folder.

    Returns:
        dict: the template converted from the yaml file.
    """

    if not name.endswith(".yaml"):
        name += ".yaml"

    filepath = os.path.join(os.path.abspath(os.path.dirname(__file__)), "templates", name)

    if not os.path.isfile(filepath):
        raise ValueError("Could not find template file {}".format(filepath))

    return loadfn(filepath)


class ParametersValidationError(DefineError):
    """
    Exception raised when the validation of the parameters fails due
    with respect to the defined schema.
    """
    pass


def ex_setter(value):
    """
    Helper function for setting excited states in the cerberus schema
    of DefineRunner.
    """
    def setter(doc):
        if doc.get("ex_method", None):
            return value
        else:
            return None

    return setter


ex_method = ["rpa", "cis", "dynpol", "polly"]
ex_multi = ["singlet", "doublet", "triplet"]
mp2_calc = ["mp2", "adc(2)", "ccsd(t)", "ccsdt"]


#: The schema for the validation of the DefineRunner parameters. Based on cerberus.
schema_define_params = {
    "title": {"type": "string",
              "nullable": True,
              "default": ""},
    "metric": {"type": "integer",
               "nullable": True,
               "default": None},
    "copymo": {"type": "string",
               "nullable": True,
               "default": None},
    "sym": {"type": "string",
            "nullable": True,
            "default": None},
    "sym_eps": {"type": "float",
                "nullable": True,
                "default": None},
    "desy": {"type": "boolean",
             "default": True},
    "desy_eps": {"type": "float",
                 "nullable": True,
                 "dependencies": {"desy": True}, "default": None},
    "ired": {"type": "boolean",
             "default": True},
    "usemo":  {"type": "string",
               "nullable": True,
               "default": None},
    "ex_method":  {"type": "string",
                   "nullable": True,
                   "allowed": ex_method,
                   "dependencies": {"ex_multi": ex_multi}, "default": None},
    "ex_multi": {"type": "string",
                 "nullable": True,
                 "allowed": ex_multi,
                 "dependencies": {"ex_method": ex_method}, "default": None},
    "ex_all_states": {"type": "integer",
                      "nullable": True,
                      "min": 0,
                      "dependencies": ["ex_method", "ex_multi"],
                      "default_setter": ex_setter(10)},
    "ex_irrep_states": {"type": "dict",
                        "nullable": True,
                        "keysrules": {"type": "string"},
                        "valuesrules": {"type": "integer"},
                        "dependencies": ["ex_method", "ex_multi"]},
    "ex_mp2": {"type": "dict",
               "nullable": True,
               "keysrules": {"type": "string"},
               "valuesrules": {"type": "list", "valuesrules": {"type": "integer", "min": 0}},
               "dependencies": ["ex_method", "ex_multi"]},
    "ex_frequency": {"type": "float",
                     "nullable": True,
                     "default_setter": ex_setter(589)},
    "ex_frequency_unit": {"type": "string",
                          "nullable": True,
                          "default_setter": ex_setter("nm")},
    "ex_exopt": {"type": "integer",
                 "nullable": True,
                 "default": None,
                 "min": 0},
    "method":  {"type": "string",
                "default": "dft",
                "allowed": ["dft", "hf", "mp2", "adc(2)", "ccsd(t)", "ccsdt"]},
    "mp2energy": {"type": "integer",
                  "nullable": True,
                  "default": None},
    "basis": {"type": "string",
              "default": None,
              "nullable": True},
    "basis_atom": {"type": "dict",
                   "nullable": True,
                   "keysrules": {"type": "string"},
                   "valuesrules": {"type": "string"}},
    "charge": {"type": "integer",
               "nullable": True,
               "default": None},
    "unpaired_electrons": {"type": "integer",
                           "nullable": True,
                           "default": None},
    "rijk": {"type": "boolean",
            "default": False,
            "oneof": [{"allowed": [True], "dependencies": {"ri": False}},
                      {"allowed": [False]}]},
    "ri": {"type": "boolean",
           "default": False,
           "oneof": [{"allowed": [True], "dependencies": {"rijk": False}},
                     {"allowed": [False]}]},
    "marij": {"type": "boolean",
              "default": False,
              "oneof": [{"allowed": [True], "dependencies": {"rijk": False, "ri": True}},
                        {"allowed": [False]}]},
    "functional": {"type": "string",
                   "nullable": True,
                   "default": None},
    "gridsize": {"type": "string",
                 "nullable": True,
                 "default": None},
    "maxcor": {"type": "float",
               "nullable": True,
               "default": None,
               "dependencies": {"method": mp2_calc}},
    "use_f12": {"type": "boolean",
                "default": False,
                "oneof": [{"allowed": [True], "dependencies": {"method": mp2_calc}},
                          {"allowed": [False]}]},
    "use_f12*": {"type": "boolean",
                 "default": False,
                 "oneof": [{"allowed": [True], "dependencies": {"method": mp2_calc, "use_f12": True}},
                           {"allowed": [False]}]},
    "maxiter": {"type": "integer",
                "nullable": True,
                "default": None,
                "dependencies": {"method": mp2_calc}},
    "scfiterlimit": {"type": "integer",
                     "default": 200},
    "scfconv": {"type": "integer",
                "nullable": True,
                "min":4,
                "max": 9},
    "coord_file": {"type": "string",
                   "nullable": True,
                   "default": None},
    "disp": {"type": "string",
             "nullable": True,
             "allowed": ["DFT-D1", "DFT-D2", "DFT-D3", "DFT-D3 BJ"],
             "default": None},
    "use_cosmo": {"type": "boolean",
                  "default": False},
    "epsilon": {"type": "float",
                "nullable": True,
                "default": None,
                "dependencies": {"use_cosmo": True}},
    "nppa": {"type": "integer",
             "nullable": True,
             "default": None,
             "min": 0,
             "dependencies": {"use_cosmo": True}},
    "nspa": {"type": "integer",
             "nullable": True,
             "default": None,
             "min": 0,
             "dependencies": {"use_cosmo": True}},
    "disex": {"type": "float",
              "nullable": True,
              "default": None,
              "min": 0,
              "dependencies": {"use_cosmo": True}},
    "rsolv": {"type": "float",
              "nullable": True,
              "default": None,
              "min": 0,
             "dependencies": {"use_cosmo": True}},
    "routf": {"type": "float",
              "nullable": True,
              "default": None,
              "dependencies": {"use_cosmo": True}},
    "cavity": {"type": "string",
               "nullable": True,
               "default": None,
               "dependencies": {"use_cosmo": True}},
    "use_old_amat": {"type": "boolean",
                     "nullable": True,
                     "default": None,
                     "dependencies": {"use_cosmo": True}}
}


#: An instance of cerberus.Validator that can be used to validate the input for DefineRunner.
define_parameters_validator = cerberus.Validator(schema=schema_define_params, ignore_none_values=True)


def validate_parameters(parameters):
    """
    Function using the cerberus tool and the schema_define_params dictionary
    to validate the dictionary that should be passed to DefineRunner in the
    parameters argument.
    Consider it as an experimental feature.

    Args:
        parameters (dict): the dictionary that should be validate.

    Returns:
        bool: True if valid according to the schema defined, False otherwise.
    """
    return define_parameters_validator.validate(parameters)
