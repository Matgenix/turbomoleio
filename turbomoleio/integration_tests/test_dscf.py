import os
import pytest

from turbomoleio.testfiles.utils import run_itest, get_tfp
from turbomoleio.input.utils import get_define_template
from turbomoleio.output.files import ScfOutput

# structures = ['aceton', 'ch4', 'h2o', 'h3cbr', 'methanol', 'nh3', 'phenol', 'sf4', 'sih4']
structures = ['h2o', 'nh3']

disp_corrections = ["DFT-D1", "DFT-D2", "DFT-D3", "DFT-D3 BJ"]


@pytest.mark.integration
class TestDscf:

    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf(self, structure):

        assert run_itest("dscf", get_define_template("dscf"), structure, "dscf_{}_std".format(structure), ScfOutput)

    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf_ired(self, structure):
        define_opt = get_define_template("dscf")
        define_opt["ired"] = True

        assert run_itest("dscf", define_opt, structure, "dscf_{}_ired".format(structure), ScfOutput)

    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf_hf(self, structure):

        define_opt = get_define_template("dscf")
        define_opt["method"] = "hf"
        define_opt["desy"] = False
        define_opt["sym"] = "c1"
        define_opt["sym_eps"] = 0.001
        define_opt["scfiterlimit"] = 300
        define_opt["scfconv"] = 6
        define_opt["title"] = None

        assert run_itest("dscf", define_opt, structure, "dscf_{}_hf".format(structure), ScfOutput)

    def test_run_dscf_usemo(self):
        """
        Tests the usemo functionalities with mos file
        """

        define_opt = get_define_template("dscf")

        define_opt["basis"] = "def2-SV(P)"
        define_opt["usemo"] = os.path.join(get_tfp(), "mo", "mos_nh3_nosym")

        assert run_itest("dscf", define_opt, "nh3", "dscf_nh3_usemo", ScfOutput)

    # TODO check if alpha/beta for usemo are needed
    # def test_run_dscf_usemo_alpha(self):
    #     """
    #     Tests the usemo functionalities with alpha and beta files in the case of unpaired electrons
    #     """
    #
    #     define_opt = get_define_template("dscf")
    #
    #     define_opt["basis"] = "def2-SV(P)"
    #     define_opt["charge"] = 1
    #     define_opt["unpaired_electrons"] = 1
    #     define_opt["usemo"] = os.path.join(get_tfp(), "mo", "alpha_beta_nh3")
    #
    #     assert run_itest("dscf", define_opt, "nh3", "dscf_nh3_usemo_alpha", {})

    def test_run_dscf_copymo(self):
        """
        Tests the copy functionalities with mos file
        """

        define_opt = get_define_template("dscf")

        define_opt["basis"] = "def2-SV(P)"
        define_opt["copymo"] = os.path.join(get_tfp(), "mo", "mos_nh3")

        assert run_itest("dscf", define_opt, "nh3", "dscf_nh3_copymo_mos", ScfOutput)

    def test_run_dscf_copymo_mos_unpaired(self):
        """
        Tests the copy functionalities with mos file
        """

        define_opt = get_define_template("dscf")

        define_opt["basis"] = "def2-SV(P)"
        define_opt["charge"] = 1
        define_opt["unpaired_electrons"] = 1
        define_opt["copymo"] = os.path.join(get_tfp(), "mo", "mos_nh3")

        assert run_itest("dscf", define_opt, "nh3", "dscf_nh3_copymo_mos_unpaired", ScfOutput)

    def test_run_dscf_copymo_alpha_unpaired(self):
        """
        Tests the copy functionalities with mos file
        """

        define_opt = get_define_template("dscf")

        define_opt["basis"] = "def2-SV(P)"
        define_opt["charge"] = 1
        define_opt["unpaired_electrons"] = 1
        define_opt["copymo"] = os.path.join(get_tfp(), "mo", "alpha_beta_nh3")

        assert run_itest("dscf", define_opt, "nh3", "dscf_nh3_copymo_alpha_unpaired", ScfOutput)

    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf_cosmo(self, structure):

        define_opt = get_define_template("dscf")
        define_opt["use_cosmo"] = True
        define_opt.update({"epsilon": 60.0, "nppa": 92, "nspa": None, "disex": 0, "rsolv": 1.3,
         "routf": None, "cavity": "open"})
        assert run_itest("dscf", define_opt, structure, "dscf_{}_cosmo".format(structure), ScfOutput)

    @pytest.mark.parametrize("disp", disp_corrections)
    def test_run_dscf_disp(self, disp):

        define_opt = get_define_template("dscf")
        define_opt["disp"] = disp
        disp_fn = "dscf_nh3_disp_{}".format("_".join(disp.split()))
        assert run_itest("dscf", define_opt, "nh3", disp_fn, ScfOutput)
