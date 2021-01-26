import pytest

from turbomoleio.testfiles.utils import run_itest
from turbomoleio.input.utils import get_define_template
from turbomoleio.output.files import JobexOutput

structures = ['h2o', 'nh3']


@pytest.mark.integration
class TestJobex:

    @pytest.mark.parametrize("structure", structures)
    def test_run_jobex_dscf(self, structure):
        dp = get_define_template("dscf")
        dp["desy"] = True
        dp["ired"] = True

        assert run_itest("jobex", dp, structure, "jobex_dscf_{}_sym".format(structure),
                         JobexOutput, arguments="-c 2")

    # @pytest.mark.parametrize("structure", structures)
    # def test_run_ridft_rdgrad_relax(self, structure):
    #     dp = get_define_template("ridft")
    #     dp["desy"] = True
    #     dp["ired"] = True
    #
    #     assert run_itest(["ridft", "rdgrad", "relax"], dp,
    #                      structure, "ridft_rdgrad_relax_{}_sym".format(structure),
    #                      [ScfOutput, GradOutput, StatptOutput])
    #
    # @pytest.mark.parametrize("structure", structures)
    # def test_run_dscf_egrad_relax(self, structure):
    #     dp = get_define_template("dscf_escf")
    #     dp["desy"] = False
    #     dp["ired"] = False
    #
    #     assert run_itest(["dscf", "egrad", "relax"], dp,
    #                      structure, "dscf_escf_relax_{}_nosym".format(structure),
    #                      [ScfOutput, EgradOutput, RelaxOutput])
    #
    # @pytest.mark.parametrize("structure", structures)
    # def test_run_ridft_egrad_statpt_ex_irrep(self, structure):
    #     dp = get_define_template("ridft_escf")
    #     dp["desy"] = True
    #     dp["ired"] = True
    #     dp["ex_all_states"] = None
    #     dp["ex_irrep_states"] = {"a1": 1}
    #
    #     assert run_itest(["ridft", "egrad", "statpt"], dp,
    #                      structure, "ridft_egrad_statpt_{}_ex_irrep".format(structure),
    #                      [ScfOutput, EgradOutput, RelaxOutput])
