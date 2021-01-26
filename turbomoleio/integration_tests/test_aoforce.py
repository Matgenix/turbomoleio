import pytest

from turbomoleio.testfiles.utils import run_itest
from turbomoleio.input.utils import get_define_template
from turbomoleio.output.files import ScfOutput, AoforceOutput

structures = ['h2o', 'nh3']


@pytest.mark.integration
class TestAoforce:

    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft_aoforce_nosym(self, structure):
        dp = get_define_template("ridft")
        dp["desy"] = False

        assert run_itest(["ridft", "aoforce"], dp,
                         structure, "ridft_aoforce_{}_nosym".format(structure), [ScfOutput, AoforceOutput])

    @pytest.mark.parametrize("structure", structures)
    def test_run_dscf_aoforce_sym(self, structure):
        dp = get_define_template("dscf")
        dp["desy"] = True

        assert run_itest(["dscf", "aoforce"], dp,
                         structure, "dscf_aoforce_{}_sym".format(structure), [ScfOutput, AoforceOutput])
