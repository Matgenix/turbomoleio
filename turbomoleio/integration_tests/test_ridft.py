import pytest

from turbomoleio.testfiles.utils import run_itest
from turbomoleio.input.utils import get_define_template
from turbomoleio.output.files import ScfOutput

# structures = ['aceton', 'ch4', 'h2o', 'h3cbr', 'methanol', 'nh3', 'phenol', 'sf4', 'sih4']
structures = ['h2o', 'nh3']


@pytest.mark.integration
class TestRidft:

    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft(self, structure):

        assert run_itest("ridft", get_define_template("ridft"), structure, "ridft_{}_std".format(structure), ScfOutput)

    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft_rijk(self, structure):
        define_opt = get_define_template("ridft")
        define_opt["ri"] = False
        define_opt["rijk"] = True

        assert run_itest("ridft", define_opt, structure, "ridft_{}_rik".format(structure), ScfOutput)

    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft_marij(self, structure):
        define_opt = get_define_template("ridft")
        define_opt["ri"] = True
        define_opt["marij"] = True

        assert run_itest("ridft", define_opt, structure, "ridft_{}_marij".format(structure), ScfOutput)

    @pytest.mark.parametrize("structure", structures)
    def test_run_ridft_fermi(self, structure):
        define_opt = get_define_template("ridft")
        define_opt["ri"] = False
        define_opt["rijk"] = True
        datagroups_options = {"fermi": "tmstrt=500 tmend=50 tmfac=0.9 hlcrt=0.2"}

        assert run_itest("ridft", define_opt, structure, "ridft_{}_fermi".format(structure), ScfOutput,
                         datagroups_options=datagroups_options)
