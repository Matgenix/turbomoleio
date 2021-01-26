"""
The parsing of all the Data objects is done in the test_files.py module and
is not repeated here. In this module only specific functions of the Data
objects are tested.
"""
import os
import pytest

from turbomoleio.output.data import TurbomoleData, ScfIterationData, AoforceVibrationalData
from turbomoleio.testfiles.utils import has_matplotlib


def test_str(testdir):

    path = os.path.join(testdir, "outputs", "dscf", "h2o_std.log")
    with open(path) as f:
        log = f.read()

    td = TurbomoleData.from_string(log)
    s = str(td)
    assert "Content of the TurbomoleData object" in s
    assert "dscf" in s


def test_ScfIterationData(testdir):
    path = os.path.join(testdir, "outputs", "dscf", "h2o_std.log")
    sid = ScfIterationData.from_file(path)

    if has_matplotlib():
        assert sid.plot_energies(show=False)


def test_AoforceVibrationalData(testdir):
    path = os.path.join(testdir, "outputs", "aoforce", "aceton_full.log")
    avd = AoforceVibrationalData.from_file(path)

    assert avd.n_negative_freqs(tol=0.1) == 1
    assert avd.n_zero_freqs(tol=0.1) == 6
    assert avd.n_positive_freqs(tol=0.1) == 23

    df = avd.get_freqs_df()
    assert df["frequency"][0] == pytest.approx(-113.79)
