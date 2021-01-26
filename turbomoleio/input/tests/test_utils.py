import pytest
from turbomoleio.input.utils import get_define_template, validate_parameters


def test_get_define_template():
    """Testing get define template util function."""

    dscf_dict = get_define_template("dscf")
    assert dscf_dict["basis"] == "def-SV(P)"

    with pytest.raises(ValueError,
                       match=r'^Could not find template file '
                             r'\S*non_existing_template.yaml$'):
        get_define_template("non_existing_template.yaml")


def test_validate_parameters():

    assert validate_parameters({})

    assert validate_parameters(get_define_template("dscf"))
    assert validate_parameters(get_define_template("ridft"))
    assert validate_parameters(get_define_template("dscf_escf"))
    assert validate_parameters(get_define_template("ridft_escf"))
    assert validate_parameters(get_define_template("ridft_rimp2"))

    assert not validate_parameters({"fake_key": 1})

    # check dependency
    d = {"use_f12*": True}
    assert not validate_parameters(d)
    d["use_f12"] = True
    assert not validate_parameters(d)
    d["method"] = "mp2"
    assert validate_parameters(d)

    #wrong type
    assert not validate_parameters({"method": 1})

    assert validate_parameters({"disp": "DFT-D1"})
    # value not in possible options
    assert not validate_parameters({"disp": "wrong_value"})
