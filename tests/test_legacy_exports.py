import dynbem
from dynbem import bem, polar, rotor_definition, rotor_state


def test_backcompat_exports_are_not_public():
    assert not hasattr(dynbem, "AirfoilPolar")
    assert not hasattr(polar, "AirfoilPolar")
    assert not hasattr(dynbem, "BEMModel")
    assert not hasattr(bem, "BEMModel")
    assert not hasattr(rotor_state, "RotorState")
    assert not hasattr(dynbem, "KamanFlap")
    assert not hasattr(rotor_definition, "KamanFlap")
