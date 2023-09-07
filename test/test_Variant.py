import pytest
from tbprofiler_parser import Variant

def test_get_position_nucleotide():
    assert Variant.get_position("c.1234A>T") == 1234

def test_get_position_aminoacid():
    assert Variant.get_position("p.Met291Ile") == 291