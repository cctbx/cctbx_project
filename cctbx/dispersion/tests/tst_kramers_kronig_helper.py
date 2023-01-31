"""Tests for kramkron core_functions"""

import sys

import pytest
import numpy as np

sys.path.append("..")
import kramers_kronig.kramers_kronig_helper as kramers_kronig_helper

@pytest.fixture
def Fe3():
    path = kramers_kronig_helper.SAMPLE_DATA_PATH + "/pf-rd-ox_fftkk.out"
    return(kramers_kronig_helper.parse_data(path))


@pytest.fixture
def Fe2():
    path = kramers_kronig_helper.SAMPLE_DATA_PATH + "/pf-rd-red_fftkk.out"
    return(kramers_kronig_helper.parse_data(path))



def test_parse_Fe3_beginning(Fe3):
    """Test that input is parsed properly."""
    np.testing.assert_equal(Fe3[:4], np.array([[1006.0,-1.95442043311, 5.92009170594],
                                               [1007.0, -2.81223119888, 8.52503033764],
                                               [1008.0, -3.33235759037, 10.1156545543],
                                               [1009.0, -3.56113395273, 10.8295300422]]))

    
def test_parse_Fe3_end(Fe3):
    """Test that input is parsed properly."""
    np.testing.assert_equal(Fe3[-4:], np.array([[24902.0,0.244888423888, 0.418575531827],
                                                [24903.0,0.237396106135, 0.405759333696],
                                                [24904.0,0.220802089543, 0.377388435958],
                                                [24905.0,0.18454034239, 0.315405661584]]))
