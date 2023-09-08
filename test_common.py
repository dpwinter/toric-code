from toric_code.common import *

def test_pdist_wrap():
    assert pdist(1,5,5)==1

def test_pdist_nowrap():
    assert pdist(3,6,8)==3

def test_pdist_same():
    assert pdist(2,4,4)==2

def test_pdir_wrap():
    assert pdir(1,5,5)==-1

def test_pdir_nowrap():
    assert pdir(3,6,8)==1

def test_pdir_same():
    assert pdir(2,4,4)==1

def test_pmanhatten_rowwrap():
    assert pmanhatten((1,2), (5,2), 5) == 1

def test_pmanhatten_colwrap():
    assert pmanhatten((2,1), (2,5), 5) == 1

def test_pmanhatten_rowcolwrap():
    assert pmanhatten((1,1), (5,5), 5) == 2

def test_pmanhatten_nowrap():
    assert pmanhatten((2,2), (3,3), 5) == 2
