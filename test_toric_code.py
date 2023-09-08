from toric_code import ToricCode
import numpy as np

def test_apply_errors():
    np.random.seed(33)
    t = ToricCode(L=3)
    t.apply_errors(p=0.1)
    assert t.qubits[1,0,1] and t.qubits[1,2,0] and t.qubits[2,1,1]
