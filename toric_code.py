import numpy as np
import matplotlib.pyplot as plt

from common import pdist, pdir
from mwpm import MWPM

class ToricCode:

    def __init__(self, L):
        self.L = L # lattice size
        self.stabs = np.zeros((L,L)).astype(np.int8) # Z-stabilizer measurements (meas errs)
        self.qubits = np.zeros((L,L,2)).astype(np.int8) # unit cell: 0: left, 1: top

    def update_stabs(self):
        """Parity of neighboring qubits"""
        for i in range(self.L):
            for j in range(self.L):
                l = self.qubits[i,j,0]
                u = self.qubits[i,j,1]
                r = self.qubits[i,(j+1)%self.L,0]
                d = self.qubits[(i+1)%self.L,j,1]
                self.stabs[i,j] = (l+r+u+d) % 2

    def apply_errors(self, p):
        """Apply errors to qubits"""
        self.qubits ^= np.random.binomial(n=1, p=p, size=self.qubits.shape).astype(np.int8)

    def step(self, p):
        """public interface"""
        self.apply_errors(p)
        self.update_stabs()

    def __str__(self):
        out = ""
        n_spaces_in_tab = '\t'.expandtabs().count(' ')

        for i in range(self.L):
            for j in range(self.L):
                out += f"\t{'X' if self.qubits[i,j,1] else '-'}\t"
            out += "\n\n"

            for j in range(self.L):
                out += f"{'X' if self.qubits[i,j,0] else '-'}\t\033[1m{self.stabs[i,j]}\033[0m\t"

            out += f"({'X' if self.qubits[i,0,0] else '-'})" # PBC wrap
            out += "\n\n"

        for i in range(self.L):
            out += f"{' ' * (n_spaces_in_tab - 1)}({'X' if self.qubits[0,i,1] else '-'})\t" # PBC wrap

        return out

    def join(self, a, b):
        """Join points `a` and `b` via a minimal distance X-error
        correction chain from `a` to `b`."""
        dists = [pdist(x1,x2,self.L) for x1,x2 in zip(a,b)]
        dirs =  [pdir(x1,x2,self.L) for x1,x2 in zip(a,b)]

        # 1. row correction via top qubits in unit cell
        c0 = a[1] # const.
        r0 = a[0] if dirs[0] == -1 else a[0]+1 # dir -1: up, +1: down
        for i in range(dists[0]):
            r = r0 + i * dirs[0]
            self.qubits[r%self.L,c0,1] ^= 1

        # 2. col correction via left qubits in unit cell
        c0 = a[1] if dirs[1] == -1 else a[1]+1 # dir -1: left, +1: right
        r0 = b[0] # const.
        for i in range(dists[1]):
            c = c0 + i * dirs[1]
            self.qubits[r0,c%self.L,0] ^= 1

    def check_log_err(self):
        """Logical error occurred when odd number of crossings of
        either left or top boundary. Assumes after MWPM correction
        only either trivial loops or logical operators left."""
        return self.qubits[0,:,1].sum() % 2 == 1 or self.qubits[:,0,0].sum() % 2 == 1

if __name__ == '__main__':

    # test threshold

    Ls = [4,8,12]
    ps = np.linspace(0.08, 0.12, 5)
    N = 5000

    counts = np.zeros((len(Ls),len(ps)))
    for i,L in enumerate(Ls):
        print(f"- L={L}")
        for j,p in enumerate(ps):
            print(f"-- p={p}")
            for n in range(N):
                t = ToricCode(L)
                t.step(p)
                matchings = MWPM.decode(t.stabs)
                for a,b in matchings:
                    t.join(a,b)
                counts[i,j] += int(t.check_log_err())

    ys = counts / N # mean
    ys_err = np.sqrt( ys * (1-ys) / N ) # MC std (Wald) of Bernoulli r.v.

    for y,yerr in zip(ys,ys_err):
        plt.errorbar(ps,y,yerr=yerr)

    plt.legend(Ls)
    plt.yscale('log')
    plt.xlabel(r'phys. error rate $p$')
    plt.ylabel(r'log. error rate $p_L$')
    plt.title(f'Toric code (L={L}), N={N} samples per point')
    plt.savefig('./thresholds.png', dpi=300)

# if __name__ == '__main__':
#     t = ToricCode(3)
#     t.step(p=0.2)
#     print(t)
#     print('\n')
#     matchings = MWPM.decode(t.stabs)
#     for a,b in matchings:
#         t.join(a,b)
#     t.update_stabs()
#     print(t)
#     print(t.check_log_err())
