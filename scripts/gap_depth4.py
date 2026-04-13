#!/usr/bin/env python3
"""Quick depth-4 solver for unequal-norm case."""
import sys
from fractions import Fraction
import time

class Q:
    __slots__ = ('c', '_h')
    def __init__(self, r, i=0, j=0, k=0):
        self.c = (Fraction(r), Fraction(i), Fraction(j), Fraction(k))
        self._h = hash(self.c)
    def __mul__(self, o):
        a,b,c,d = self.c; e,f,g,h = o.c
        return Q(a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g,
                 a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e)
    def norm_sq(self): return sum(x*x for x in self.c)
    def __eq__(self, o): return isinstance(o, Q) and self.c == o.c
    def __hash__(self): return self._h
    def __repr__(self):
        r,i,j,k = self.c
        parts = []
        if r: parts.append(f'{float(r):.4g}')
        if i: parts.append(f'{float(i):.4g}i')
        if j: parts.append(f'{float(j):.4g}j')
        if k: parts.append(f'{float(k):.4g}k')
        return '(' + '+'.join(parts) + ')' if parts else '0'

ZERO = Q(0); ONE = Q(1)

class FastSolver:
    """Solver with propagation, monotonicity, and lattice feasibility."""
    def __init__(self, elts, prod):
        self.n = len(elts)
        self.elts = elts
        self.prod = prod
        n = self.n
        self.le = [[None]*n for _ in range(n)]
        for i in range(n): self.le[i][i] = True
        self.trail = []
        self.conflict = False
        self._q = []
        # Precompute contrapositives
        self.left_pre = [[[] for _ in range(n)] for _ in range(n)]
        self.right_pre = [[[] for _ in range(n)] for _ in range(n)]
        for z in range(n):
            for x in range(n):
                zx = prod[z][x]
                if zx < 0: continue
                for y in range(n):
                    zy = prod[z][y]
                    if zy < 0 or zx == zy: continue
                    self.left_pre[zx][zy].append((x, y))
            for x in range(n):
                xz = prod[x][z]
                if xz < 0: continue
                for y in range(n):
                    yz = prod[y][z]
                    if yz < 0 or xz == yz: continue
                    self.right_pre[xz][yz].append((x, y))
        self.idx0 = next(i for i,e in enumerate(elts) if e == ZERO)
        self.idx1 = next(i for i,e in enumerate(elts) if e == ONE)

    def _set(self, i, j, val):
        if i == j:
            if val is False: self.conflict = True
            return
        cur = self.le[i][j]
        if cur is not None:
            if cur != val: self.conflict = True
            return
        self.trail.append((i, j))
        self.le[i][j] = val
        self._q.append((i, j, val))

    def propagate(self):
        n = self.n; prod = self.prod
        while self._q and not self.conflict:
            i, j, val = self._q.pop(0)
            if val:
                self._set(j, i, False)
                for k in range(n):
                    if self.le[j][k] is True: self._set(i, k, True)
                    if self.le[k][i] is True: self._set(k, j, True)
                for z in range(n):
                    zi,zj = prod[z][i],prod[z][j]
                    if zi >= 0 and zj >= 0: self._set(zi, zj, True)
                    iz,jz = prod[i][z],prod[j][z]
                    if iz >= 0 and jz >= 0: self._set(iz, jz, True)
            else:
                for k in range(n):
                    if k != i and k != j:
                        if self.le[i][k] is True: self._set(k, j, False)
                        if self.le[k][j] is True: self._set(i, k, False)
                for (x, y) in self.left_pre[i][j]: self._set(x, y, False)
                for (x, y) in self.right_pre[i][j]: self._set(x, y, False)

    def check_lattice_feasibility(self):
        n = self.n
        for x in range(n):
            for y in range(x+1, n):
                # Check join feasibility
                has_join = False
                for z in range(n):
                    if self.le[x][z] is not False and self.le[y][z] is not False:
                        is_lub = True
                        for w in range(n):
                            if w != z and self.le[x][w] is True and self.le[y][w] is True:
                                if self.le[z][w] is False:
                                    is_lub = False; break
                        if is_lub: has_join = True; break
                if not has_join: return False
                # Check meet feasibility
                has_meet = False
                for z in range(n):
                    if self.le[z][x] is not False and self.le[z][y] is not False:
                        is_glb = True
                        for w in range(n):
                            if w != z and self.le[w][x] is True and self.le[w][y] is True:
                                if self.le[w][z] is False:
                                    is_glb = False; break
                        if is_glb: has_meet = True; break
                if not has_meet: return False
        return True

    def check_lattice_strict(self):
        n = self.n
        for x in range(n):
            for y in range(x+1, n):
                lower = [z for z in range(n) if self.le[z][x] is True and self.le[z][y] is True]
                if not any(all(self.le[w][z] is True for w in lower) for z in lower):
                    return False
                upper = [z for z in range(n) if self.le[x][z] is True and self.le[y][z] is True]
                if not any(all(self.le[z][w] is True for w in upper) for z in upper):
                    return False
        return True

    def init_boundary(self):
        for x in range(self.n):
            self._set(self.idx0, x, True)
            self._set(x, self.idx1, True)
        self.propagate()

    def checkpoint(self): return len(self.trail)
    def restore(self, cp):
        while len(self.trail) > cp:
            i, j = self.trail.pop()
            self.le[i][j] = None
        self.conflict = False; self._q.clear()

    def pick_unknown(self):
        best = None; best_score = -1
        for i in range(self.n):
            for j in range(i+1, self.n):
                if self.le[i][j] is None:
                    score = sum(1 for k in range(self.n) if self.le[i][k] is not None) + \
                            sum(1 for k in range(self.n) if self.le[j][k] is not None)
                    if score > best_score: best_score = score; best = (i, j)
        return best

    def solve(self, max_nodes=2000000):
        self.nodes = 0
        self.last_print = time.time()
        def search():
            self.nodes += 1
            if self.nodes > max_nodes: return 'TIMEOUT'
            if self.conflict: return 'UNSAT'
            if not self.check_lattice_feasibility(): return 'UNSAT'
            now = time.time()
            if now - self.last_print > 5:
                unk = sum(1 for i in range(self.n) for j in range(i+1,self.n) if self.le[i][j] is None)
                print(f'  ... {self.nodes} nodes, {unk} unknowns remaining', flush=True)
                self.last_print = now
            pair = self.pick_unknown()
            if pair is None:
                if self.check_lattice_strict(): return 'SAT'
                return 'UNSAT'
            i, j = pair
            for (iv, jv) in [(True, False), (False, True), (False, False)]:
                cp = self.checkpoint()
                if iv: self._set(i, j, True)
                elif jv: self._set(j, i, True)
                else:
                    self._set(i, j, False)
                    if not self.conflict: self._set(j, i, False)
                self.propagate()
                if not self.conflict:
                    result = search()
                    if result in ('SAT', 'TIMEOUT'): return result
                self.restore(cp)
            return 'UNSAT'
        return search()


a = Q(0, Fraction(1,2)); b = Q(0, 0, Fraction(1,3))

for max_d in [3, 4, 5, 6]:
    D = {ZERO, ONE, a, b}
    current = {a, b}
    for depth in range(2, max_d + 1):
        nxt = set()
        for w in current:
            for g in [a, b]:
                for p in [w*g, g*w]:
                    if p not in D: nxt.add(p); D.add(p)
        current = nxt

    elts = sorted(D, key=lambda q: (q.norm_sq(), q.c))
    n = len(elts)
    idx = {q: i for i, q in enumerate(elts)}
    prod = [[-1]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            p = elts[i]*elts[j]
            if p in idx: prod[i][j] = idx[p]

    print(f'\n=== Depth {max_d}: {n} elements ===', flush=True)
    solver = FastSolver(elts, prod)
    solver.init_boundary()
    if solver.conflict:
        print(f'  UNSAT from boundary propagation', flush=True)
        continue

    t0 = time.time()
    result = solver.solve(max_nodes=2000000)
    elapsed = time.time() - t0
    print(f'  Result: {result} ({solver.nodes} nodes, {elapsed:.1f}s)', flush=True)
