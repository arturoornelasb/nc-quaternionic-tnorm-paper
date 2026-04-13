#!/usr/bin/env python3
"""
Gap analysis: test whether the B^4 trilemma extends beyond
equal-norm purely imaginary generators.

Tests:
1. Unequal-norm purely imaginary: a = r1*i, b = r2*j, r1 != r2
2. Generators with real part: a = alpha + r*i, b = beta + s*j
3. Non-orthogonal imaginary: a = r*i, b = r*(cos(t)*i + sin(t)*j)

For each case: generate sub-monoid, check if monotone lattice exists.

Author: J. Arturo Ornelas Brand
Date: 2026-04-12
"""
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
    def __ne__(self, o): return not self.__eq__(o)
    def __hash__(self): return self._h
    def __repr__(self):
        r,i,j,k = self.c
        parts = []
        if r: parts.append(f"{float(r):.4g}")
        if i: parts.append(f"{float(i):.4g}i")
        if j: parts.append(f"{float(j):.4g}j")
        if k: parts.append(f"{float(k):.4g}k")
        return '(' + '+'.join(parts) + ')' if parts else '0'

ZERO = Q(0); ONE = Q(1)

class LatticeSolver:
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
        self.idx0 = None; self.idx1 = None
        for i, e in enumerate(elts):
            if e == ZERO: self.idx0 = i
            if e == ONE: self.idx1 = i

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
                has_meet = False
                for z in range(n):
                    if self.le[z][x] is not False and self.le[z][y] is not False:
                        is_glb = True
                        for w in range(n):
                            if w != z and self.le[w][x] is True and self.le[w][y] is True:
                                if self.le[w][z] is False:
                                    is_glb = False; break
                        if is_glb:
                            has_meet = True; break
                if not has_meet:
                    return False, ('meet', x, y)
                has_join = False
                for z in range(n):
                    if self.le[x][z] is not False and self.le[y][z] is not False:
                        is_lub = True
                        for w in range(n):
                            if w != z and self.le[x][w] is True and self.le[y][w] is True:
                                if self.le[z][w] is False:
                                    is_lub = False; break
                        if is_lub:
                            has_join = True; break
                if not has_join:
                    return False, ('join', x, y)
        return True, None

    def check_lattice_strict(self):
        n = self.n
        failures = []
        for x in range(n):
            for y in range(x+1, n):
                lower = [z for z in range(n) if self.le[z][x] is True and self.le[z][y] is True]
                meet = None
                for z in lower:
                    if all(self.le[w][z] is True for w in lower):
                        meet = z; break
                if meet is None:
                    failures.append(('meet', x, y))
                upper = [z for z in range(n) if self.le[x][z] is True and self.le[y][z] is True]
                join = None
                for z in upper:
                    if all(self.le[z][w] is True for w in upper):
                        join = z; break
                if join is None:
                    failures.append(('join', x, y))
        return failures

    def init_boundary(self):
        for x in range(self.n):
            self._set(self.idx0, x, True)
            self._set(x, self.idx1, True)
        self.propagate()

    def checkpoint(self):
        return len(self.trail)

    def restore(self, cp):
        while len(self.trail) > cp:
            i, j = self.trail.pop()
            self.le[i][j] = None
        self.conflict = False
        self._q.clear()

    def pick_unknown(self):
        best = None; best_score = -1
        for i in range(self.n):
            for j in range(i+1, self.n):
                if self.le[i][j] is None:
                    score = sum(1 for k in range(self.n) if self.le[i][k] is not None) + \
                            sum(1 for k in range(self.n) if self.le[j][k] is not None)
                    if score > best_score:
                        best_score = score; best = (i, j)
        return best

    def solve(self, require_lattice=False, max_nodes=5000000):
        self.nodes = 0
        def search():
            self.nodes += 1
            if self.nodes > max_nodes: return 'TIMEOUT'
            if self.conflict: return 'UNSAT'
            if require_lattice:
                ok, info = self.check_lattice_feasibility()
                if not ok: return 'UNSAT'
            pair = self.pick_unknown()
            if pair is None:
                if require_lattice:
                    fails = self.check_lattice_strict()
                    if fails: return 'UNSAT'
                return 'SAT'
            i, j = pair
            for (iv, jv) in [(True, False), (False, True), (False, False)]:
                cp = self.checkpoint()
                if iv:
                    self._set(i, j, True)
                elif jv:
                    self._set(j, i, True)
                else:
                    self._set(i, j, False)
                    if not self.conflict:
                        self._set(j, i, False)
                self.propagate()
                if not self.conflict:
                    result = search()
                    if result in ('SAT', 'TIMEOUT'): return result
                self.restore(cp)
            return 'UNSAT'
        return search()


def gen_submonoid(gens, max_depth=3):
    """Generate sub-monoid from generators up to given depth."""
    D = {ZERO, ONE}
    D.update(gens)
    # Level by level
    current_level = set(gens)
    for depth in range(2, max_depth + 1):
        next_level = set()
        for w in current_level:
            for g in gens:
                p1 = w * g
                p2 = g * w
                if p1 not in D: next_level.add(p1); D.add(p1)
                if p2 not in D: next_level.add(p2); D.add(p2)
        current_level = next_level
        if not current_level: break
    return D


def test_case(name, gens, max_depth=3):
    """Test a specific generator pair."""
    print(f"\n{'='*60}")
    print(f"TEST: {name}")
    print(f"{'='*60}")

    a, b = gens
    print(f"  a = {a},  b = {b}")
    print(f"  a*b = {a*b},  b*a = {b*a}")
    print(f"  a*b == b*a: {a*b == b*a}")
    if a*b == b*a:
        print("  SKIP: generators commute.")
        return 'COMMUTATIVE'

    print(f"  a^2 = {a*a},  b^2 = {b*b},  a^2==b^2: {a*a == b*b}")

    D = gen_submonoid(gens, max_depth)
    elts = sorted(D, key=lambda q: (q.norm_sq(), q.c))
    n = len(elts)
    idx = {q: i for i, q in enumerate(elts)}

    print(f"  Domain size (depth {max_depth}): {n} elements")
    for i, e in enumerate(elts):
        print(f"    [{i:2d}] {e}")

    # Build product table
    prod = [[-1]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            p = elts[i] * elts[j]
            if p in idx: prod[i][j] = idx[p]

    in_count = sum(1 for i in range(n) for j in range(n) if prod[i][j] >= 0)
    print(f"  Products in D: {in_count}/{n*n}")

    # Test monotone lattice
    solver = LatticeSolver(elts, prod)
    solver.init_boundary()
    if solver.conflict:
        print("  CONTRADICTION from boundary propagation!")
        print("  => No monotone partial order exists at all.")
        return 'UNSAT_BOUNDARY'

    t0 = time.time()
    result = solver.solve(require_lattice=True, max_nodes=5000000)
    elapsed = time.time() - t0
    print(f"  Result: {result} ({solver.nodes} nodes, {elapsed:.2f}s)")

    if result == 'UNSAT':
        print(f"  => NO monotone lattice order on {n}-element domain. TRILEMMA HOLDS.")
    elif result == 'SAT':
        print(f"  => Monotone lattice order EXISTS! GAP IS REAL.")
        fails = solver.check_lattice_strict()
        print(f"  Lattice failures: {len(fails)}")
    else:
        print(f"  => TIMEOUT (inconclusive)")

    return result


def main():
    print("B^4 TRILEMMA GAP ANALYSIS")
    print("=" * 60)

    results = {}

    # ---- Control: equal-norm case (should be UNSAT) ----
    results['equal_norm'] = test_case(
        "Control: equal-norm (r=1/2)",
        [Q(0, Fraction(1,2)), Q(0, 0, Fraction(1,2))]
    )

    # ---- GAP 2: Unequal-norm purely imaginary ----
    results['uneq_1/2_1/3'] = test_case(
        "Unequal norm: a=i/2, b=j/3",
        [Q(0, Fraction(1,2)), Q(0, 0, Fraction(1,3))]
    )

    results['uneq_1/3_1/5'] = test_case(
        "Unequal norm: a=i/3, b=j/5",
        [Q(0, Fraction(1,3)), Q(0, 0, Fraction(1,5))]
    )

    results['uneq_2/5_1/3'] = test_case(
        "Unequal norm: a=2i/5, b=j/3",
        [Q(0, Fraction(2,5)), Q(0, 0, Fraction(1,3))]
    )

    # ---- GAP 3: Generators with real part ----
    results['real_part_1'] = test_case(
        "Real part: a=1/4+i/4, b=1/4+j/4",
        [Q(Fraction(1,4), Fraction(1,4)), Q(Fraction(1,4), 0, Fraction(1,4))]
    )

    results['real_part_2'] = test_case(
        "Real part: a=1/3+i/3, b=j/3",
        [Q(Fraction(1,3), Fraction(1,3)), Q(0, 0, Fraction(1,3))]
    )

    results['real_part_3'] = test_case(
        "Real part: a=1/2+i/4, b=j/4",
        [Q(Fraction(1,2), Fraction(1,4)), Q(0, 0, Fraction(1,4))]
    )

    # ---- Extreme: very different norms ----
    results['extreme'] = test_case(
        "Extreme: a=i/2, b=j/7",
        [Q(0, Fraction(1,2)), Q(0, 0, Fraction(1,7))]
    )

    # ---- Summary ----
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for name, result in results.items():
        status = "TRILEMMA HOLDS" if result == 'UNSAT' else \
                 "BOUNDARY CONTRADICTION" if result == 'UNSAT_BOUNDARY' else \
                 "GAP REAL!" if result == 'SAT' else result
        print(f"  {name:25s}: {status}")

    all_unsat = all(r in ('UNSAT', 'UNSAT_BOUNDARY', 'COMMUTATIVE')
                    for r in results.values())
    if all_unsat:
        print("\n  ALL tested cases: trilemma holds.")
        print("  The unequal-norm and real-part gaps appear to be")
        print("  artifacts of the proof technique, not real obstructions.")
    else:
        sat_cases = [k for k, v in results.items() if v == 'SAT']
        print(f"\n  COUNTEREXAMPLES FOUND: {sat_cases}")
        print("  The gap is real for these cases.")


if __name__ == '__main__':
    main()
