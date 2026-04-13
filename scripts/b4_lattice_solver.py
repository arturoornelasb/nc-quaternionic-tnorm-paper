#!/usr/bin/env python3
"""
B^4 Trilemma: Lattice-integrated solver.

Key finding from b4_trilemma_solver.py: monotone partial orders EXIST
on truncated domains, but none are lattices. This script integrates
the lattice constraint directly into backtracking to determine whether
the failure is inherent or an artifact of truncation.

Strategy:
1. Work on 7-element domain {0, 1, i/2, j/2, k/4, -k/4, -1/4}
2. After each assignment, check if every pair still has a POSSIBLE
   meet and join. If not -> contradiction.
3. Print the found order (if any) or prove UNSAT.

Author: J. Arturo Ornelas Brand
Date: 2026-04-12
"""
from fractions import Fraction
import time

# ---- Quaternion arithmetic ----
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
        if r: parts.append(f"{float(r):.4g}")
        if i: parts.append(f"{float(i):.4g}i")
        if j: parts.append(f"{float(j):.4g}j")
        if k: parts.append(f"{float(k):.4g}k")
        return '(' + '+'.join(parts) + ')' if parts else '0'

ZERO = Q(0); ONE = Q(1)

# ---- Solver with lattice integration ----

class LatticeSolver:
    def __init__(self, elts, prod):
        self.n = len(elts)
        self.elts = elts
        self.prod = prod  # prod[i][j] = index or -1
        n = self.n
        # le[i][j]: True/False/None
        self.le = [[None]*n for _ in range(n)]
        for i in range(n): self.le[i][i] = True
        self.trail = []
        self.conflict = False
        self._q = []
        # Precompute product preimages for contrapositive
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
        # Find 0 and 1
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
        """Check if every pair COULD have a meet and join.
        Returns True if feasible, False if some pair has no possible meet/join."""
        n = self.n
        for x in range(n):
            for y in range(x+1, n):
                # Meet: find z such that z<=x, z<=y possible,
                # and for all w with w<=x, w<=y possible: w<=z possible
                has_meet_candidate = False
                has_join_candidate = False
                for z in range(n):
                    # Check if z can be meet(x,y)
                    if self.le[z][x] is not False and self.le[z][y] is not False:
                        # z is a possible lower bound. Is it a possible GLB?
                        # For all w that MUST be a lower bound (le[w][x]=True AND le[w][y]=True),
                        # we need le[w][z] not False
                        is_glb = True
                        for w in range(n):
                            if w != z and self.le[w][x] is True and self.le[w][y] is True:
                                if self.le[w][z] is False:
                                    is_glb = False
                                    break
                        if is_glb:
                            has_meet_candidate = True
                            break
                    # Check if z can be join(x,y)
                if not has_meet_candidate:
                    return False, ('meet', x, y)
                for z in range(n):
                    if self.le[x][z] is not False and self.le[y][z] is not False:
                        is_lub = True
                        for w in range(n):
                            if w != z and self.le[x][w] is True and self.le[y][w] is True:
                                if self.le[z][w] is False:
                                    is_lub = False
                                    break
                        if is_lub:
                            has_join_candidate = True
                            break
                if not has_join_candidate:
                    return False, ('join', x, y)
        return True, None

    def check_lattice_strict(self):
        """Strict lattice check: every pair has a UNIQUE meet and join."""
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

    def count_unknowns(self):
        c = 0
        for i in range(self.n):
            for j in range(self.n):
                if i != j and self.le[i][j] is None: c += 1
        return c

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
        """Backtracking search.
        If require_lattice=True, prune branches where lattice is infeasible."""
        self.nodes = 0
        self.solutions = []

        def search():
            self.nodes += 1
            if self.nodes > max_nodes: return 'TIMEOUT'
            if self.conflict: return 'UNSAT'

            if require_lattice:
                ok, info = self.check_lattice_feasibility()
                if not ok:
                    return 'UNSAT'

            pair = self.pick_unknown()
            if pair is None:
                # All assigned
                if require_lattice:
                    fails = self.check_lattice_strict()
                    if fails:
                        return 'UNSAT'
                return 'SAT'

            i, j = pair
            # Try 3 options: i<j, j<i, incomparable
            for (iv, jv) in [(True, False), (False, True), (False, False)]:
                cp = self.checkpoint()
                if iv: # i <= j
                    self._set(i, j, True)
                elif jv: # j <= i
                    self._set(j, i, True)
                else: # incomparable
                    self._set(i, j, False)
                    if not self.conflict:
                        self._set(j, i, False)
                self.propagate()
                if not self.conflict:
                    result = search()
                    if result == 'SAT':
                        return 'SAT'
                    if result == 'TIMEOUT':
                        return 'TIMEOUT'
                self.restore(cp)
            return 'UNSAT'

        return search()

    def print_order(self):
        n = self.n
        print("\n  Hasse diagram (direct covers):")
        covers = []
        for i in range(n):
            for j in range(n):
                if i != j and self.le[i][j] is True:
                    # Check if it's a cover (no k with i < k < j)
                    is_cover = True
                    for k in range(n):
                        if k != i and k != j and self.le[i][k] is True and self.le[k][j] is True:
                            is_cover = False; break
                    if is_cover:
                        covers.append((i, j))
        for (i, j) in covers:
            print(f"    {self.elts[i]} < {self.elts[j]}")
        # Print levels
        print("\n  Order levels (bottom to top):")
        remaining = set(range(n))
        level = 0
        while remaining:
            # Find minimal elements in remaining
            minimals = []
            for x in remaining:
                is_min = True
                for y in remaining:
                    if y != x and self.le[y][x] is True:
                        is_min = False; break
                if is_min:
                    minimals.append(x)
            print(f"    Level {level}: {[str(self.elts[m]) for m in minimals]}")
            remaining -= set(minimals)
            level += 1
        # Print incomparable pairs
        print("\n  Incomparable pairs:")
        inc = []
        for i in range(n):
            for j in range(i+1, n):
                if self.le[i][j] is False and self.le[j][i] is False:
                    inc.append((i, j))
        for (i, j) in inc:
            print(f"    {self.elts[i]} || {self.elts[j]}")

# ---- Main ----

def main():
    print("=" * 60)
    print("B^4 TRILEMMA: LATTICE-INTEGRATED SOLVER")
    print("=" * 60)

    # Domain: {0, 1, i/2, j/2, k/4, -k/4, -1/4}
    elts = [
        ZERO,
        Q(Fraction(-1,4)),           # -1/4
        Q(0,0,0,Fraction(-1,4)),     # -k/4
        Q(0,0,0,Fraction(1,4)),      # k/4
        Q(0,0,Fraction(1,2)),        # j/2
        Q(0,Fraction(1,2)),          # i/2
        ONE,
    ]
    names = ['0', '-1/4', '-k/4', 'k/4', 'j/2', 'i/2', '1']
    n = len(elts)
    idx = {q: i for i, q in enumerate(elts)}

    print(f"\nDomain: {n} elements")
    for i, (e, nm) in enumerate(zip(elts, names)):
        print(f"  [{i}] {nm}  =  {e}")

    # Build product table
    prod = [[-1]*n for _ in range(n)]
    print("\nProduct table (indices, -1 = outside D):")
    for i in range(n):
        row = []
        for j in range(n):
            p = elts[i] * elts[j]
            if p in idx:
                prod[i][j] = idx[p]
            row.append(f"{prod[i][j]:3d}")
        print(f"  [{names[i]:>5s}] " + " ".join(row))

    in_count = sum(1 for i in range(n) for j in range(n) if prod[i][j] >= 0)
    print(f"\n  {in_count}/{n*n} products land in D ({100*in_count/n/n:.0f}%)")

    # ---- Phase 1: Monotone order WITHOUT lattice ----
    print("\n" + "=" * 60)
    print("PHASE 1: Monotone partial order (no lattice requirement)")
    print("=" * 60)
    solver = LatticeSolver(elts, prod)
    solver.init_boundary()
    if solver.conflict:
        print("CONTRADICTION from boundary alone!")
        return

    unk = solver.count_unknowns()
    print(f"After propagation: {n*(n-1)-unk}/{n*(n-1)} determined, {unk} unknown")

    # Print forced relations
    print("\nForced relations (after boundary + propagation):")
    for i in range(n):
        for j in range(n):
            if i != j and solver.le[i][j] is True:
                print(f"  {names[i]} <= {names[j]}")
    print("\nForced non-relations:")
    for i in range(n):
        for j in range(n):
            if i != j and solver.le[i][j] is False:
                print(f"  NOT {names[i]} <= {names[j]}")
    print(f"\nUnknown pairs:")
    for i in range(n):
        for j in range(n):
            if i != j and solver.le[i][j] is None:
                print(f"  {names[i]} ? {names[j]}")

    result = solver.solve(require_lattice=False)
    print(f"\nResult: {result} ({solver.nodes} nodes)")
    if result == 'SAT':
        solver.print_order()
        fails = solver.check_lattice_strict()
        if fails:
            print(f"\n  Lattice check: FAIL ({len(fails)} failures)")
            for f in fails:
                op, x, y = f
                print(f"    {op} missing: ({names[x]}, {names[y]})")
        else:
            print("\n  Lattice check: PASS")

    # ---- Phase 2: Monotone LATTICE order ----
    print("\n" + "=" * 60)
    print("PHASE 2: Monotone LATTICE order (lattice required)")
    print("=" * 60)
    solver2 = LatticeSolver(elts, prod)
    solver2.init_boundary()
    if solver2.conflict:
        print("CONTRADICTION from boundary alone!")
        return
    result2 = solver2.solve(require_lattice=True)
    print(f"\nResult: {result2} ({solver2.nodes} nodes)")
    if result2 == 'SAT':
        solver2.print_order()
        print("\n  *** MONOTONE LATTICE ORDER EXISTS on 7-element domain! ***")
    elif result2 == 'UNSAT':
        print("\n  *** NO monotone lattice order on 7-element domain. ***")
        print("  (But meets/joins might exist in larger domain.)")

    # ---- Phase 3: Larger domains ----
    print("\n" + "=" * 60)
    print("PHASE 3: 15-element domain with lattice constraint")
    print("=" * 60)

    # Generate 15-element closure
    a = Q(0, Fraction(1,2)); b_gen = Q(0, 0, Fraction(1,2))
    D = {ZERO, ONE, a, b_gen}
    for _ in range(3):  # rounds of closure
        new = set()
        for x in list(D):
            for y in list(D):
                p = x * y
                if p not in D and p not in new and p.norm_sq() <= 1:
                    new.add(p)
        if not new or len(D) + len(new) > 20: break
        D |= new

    elts15 = sorted(D, key=lambda q: (q.norm_sq(), q.c))
    n15 = len(elts15)
    idx15 = {q: i for i, q in enumerate(elts15)}
    names15 = [str(e) for e in elts15]

    print(f"Domain: {n15} elements")
    prod15 = [[-1]*n15 for _ in range(n15)]
    for i in range(n15):
        for j in range(n15):
            p = elts15[i] * elts15[j]
            if p in idx15: prod15[i][j] = idx15[p]
    in15 = sum(1 for i in range(n15) for j in range(n15) if prod15[i][j] >= 0)
    print(f"Products in D: {in15}/{n15*n15} ({100*in15/n15/n15:.0f}%)")

    solver3 = LatticeSolver(elts15, prod15)
    solver3.init_boundary()
    if solver3.conflict:
        print("CONTRADICTION from boundary!")
        return
    unk3 = solver3.count_unknowns()
    print(f"After propagation: {unk3} unknowns")

    t0 = time.time()
    result3 = solver3.solve(require_lattice=True, max_nodes=5000000)
    elapsed = time.time() - t0
    print(f"\nResult: {result3} ({solver3.nodes} nodes, {elapsed:.1f}s)")
    if result3 == 'SAT':
        solver3.print_order()
        print("\n  *** MONOTONE LATTICE ORDER EXISTS on 15-element domain! ***")
    elif result3 == 'UNSAT':
        print("\n  *** NO monotone lattice order on 15-element domain. ***")

    # ---- Summary ----
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  7-element monotone order: SAT (lattice: {result2})")
    print(f"  15-element monotone lattice: {result3}")
    if result2 == 'UNSAT' and result3 == 'UNSAT':
        print("\n  Lattice requirement is the key obstruction.")
        print("  But NOTE: meets/joins might exist outside truncation.")
    print("  The B^4 conjecture status depends on whether the full")
    print("  infinite monoid can provide the missing meets/joins.")

if __name__ == '__main__':
    main()
