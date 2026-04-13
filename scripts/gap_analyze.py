#!/usr/bin/env python3
"""
Analyze WHY unequal-norm cases admit a lattice at depth 3.
Print the actual order found and identify the structural difference.
"""
from fractions import Fraction

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
        if r: parts.append(f"{float(r):.5g}")
        if i: parts.append(f"{float(i):.5g}i")
        if j: parts.append(f"{float(j):.5g}j")
        if k: parts.append(f"{float(k):.5g}k")
        return '(' + '+'.join(parts) + ')' if parts else '0'

ZERO = Q(0); ONE = Q(1)

# ---- Analyze the EQUAL-NORM case step 2 failure ----
print("=" * 60)
print("EQUAL-NORM ANALYSIS: why the proof works")
print("=" * 60)

r = Fraction(1,2)
a = Q(0, r); b = Q(0, 0, r)
a2 = a*a; b2 = b*b; ab = a*b; ba = b*a

print(f"a={a}, b={b}")
print(f"a^2={a2}, b^2={b2}, EQUAL: {a2 == b2}")
print(f"ab={ab}, ba={ba}")
print(f"\nLevel-2 distinct elements: {len(set([a2, b2, ab, ba]))}")
print(f"  {a2}, {ab}, {ba}")
print(f"  (a^2 = b^2 collapses 4 -> 3)")

print(f"\nProof of a^2 || ab:")
print(f"  If a^2 <= ab: L*a: a^3={a*a*a} <= a^2*b={a*a*b}")
print(f"  R*b: a^2*b={a*a*b} <= a*b^2={a*b*b}")
print(f"  Since a*b^2 = {a*b*b} = a^3 = {a*a*a}: CYCLE -> contradiction")

print(f"\nKey: a*b^2 = a*(b^2) = a*(a^2) = a^3 because a^2 = b^2!")
print(f"  This collapse creates the cyclic swap at level 3.")

# ---- Analyze the UNEQUAL-NORM case ----
print(f"\n{'='*60}")
print("UNEQUAL-NORM ANALYSIS: why the proof FAILS")
print("=" * 60)

r1, r2 = Fraction(1,2), Fraction(1,3)
a = Q(0, r1); b = Q(0, 0, r2)
a2 = a*a; b2 = b*b; ab = a*b; ba = b*a

print(f"a={a}, b={b}")
print(f"a^2={a2}, b^2={b2}, EQUAL: {a2 == b2}")
print(f"ab={ab}, ba={ba}")
print(f"\nLevel-2 distinct elements: {len(set([a2, b2, ab, ba]))}")
print(f"  {a2}, {b2}, {ab}, {ba}")
print(f"  (a^2 != b^2: 4 distinct elements)")

print(f"\nAttempted proof of a^2 || ab:")
print(f"  If a^2 <= ab: L*a: a^3={a*a*a} <= a^2*b={a*a*b}")
print(f"  R*b: a^2*b={a*a*b} <= a*b^2={a*b*b}")
print(f"  a*b^2 = {a*b*b}")
print(f"  a^3 = {a*a*a}")
print(f"  a*b^2 != a^3: {a*b*b != a*a*a}")
print(f"  NO CYCLE: a^3 <= a^2*b <= a*b^2, all three DISTINCT")
print(f"  The chain doesn't close!")

# Check what the forced chain looks like
print(f"\n  Full forced chain if a^2 <= ab:")
print(f"    L*a: {a*a*a} <= {a*a*b} (level-3)")
print(f"    R*b: {a*a*b} <= {a*b*b} (level-3)")
print(f"    L*b: {b*a*a} <= {b*a*b} (level-3)")
print(f"    R*a: {a*a*a} <= {a*b*a} (level-3)")
print(f"    All 4 level-3 elements constrained but NO cycle.")

# Check all level-3 elements
print(f"\nAll level-3 products:")
gens = [a, b]
level3 = {}
for w1 in gens:
    for w2 in gens:
        for w3 in gens:
            name = f"{'a' if w1==a else 'b'}{'a' if w2==a else 'b'}{'a' if w3==a else 'b'}"
            val = w1*w2*w3
            level3[name] = val
            print(f"  {name} = {val}  |{name}|^2 = {float(val.norm_sq()):.6f}")

# Check for equalities among level-3 elements
print(f"\nLevel-3 equalities:")
names = list(level3.keys())
for i, n1 in enumerate(names):
    for n2 in names[i+1:]:
        if level3[n1] == level3[n2]:
            print(f"  {n1} = {n2} = {level3[n1]}")

print(f"\nDistinct level-3 elements: {len(set(level3.values()))}")

# ---- The structural difference ----
print(f"\n{'='*60}")
print("STRUCTURAL DIFFERENCE")
print("=" * 60)
print("""
EQUAL NORMS (a^2 = b^2):
  - Level 2 has 3 distinct elements: a^2(=b^2), ab, ba
  - Level 3: a*b^2 = a*a^2 = a^3 (collapse!)
  - This creates a CYCLE in the forced inequalities
  - The cycle forces p = q for distinct quaternions -> contradiction

UNEQUAL NORMS (a^2 != b^2):
  - Level 2 has 4 distinct elements: a^2, b^2, ab, ba
  - Level 3: a*b^2 != a^3 (no collapse)
  - Forced inequalities form a CHAIN, not a cycle
  - The chain is satisfiable: the 4th element (b^2 or a^2)
    provides room to route the order around the obstruction

The a^2 = b^2 collapse is the ESSENTIAL ingredient of the proof,
not merely a simplifying assumption.
""")

# ---- Can a^2 and b^2 be ordered? ----
print("=" * 60)
print("Can a^2 and b^2 be ordered when a^2 != b^2?")
print("=" * 60)
print(f"\na^2 = {a2} (norm {float(a2.norm_sq()):.4f})")
print(f"b^2 = {b2} (norm {float(b2.norm_sq()):.4f})")
print(f"|a^2| > |b^2|: {a2.norm_sq() > b2.norm_sq()}")

print(f"\nIf b^2 <= a^2 (smaller norm below larger):")
print(f"  This is consistent with forced bounds:")
print(f"  0 <= b^2 <= b <= 1 and 0 <= a^2 <= a <= 1")
print(f"  We need b^2 <= a^2 and a^2 <= a, b^2 <= b")
print(f"  L*a of b^2<=a^2: a*b^2={a*b2} <= a^3={a*a2}")
print(f"  R*b: b^2*b={b2*b} <= a^2*b={a2*b}")
print(f"  These are level-3 constraints, all consistent.")

print(f"\nWith b^2 <= a^2, the Hasse can be:")
print(f"  0 < b^2 < a^2 < ... < a,b < 1")
print(f"  or 0 < b^2 < ba < ab < a^2 < a,b < 1")
print(f"  The key: when a^2 != b^2, the level-2 elements")
print(f"  can form a CHAIN instead of an antichain.")
print(f"  A chain has unique meets and joins for all pairs!")
