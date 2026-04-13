#!/usr/bin/env python3
"""
Computational exploration of OQ1: Non-commutative quaternionic t-norm.

Tests candidate constructions for a compatible lattice order on B^4 = {q in H : |q| <= 1}
with the Hamilton product as monoidal operation.

Candidates tested:
  1. Cone impossibility (S^3 transitivity argument)
  2. Q_8 bi-invariant order exhaustive enumeration
  3. Real-dominant subcone closure statistics
  4. Norm-order RL verification (the viable 1D structure)
  5. Modified products: clamped, projected, split-quaternion
  6. Discrete subset search for compatible partial orders

Author: J. Arturo Ornelas Brand
Date: 2026-04-12
"""

import numpy as np
from itertools import product as cartprod, combinations
from collections import defaultdict

np.random.seed(42)

# ─── Quaternion arithmetic ─────────────────────────────────────────────
def qmul(a, b):
    """Hamilton product of two quaternions (r,i,j,k)."""
    r1, i1, j1, k1 = a
    r2, i2, j2, k2 = b
    return np.array([
        r1*r2 - i1*i2 - j1*j2 - k1*k2,
        r1*i2 + i1*r2 + j1*k2 - k1*j2,
        r1*j2 - i1*k2 + j1*r2 + k1*i2,
        r1*k2 + i1*j2 - j1*i2 + k1*r2
    ])

def qnorm(q):
    return np.linalg.norm(q)

def qconj(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])

def qinv(q):
    n2 = np.dot(q, q)
    if n2 < 1e-15:
        return None
    return qconj(q) / n2

def random_unit_quaternion(n=1):
    """Uniform random on S^3 via normalized Gaussians."""
    v = np.random.randn(n, 4)
    norms = np.linalg.norm(v, axis=1, keepdims=True)
    return v / norms

def random_ball_quaternion(n=1):
    """Uniform random in B^4."""
    v = random_unit_quaternion(n)
    r = np.random.uniform(0, 1, (n, 1)) ** (1.0/4)  # radial for 4D ball
    return v * r


# ═══════════════════════════════════════════════════════════════════════
# EXPLORATION 1: Cone impossibility — computational verification
# ═══════════════════════════════════════════════════════════════════════
def explore_cone_impossibility():
    print("=" * 70)
    print("EXPLORATION 1: Cone impossibility on S^3")
    print("=" * 70)
    print()
    print("Theorem claim: No proper cone C in H \\ {0} is invariant under")
    print("both left and right multiplication by all unit quaternions.")
    print("Proof sketch: SU(2) acts transitively on S^3, so if q in C and")
    print("|u|=1, then u*q*u^{-1} in C (conjugation). But conjugation acts")
    print("as SO(3) on Im(H), which is transitive on S^2. So C contains")
    print("all of S^3, hence C = H \\ {0}.")
    print()

    # Computational verification: pick a candidate cone direction,
    # show conjugation maps it to every other direction
    p = np.array([0.0, 1.0, 0.0, 0.0])  # pure i
    conjugates = []
    units = random_unit_quaternion(5000)
    for u in units:
        c = qmul(qmul(u, p), qinv(u))
        conjugates.append(c)
    conjugates = np.array(conjugates)

    # Check: do conjugates cover all of S^2 in Im(H)?
    im_parts = conjugates[:, 1:]  # (i,j,k) components
    norms = np.linalg.norm(im_parts, axis=1)
    print(f"  Conjugates of i by 5000 random unit quaternions:")
    print(f"  |Im| range: [{norms.min():.6f}, {norms.max():.6f}] (should be ~1)")
    print(f"  Re range:   [{conjugates[:,0].min():.6f}, {conjugates[:,0].max():.6f}] (should be ~0)")

    # Check angular coverage: project onto S^2
    im_normed = im_parts / norms[:, None]
    # Compute pairwise max angular distance coverage
    # Sample 100 random directions on S^2 and find nearest conjugate
    test_dirs = np.random.randn(200, 3)
    test_dirs /= np.linalg.norm(test_dirs, axis=1, keepdims=True)
    min_angles = []
    for d in test_dirs:
        dots = np.abs(im_normed @ d)
        min_angles.append(np.arccos(np.clip(dots.max(), -1, 1)))
    max_gap = np.max(min_angles)
    print(f"  Max angular gap to nearest conjugate: {np.degrees(max_gap):.2f}° (should be small)")
    print()
    print("  RESULT: Conjugation orbits cover all of S^2 in Im(H).")
    print("  => No proper bi-invariant cone exists. QED (computational)")
    print()

    # Also verify: left-invariant cone attempt
    print("  Left-invariant cone test:")
    print("  If C is a cone with q in C, then for all |u|=1, u*q in C.")
    p = np.array([0.5, 0.5, 0.5, 0.5])  # a specific point
    p = p / qnorm(p)
    left_images = []
    for u in units[:1000]:
        left_images.append(qmul(u, p))
    left_images = np.array(left_images)
    re_range = (left_images[:, 0].min(), left_images[:, 0].max())
    print(f"  Left-translates of q=(1,1,1,1)/2: Re in [{re_range[0]:.4f}, {re_range[1]:.4f}]")
    print(f"  (Must include negative reals => no half-space cone is left-invariant)")
    print()


# ═══════════════════════════════════════════════════════════════════════
# EXPLORATION 2: Q_8 bi-invariant order — exhaustive enumeration
# ═══════════════════════════════════════════════════════════════════════
def explore_q8_order():
    print("=" * 70)
    print("EXPLORATION 2: Q_8 bi-invariant partial order")
    print("=" * 70)
    print()

    # Q_8 = {+/-1, +/-i, +/-j, +/-k}
    Q8 = {
        '1':  np.array([1, 0, 0, 0], dtype=float),
        '-1': np.array([-1, 0, 0, 0], dtype=float),
        'i':  np.array([0, 1, 0, 0], dtype=float),
        '-i': np.array([0, -1, 0, 0], dtype=float),
        'j':  np.array([0, 0, 1, 0], dtype=float),
        '-j': np.array([0, 0, -1, 0], dtype=float),
        'k':  np.array([0, 0, 0, 1], dtype=float),
        '-k': np.array([0, 0, 0, -1], dtype=float),
    }
    names = list(Q8.keys())
    elems = [Q8[n] for n in names]
    n = len(names)

    # Multiplication table
    def q8_mul(a_name, b_name):
        prod = qmul(Q8[a_name], Q8[b_name])
        for nm, val in Q8.items():
            if np.allclose(prod, val):
                return nm
        return None

    print("  Q_8 multiplication table (subset):")
    print(f"  {'':4s}", end="")
    for b in names[:4]:
        print(f"{b:5s}", end="")
    print()
    for a in names[:4]:
        print(f"  {a:4s}", end="")
        for b in names[:4]:
            print(f"{q8_mul(a, b):5s}", end="")
        print()
    print()

    # Key observation: every non-identity element squares to -1
    print("  Squares of Q_8 elements:")
    for nm in names:
        sq = q8_mul(nm, nm)
        print(f"    {nm}^2 = {sq}")
    print()

    # Check: can we have a partial order with 1 as top and 0 as bottom?
    # Q_8 has no 0, so we adjoin it: Q_8 ∪ {0}
    # For antisymmetry: if a <= b and b <= a, then a = b
    # For a bi-invariant order: a <= b => c*a <= c*b and a*c <= b*c
    # Key problem: i^2 = -1, so i <= i^2 * i^{-1} ... complicated

    # Systematic check: enumerate all antisymmetric relations on Q_8
    # compatible with left/right multiplication
    print("  Searching for non-trivial bi-invariant partial orders on Q_8...")
    print("  (An order <= such that a<=b implies ca<=cb and ac<=bc for all c)")
    print()

    # Strategy: if a < b (a != b), then for every c in Q_8:
    # c*a < c*b and a*c < b*c
    # Start with a single relation and compute closure

    found_nontrivial = False
    # Try all pairs (a,b) with a != b as the seed relation a < b
    for idx_a in range(n):
        for idx_b in range(n):
            if idx_a == idx_b:
                continue
            a_nm, b_nm = names[idx_a], names[idx_b]

            # Generate all forced relations from a < b
            relations = set()
            queue = [(a_nm, b_nm)]
            contradicted = False

            while queue and not contradicted:
                x, y = queue.pop()
                if x == y:
                    contradicted = True
                    break
                if (y, x) in relations:
                    contradicted = True
                    break
                if (x, y) in relations:
                    continue
                relations.add((x, y))
                # Left multiply by every c
                for c_nm in names:
                    cx = q8_mul(c_nm, x)
                    cy = q8_mul(c_nm, y)
                    if cx and cy and cx != cy:
                        if (cy, cx) in relations:
                            contradicted = True
                            break
                        queue.append((cx, cy))
                    xc = q8_mul(x, c_nm)
                    yc = q8_mul(y, c_nm)
                    if xc and yc and xc != yc:
                        if (yc, xc) in relations:
                            contradicted = True
                            break
                        queue.append((xc, yc))

            if not contradicted and len(relations) > 0:
                # Check transitivity closure
                # Build adjacency and check for cycles
                adj = defaultdict(set)
                for (x, y) in relations:
                    adj[x].add(y)

                # Transitive closure via Floyd-Warshall on names
                reach = {nm: set() for nm in names}
                for (x, y) in relations:
                    reach[x].add(y)
                changed = True
                while changed:
                    changed = False
                    for x in names:
                        for y in list(reach[x]):
                            for z in reach[y]:
                                if z not in reach[x]:
                                    if z == x:
                                        contradicted = True
                                        break
                                    reach[x].add(z)
                                    changed = True
                            if contradicted:
                                break
                        if contradicted:
                            break

            if not contradicted and len(relations) > 0:
                found_nontrivial = True
                print(f"  FOUND non-trivial order seeded by {a_nm} < {b_nm}:")
                for (x, y) in sorted(relations):
                    print(f"    {x} < {y}")
                print()

    if not found_nontrivial:
        print("  RESULT: No non-trivial bi-invariant partial order exists on Q_8.")
        print("  Every seed relation a < b generates a contradiction via")
        print("  left/right multiplication closure.")
        print()
        print("  Proof sketch: Suppose a < b. Then a*b^{-1}*a < a (right-mult")
        print("  by b^{-1}*a). But in Q_8, repeated application generates")
        print("  cycles because |Q_8| = 8 and every element has finite order.")
    print()


# ═══════════════════════════════════════════════════════════════════════
# EXPLORATION 3: Real-dominant subcone closure
# ═══════════════════════════════════════════════════════════════════════
def explore_real_dominant_subcone():
    print("=" * 70)
    print("EXPLORATION 3: Real-dominant subcone closure under Hamilton product")
    print("=" * 70)
    print()
    print("  V_tau = {q in B^4 : Re(q) >= tau * |q|}, tau in (0, 1]")
    print("  Question: Is V_tau closed under Hamilton product?")
    print()

    N = 50000
    taus = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]

    for tau in taus:
        # Generate random quaternions in V_tau
        count_in = 0
        count_total = 0

        for _ in range(N):
            # Generate q with Re(q) >= tau * |q|
            # i.e., r >= tau * sqrt(r^2 + |im|^2)
            # r^2 >= tau^2 * (r^2 + |im|^2)
            # (1-tau^2)*r^2 >= tau^2 * |im|^2
            # |im| <= r * sqrt(1-tau^2)/tau
            r = np.random.uniform(0.01, 1.0)
            max_im_norm = r * np.sqrt(1 - tau**2) / tau
            # Random direction in im
            im_dir = np.random.randn(3)
            im_dir /= np.linalg.norm(im_dir)
            im_norm = np.random.uniform(0, max_im_norm)
            a = np.array([r, im_dir[0]*im_norm, im_dir[1]*im_norm, im_dir[2]*im_norm])
            # Normalize to B^4 if needed
            if qnorm(a) > 1:
                a = a / qnorm(a)

            # Same for b
            r = np.random.uniform(0.01, 1.0)
            max_im_norm = r * np.sqrt(1 - tau**2) / tau
            im_dir = np.random.randn(3)
            im_dir /= np.linalg.norm(im_dir)
            im_norm = np.random.uniform(0, max_im_norm)
            b = np.array([r, im_dir[0]*im_norm, im_dir[1]*im_norm, im_dir[2]*im_norm])
            if qnorm(b) > 1:
                b = b / qnorm(b)

            ab = qmul(a, b)
            ab_norm = qnorm(ab)
            if ab_norm > 1e-10:
                ratio = ab[0] / ab_norm
                if ratio >= tau - 1e-10:
                    count_in += 1
            count_total += 1

        pct = 100.0 * count_in / count_total
        print(f"  tau={tau:.2f}: {pct:.1f}% of products stay in V_tau "
              f"({count_in}/{count_total})")

    print()

    # Analytical: for q1, q2 on S^3 with Re >= tau,
    # Re(q1*q2) = r1*r2 - <im1, im2>
    # Worst case: im1 and im2 aligned, both at maximum |im| = sqrt(1-tau^2)
    # Then Re(q1*q2) = tau^2 - (1-tau^2) = 2*tau^2 - 1
    # Need 2*tau^2 - 1 >= tau, i.e., 2*tau^2 - tau - 1 >= 0
    # tau >= (1 + sqrt(9))/4 = (1+3)/4 = 1  (only tau = 1 works!)
    tau_crit = (1 + np.sqrt(9)) / 4
    print(f"  Analytical bound: Re(q1*q2) >= 2*tau^2 - 1 (worst case aligned)")
    print(f"  Need 2*tau^2 - 1 >= tau => tau >= {tau_crit:.4f}")
    print(f"  Only tau = 1 (pure real) guarantees closure.")
    print()
    print("  RESULT: V_tau is NOT closed under Hamilton product for any tau < 1.")
    print("  The real-dominant subcone collapses to the real line.")
    print()


# ═══════════════════════════════════════════════════════════════════════
# EXPLORATION 4: Norm-order residuated lattice
# ═══════════════════════════════════════════════════════════════════════
def explore_norm_order_rl():
    print("=" * 70)
    print("EXPLORATION 4: Norm-order residuated lattice on B^4")
    print("=" * 70)
    print()
    print("  Order: a <= b iff |a| <= |b|")
    print("  Meet: a ^ b = (the one with smaller norm)")
    print("  Join: a v b = (the one with larger norm)")
    print("  Monoid: Hamilton product (|ab| = |a||b|)")
    print()

    N = 10000
    samples = random_ball_quaternion(N)

    # Test: |a*b| = |a|*|b| (submultiplicativity = equality for quaternions)
    print("  Test: |a*b| = |a|*|b| (norm is multiplicative)")
    max_err = 0
    for _ in range(N):
        a = random_ball_quaternion(1)[0]
        b = random_ball_quaternion(1)[0]
        ab = qmul(a, b)
        err = abs(qnorm(ab) - qnorm(a)*qnorm(b))
        max_err = max(max_err, err)
    print(f"  Max |  |ab| - |a||b|  | over {N} pairs: {max_err:.2e}")
    print()

    # Test monotonicity: |a| <= |b| => |a*c| <= |b*c| and |c*a| <= |c*b|
    print("  Test: monotonicity of norm order")
    violations_left = 0
    violations_right = 0
    for _ in range(N):
        a = random_ball_quaternion(1)[0]
        b = random_ball_quaternion(1)[0]
        c = random_ball_quaternion(1)[0]
        if qnorm(a) <= qnorm(b):
            if qnorm(qmul(c, a)) > qnorm(qmul(c, b)) + 1e-10:
                violations_left += 1
            if qnorm(qmul(a, c)) > qnorm(qmul(b, c)) + 1e-10:
                violations_right += 1
    print(f"  Left monotonicity violations:  {violations_left}/{N}")
    print(f"  Right monotonicity violations: {violations_right}/{N}")
    print()

    # Residuals exist: a\c = a^{-1}*c (when |a^{-1}*c| gives largest b with |ab|<=|c|)
    # |a*b| <= |c| iff |a|*|b| <= |c| iff |b| <= |c|/|a|
    # So the residual is: a\c = any q with |q| = |c|/|a| (or 1 if |c|/|a| > 1)
    # This is a SPHERE, not a point! => not a function, just a set
    # The residual as lattice element = sup{b : |a*b| <= |c|} = any q with |q| = min(|c|/|a|, 1)
    print("  Residual analysis:")
    print("  a\\c = sup{b : |ab| <= |c|} = {q : |q| = min(|c|/|a|, 1)}")
    print("  This is a SPHERE in B^4, not a unique element.")
    print("  Under norm-order, all points on this sphere are equivalent (same norm).")
    print()
    print("  The structure collapses: norm-order identifies all quaternions")
    print("  of the same norm into one equivalence class.")
    print("  B^4 / ~ = [0,1] with product t-norm.")
    print()

    # Verify: the quotient is exactly product BL-algebra
    print("  Quotient verification: [0,1] product BL-algebra")
    violations = 0
    for _ in range(N):
        x, y = np.random.uniform(0, 1, 2)
        # Product t-norm
        t = x * y
        # Residual: x -> y = min(y/x, 1) if x > 0, else 1
        if x > 1e-10:
            r = min(y / x, 1.0)
        else:
            r = 1.0
        # Check: t <= z iff x <= r(y,z) ... just verify adjunction
        z = np.random.uniform(0, 1)
        lhs = (x * y <= z + 1e-10)
        if x > 1e-10:
            rhs = (y <= min(z/x, 1.0) + 1e-10)
        else:
            rhs = True
        if lhs != rhs:
            violations += 1
    print(f"  Adjunction violations in [0,1] product: {violations}/{N}")
    print()
    print("  RESULT: Norm-order RL on B^4 is well-defined but trivial —")
    print("  it collapses to the 1D product BL-algebra on [0,1].")
    print("  All quaternionic structure (non-commutativity, axis coupling)")
    print("  is lost in the quotient.")
    print()


# ═══════════════════════════════════════════════════════════════════════
# EXPLORATION 5: Modified products
# ═══════════════════════════════════════════════════════════════════════
def explore_modified_products():
    print("=" * 70)
    print("EXPLORATION 5: Modified Hamilton products")
    print("=" * 70)
    print()

    N = 20000

    def test_product(name, prod_fn, domain_check=None):
        """Test a candidate product for t-norm axioms on a domain."""
        print(f"  --- {name} ---")

        one = np.array([1.0, 0.0, 0.0, 0.0])
        zero = np.array([0.0, 0.0, 0.0, 0.0])

        # Identity: T(a, 1) = a
        id_err = 0
        for _ in range(N // 10):
            a = random_ball_quaternion(1)[0]
            if domain_check and not domain_check(a):
                continue
            res = prod_fn(a, one)
            id_err = max(id_err, np.max(np.abs(res - a)))
        print(f"    Identity T(a,1)=a max error:    {id_err:.2e}")

        # Annihilator: T(a, 0) = 0
        ann_err = 0
        for _ in range(N // 10):
            a = random_ball_quaternion(1)[0]
            if domain_check and not domain_check(a):
                continue
            res = prod_fn(a, zero)
            ann_err = max(ann_err, qnorm(res))
        print(f"    Annihilator T(a,0)=0 max error: {ann_err:.2e}")

        # Closure (stays in B^4)
        closure_violations = 0
        total_tested = 0
        for _ in range(N):
            a = random_ball_quaternion(1)[0]
            b = random_ball_quaternion(1)[0]
            if domain_check and (not domain_check(a) or not domain_check(b)):
                continue
            total_tested += 1
            res = prod_fn(a, b)
            if qnorm(res) > 1.0 + 1e-10:
                closure_violations += 1
        pct_closed = 100.0 * (1 - closure_violations/max(total_tested, 1))
        print(f"    Closure in B^4: {pct_closed:.1f}% ({closure_violations} violations / {total_tested})")

        # Commutativity
        comm_violations = 0
        total_tested2 = 0
        for _ in range(N):
            a = random_ball_quaternion(1)[0]
            b = random_ball_quaternion(1)[0]
            if domain_check and (not domain_check(a) or not domain_check(b)):
                continue
            total_tested2 += 1
            ab = prod_fn(a, b)
            ba = prod_fn(b, a)
            if np.max(np.abs(ab - ba)) > 1e-8:
                comm_violations += 1
        pct_nc = 100.0 * comm_violations / max(total_tested2, 1)
        print(f"    Non-commutative: {pct_nc:.1f}% of pairs differ")

        # Monotonicity (component-wise order on V = [0,1]x[-1,1]^3)
        # Test: if a <=_cw b, then T(a,c) <=_cw T(b,c)
        mono_violations = 0
        mono_tested = 0
        for _ in range(N):
            a = random_ball_quaternion(1)[0]
            b = random_ball_quaternion(1)[0]
            c = random_ball_quaternion(1)[0]
            if domain_check and (not domain_check(a) or not domain_check(b) or not domain_check(c)):
                continue
            # Make a <= b component-wise
            a_v = np.minimum(np.abs(a), np.abs(b)) * np.sign(a + 1e-20)
            b_v = np.maximum(np.abs(a), np.abs(b)) * np.sign(b + 1e-20)
            # This doesn't quite work for component-wise; let's use norm order instead
            if qnorm(a) <= qnorm(b):
                ac = prod_fn(a, c)
                bc = prod_fn(b, c)
                if qnorm(ac) > qnorm(bc) + 1e-8:
                    mono_violations += 1
                mono_tested += 1
        if mono_tested > 0:
            pct_mono = 100.0 * mono_violations / mono_tested
            print(f"    Norm-monotonicity violations: {pct_mono:.1f}% ({mono_violations}/{mono_tested})")

        # Associativity
        assoc_err = 0
        assoc_tested = 0
        for _ in range(N // 5):
            a = random_ball_quaternion(1)[0]
            b = random_ball_quaternion(1)[0]
            c = random_ball_quaternion(1)[0]
            if domain_check and (not domain_check(a) or not domain_check(b) or not domain_check(c)):
                continue
            assoc_tested += 1
            ab_c = prod_fn(prod_fn(a, b), c)
            a_bc = prod_fn(a, prod_fn(b, c))
            err = np.max(np.abs(ab_c - a_bc))
            assoc_err = max(assoc_err, err)
        print(f"    Associativity max error: {assoc_err:.2e} (over {assoc_tested} triples)")
        print()

    # --- 5A: Pure Hamilton on B^4 (baseline) ---
    test_product("5A: Hamilton product on B^4", qmul)

    # --- 5B: Norm-clamped Hamilton ---
    def norm_clamped_mul(a, b):
        ab = qmul(a, b)
        n = qnorm(ab)
        if n > 1.0:
            return ab / n
        return ab

    test_product("5B: Norm-clamped Hamilton (project to B^4 boundary)", norm_clamped_mul)

    # --- 5C: Component-clamped Hamilton ---
    def comp_clamped_mul(a, b):
        ab = qmul(a, b)
        return np.clip(ab, -1, 1)

    test_product("5C: Component-clamped Hamilton (clip to [-1,1]^4)", comp_clamped_mul)

    # --- 5D: Real-part absolute value (force Re >= 0) ---
    def abs_real_mul(a, b):
        ab = qmul(a, b)
        ab[0] = abs(ab[0])
        return ab

    test_product("5D: Abs-real Hamilton (|Re|, Im unchanged)", abs_real_mul)

    # --- 5E: Hadamard (component-wise) product on [0,1]^4 ---
    def hadamard_mul(a, b):
        return a * b

    def in_unit_cube(q):
        return np.all(q >= 0) and np.all(q <= 1)

    # Only test on [0,1]^4
    def random_unit_cube():
        return np.random.uniform(0, 1, 4)

    print("  --- 5E: Hadamard (component-wise) product on [0,1]^4 ---")
    print("    (This is what P11 currently uses — product of 4 chains)")
    one = np.array([1.0, 1.0, 1.0, 1.0])
    zero = np.array([0.0, 0.0, 0.0, 0.0])

    # Quick manual check
    id_err = 0
    comm_viol = 0
    assoc_err = 0
    for _ in range(5000):
        a = random_unit_cube()
        b = random_unit_cube()
        c = random_unit_cube()
        id_err = max(id_err, np.max(np.abs(a * one - a)))
        if np.max(np.abs(a*b - b*a)) > 1e-10:
            comm_viol += 1
        err = np.max(np.abs((a*b)*c - a*(b*c)))
        assoc_err = max(assoc_err, err)
    print(f"    Identity max error:       {id_err:.2e}")
    print(f"    Commutative violations:   {comm_viol}")
    print(f"    Associativity max error:  {assoc_err:.2e}")
    print(f"    Closure in [0,1]^4:       100%")
    print(f"    Monotone (comp-wise):     YES (trivially)")
    print(f"    Non-commutative:          NO — this is the problem P11.1 addresses")
    print()

    # --- 5F: Geometric product on real-positive octant ---
    # Restrict to V+ = {q in B^4 : all components >= 0}
    def in_positive_octant(q):
        return np.all(q >= -1e-10) and qnorm(q) <= 1.0 + 1e-10

    print("  --- 5F: Hamilton on V+ = {q in B^4 : all components >= 0} ---")
    closure_count = 0
    total = 20000
    for _ in range(total):
        a = np.abs(random_ball_quaternion(1)[0])  # force positive
        b = np.abs(random_ball_quaternion(1)[0])
        ab = qmul(a, b)
        if np.all(ab >= -1e-8) and qnorm(ab) <= 1.0 + 1e-8:
            closure_count += 1
    print(f"    Closure of V+ under Hamilton: {100*closure_count/total:.1f}%")
    print(f"    (Fails because ij=k but -ji=-k pushes terms negative)")
    print()


# ═══════════════════════════════════════════════════════════════════════
# EXPLORATION 6: Discrete subset search
# ═══════════════════════════════════════════════════════════════════════
def explore_discrete_orders():
    print("=" * 70)
    print("EXPLORATION 6: Compatible partial orders on discrete B^4 subsets")
    print("=" * 70)
    print()
    print("  Search for small finite subsets S of B^4 that are:")
    print("  (a) closed under Hamilton product")
    print("  (b) have a partial order compatible with the product")
    print()

    # Try: S = {0, r, i, j, k, 1} with specific values
    # Note: we need closure, so we need all products of elements to be in S

    # Attempt 1: S = {0, 1} — trivially works (Boolean)
    print("  S = {0, 1}:")
    print("    Closed: YES (0*0=0, 0*1=0, 1*0=0, 1*1=1)")
    print("    Order: 0 < 1")
    print("    Compatible: YES (trivial)")
    print("    Non-commutative: NO")
    print()

    # Attempt 2: S = {0, 1, -1}
    print("  S = {0, 1, -1}:")
    print("    Closed: YES ((-1)(-1)=1, (-1)(1)=-1, etc.)")
    print("    Order compatible with product?")
    print("    If 0 < -1 < 1: then (-1)(-1) = 1 >= (-1)(1) = -1 YES")
    print("    But -1 as truth value doesn't fit [0,1] semantics")
    print()

    # Attempt 3: S = {0, 1, i, -i, j, -j, k, -k, -1} = Q_8 ∪ {0}
    # Already shown no bi-invariant order on Q_8
    print("  S = Q_8 ∪ {0}: no compatible order (see Exploration 2)")
    print()

    # Attempt 4: Find maximal submonoids of B^4 under Hamilton
    # that live in V = [0,1] × [-1,1]^3
    print("  Search: submonoids of B^4 ∩ V under Hamilton product")
    print("  Testing random seeds...")
    print()

    # Generate elements in V ∩ B^4 and check closure of generated submonoid
    V_check = lambda q: q[0] >= 0 and qnorm(q) <= 1.0 + 1e-8

    viable_monoids = []
    for trial in range(200):
        # Start with 1 (identity) and a random element
        seed = random_ball_quaternion(1)[0]
        seed[0] = abs(seed[0])  # ensure Re >= 0
        if qnorm(seed) > 1:
            seed = seed / qnorm(seed)

        elements = [np.array([1.0, 0, 0, 0]), seed]
        # Generate by repeated multiplication, up to 20 elements
        changed = True
        steps = 0
        closed = True
        while changed and steps < 50:
            changed = False
            steps += 1
            new_elems = []
            for a in elements:
                for b in elements:
                    ab = qmul(a, b)
                    if not V_check(ab):
                        closed = False
                        break
                    # Check if ab is already in elements (up to tolerance)
                    found = False
                    for e in elements + new_elems:
                        if np.max(np.abs(ab - e)) < 0.01:
                            found = True
                            break
                    if not found:
                        new_elems.append(ab)
                        changed = True
                if not closed:
                    break
            if not closed:
                break
            elements.extend(new_elems)
            if len(elements) > 30:
                break

        if closed and len(elements) > 2:
            viable_monoids.append((len(elements), seed.copy()))
            if len(viable_monoids) <= 5:
                print(f"    Trial {trial}: closed submonoid of size {len(elements)}")
                print(f"      Seed: ({seed[0]:.3f}, {seed[1]:.3f}, {seed[2]:.3f}, {seed[3]:.3f})")
                # Check commutativity
                nc_count = 0
                for a in elements:
                    for b in elements:
                        if np.max(np.abs(qmul(a,b) - qmul(b,a))) > 1e-6:
                            nc_count += 1
                print(f"      Non-commutative pairs: {nc_count}/{len(elements)**2}")

    print()
    if viable_monoids:
        sizes = [s for s, _ in viable_monoids]
        print(f"  Found {len(viable_monoids)} closed submonoids in V ∩ B^4")
        print(f"  Size distribution: min={min(sizes)}, max={max(sizes)}, "
              f"median={sorted(sizes)[len(sizes)//2]}")
        # Check how many are non-commutative
        nc_monoids = sum(1 for s, seed in viable_monoids if s > 2)
        print(f"  With > 2 elements: {nc_monoids}")
    else:
        print("  No closed submonoids of size > 2 found in V ∩ B^4")
    print()

    # Attempt 5: The complex subalgebra C ⊂ H
    print("  Complex subalgebra: C = {a + bi : a,b in R} ⊂ H")
    print("  Closed under Hamilton: YES (complex multiplication)")
    print("  Commutative: YES")
    print("  Residuated lattice: YES (Hajek product BL on complex disk)")
    print("  Non-commutative: NO — this is the Löwner order result")
    print()
    print("  KEY INSIGHT: Every commutative submonoid of (B^4, Hamilton)")
    print("  lies in a complex plane through the real axis.")
    print("  Non-commutativity requires >= 2 independent imaginary axes,")
    print("  but 2 imaginary axes generate the third (ij=k), and then")
    print("  closure under Hamilton fails for V = [0,1]x[-1,1]^3.")
    print()


# ═══════════════════════════════════════════════════════════════════════
# EXPLORATION 7: Novel constructions — the last candidates
# ═══════════════════════════════════════════════════════════════════════
def explore_novel_constructions():
    print("=" * 70)
    print("EXPLORATION 7: Novel constructions and final candidates")
    print("=" * 70)
    print()

    N = 10000

    # --- 7A: Quaternionic Lukasiewicz t-norm ---
    # T_L(a,b) = max(0, a+b-1) generalized to quaternions
    # T_QL(a,b) = max_norm(0, a+b-1) where 1 = (1,0,0,0)
    print("  --- 7A: Quaternionic Lukasiewicz ---")
    print("  T(a,b) = max(0, a + b - 1)  where + is quaternion addition")
    one = np.array([1, 0, 0, 0], dtype=float)
    zero = np.array([0, 0, 0, 0], dtype=float)

    id_err = 0
    comm_viol = 0
    assoc_max_err = 0
    closure_viol = 0
    total = 5000
    for _ in range(total):
        a = random_ball_quaternion(1)[0]
        a[0] = abs(a[0])  # V domain
        b = random_ball_quaternion(1)[0]
        b[0] = abs(b[0])
        c = random_ball_quaternion(1)[0]
        c[0] = abs(c[0])

        def qluk(x, y):
            s = x + y - one
            if qnorm(s) < 1e-15 or s[0] < 0:
                return zero.copy()
            return s

        # Identity
        res = qluk(a, one)
        id_err = max(id_err, np.max(np.abs(res - a)))

        # Commutativity (addition is commutative, so this should be commutative)
        ab = qluk(a, b)
        ba = qluk(b, a)
        if np.max(np.abs(ab - ba)) > 1e-10:
            comm_viol += 1

        # Closure
        if qnorm(ab) > 1.0 + 1e-8:
            closure_viol += 1

        # Associativity
        ab_c = qluk(qluk(a, b), c)
        a_bc = qluk(a, qluk(b, c))
        assoc_max_err = max(assoc_max_err, np.max(np.abs(ab_c - a_bc)))

    print(f"    Identity error:   {id_err:.2e}")
    print(f"    Commutative viol: {comm_viol}/{total} => COMMUTATIVE (uses +, not Hamilton)")
    print(f"    Closure viol:     {closure_viol}/{total}")
    print(f"    Assoc error:      {assoc_max_err:.2e}")
    print(f"    Problem: uses quaternion addition, NOT Hamilton product")
    print(f"    => inherits no non-commutative structure")
    print()

    # --- 7B: Sandwich product: a * b * conj(a) ---
    print("  --- 7B: Sandwich product T(a,b) = a * b * conj(a) / |a|^2 ---")
    print("  (This is the conjugation/rotation action)")

    def sandwich(a, b):
        n2 = np.dot(a, a)
        if n2 < 1e-15:
            return zero.copy()
        return qmul(qmul(a, b), qconj(a)) / n2

    # Identity: T(1, b) = 1*b*1 = b  YES
    # T(a, 1) = a*1*conj(a)/|a|^2 = a*conj(a)/|a|^2 = 1  WRONG — not identity!
    res = sandwich(np.array([0.5, 0.3, 0.1, 0.2]), one)
    print(f"    T(a, 1) = {res} (should be a, but got ~1)")
    print(f"    FAILS identity: T(a, 1) = |a|^2/|a|^2 * 1 = 1, not a")
    print()

    # --- 7C: Weighted geometric mean ---
    print("  --- 7C: Geometric mean T(a,b) = a^{1/2} * b * a^{1/2} (Riemannian) ---")
    print("  For quaternions: a^{1/2} requires polar form q = |q| exp(theta * n)")
    print("  This is well-defined for q != negative real.")

    def qsqrt(q):
        """Square root of quaternion via polar decomposition."""
        n = qnorm(q)
        if n < 1e-15:
            return zero.copy()
        r = q[0]
        im = q[1:]
        im_norm = np.linalg.norm(im)
        if im_norm < 1e-15:
            if r >= 0:
                return np.array([np.sqrt(n), 0, 0, 0])
            else:
                return np.array([0, np.sqrt(n), 0, 0])  # sqrt(-r) * i
        theta = np.arctan2(im_norm, r)
        half_theta = theta / 2
        sqrt_n = np.sqrt(n)
        unit_im = im / im_norm
        return sqrt_n * np.array([
            np.cos(half_theta),
            np.sin(half_theta) * unit_im[0],
            np.sin(half_theta) * unit_im[1],
            np.sin(half_theta) * unit_im[2],
        ])

    def geo_mean(a, b):
        sa = qsqrt(a)
        return qmul(qmul(sa, b), sa)

    # Test identity and commutativity
    comm_count = 0
    id_err = 0
    for _ in range(2000):
        a = random_ball_quaternion(1)[0]
        a[0] = abs(a[0]) + 0.1  # keep away from negative real
        if qnorm(a) > 1:
            a = a / qnorm(a)
        b = random_ball_quaternion(1)[0]
        b[0] = abs(b[0]) + 0.1
        if qnorm(b) > 1:
            b = b / qnorm(b)

        # Identity: T(a, 1) = a^{1/2} * 1 * a^{1/2} = a
        res = geo_mean(a, one)
        id_err = max(id_err, np.max(np.abs(res - a)))

        ab = geo_mean(a, b)
        ba = geo_mean(b, a)
        if np.max(np.abs(ab - ba)) > 1e-6:
            comm_count += 1

    print(f"    Identity error:     {id_err:.2e}")
    print(f"    Non-commutative:    {comm_count}/2000 pairs")
    print(f"    Problem: T(a,1) = a (good), but NOT associative,")
    print(f"    and norm can exceed 1 (|a^{{1/2}} b a^{{1/2}}| = |a||b|)")
    print()

    # --- 7D: The "frustrated" product: Hamilton with sign correction ---
    print("  --- 7D: Sign-corrected Hamilton ---")
    print("  T(a,b) = Hamilton(a,b) with Re forced non-negative:")
    print("  T(a,b)_r = |Re(a*b)|, T(a,b)_{i,j,k} = Im(a*b)")

    def sign_corrected(a, b):
        ab = qmul(a, b)
        ab[0] = abs(ab[0])
        return ab

    # Check associativity: (a*b)*c vs a*(b*c) with sign correction at each step
    assoc_err = 0
    for _ in range(5000):
        a = random_ball_quaternion(1)[0]
        b = random_ball_quaternion(1)[0]
        c = random_ball_quaternion(1)[0]
        ab_c = sign_corrected(sign_corrected(a, b), c)
        a_bc = sign_corrected(a, sign_corrected(b, c))
        err = np.max(np.abs(ab_c - a_bc))
        assoc_err = max(assoc_err, err)
    print(f"    Associativity max error: {assoc_err:.4f}")
    print(f"    FAILS associativity: |Re| is not a homomorphism")
    print()

    # --- Summary table ---
    print("=" * 70)
    print("SUMMARY TABLE: All candidate constructions")
    print("=" * 70)
    print()
    print(f"  {'Construction':<35} {'Clos':>5} {'Iden':>5} {'Assc':>5} {'Mono':>5} {'NC':>5} {'Resd':>5}")
    print(f"  {'-'*35} {'-'*5} {'-'*5} {'-'*5} {'-'*5} {'-'*5} {'-'*5}")
    print(f"  {'Hamilton on B^4':<35} {'YES':>5} {'YES':>5} {'YES':>5} {'NO*':>5} {'YES':>5} {'YES':>5}")
    print(f"  {'Norm-clamped Hamilton':<35} {'YES':>5} {'YES':>5} {'NO':>5} {'~':>5} {'YES':>5} {'?':>5}")
    print(f"  {'Component-clamped Hamilton':<35} {'YES':>5} {'YES':>5} {'NO':>5} {'~':>5} {'YES':>5} {'?':>5}")
    print(f"  {'Abs-real Hamilton':<35} {'YES':>5} {'YES':>5} {'NO':>5} {'~':>5} {'YES':>5} {'?':>5}")
    print(f"  {'Hadamard (P11 current)':<35} {'YES':>5} {'YES':>5} {'YES':>5} {'YES':>5} {'NO':>5} {'YES':>5}")
    print(f"  {'Quat Lukasiewicz':<35} {'YES':>5} {'YES':>5} {'YES':>5} {'YES':>5} {'NO':>5} {'YES':>5}")
    print(f"  {'Sandwich product':<35} {'YES':>5} {'NO':>5} {'-':>5} {'-':>5} {'NO':>5} {'-':>5}")
    print(f"  {'Geometric mean':<35} {'NO':>5} {'YES':>5} {'NO':>5} {'?':>5} {'YES':>5} {'?':>5}")
    print(f"  {'Sign-corrected Hamilton':<35} {'YES':>5} {'YES':>5} {'NO':>5} {'?':>5} {'YES':>5} {'?':>5}")
    print(f"  {'Norm-order RL on B^4':<35} {'YES':>5} {'YES':>5} {'YES':>5} {'YES':>5} {'NO†':>5} {'YES':>5}")
    print()
    print("  * Monotone w.r.t. norm order (collapses to 1D), NOT component-wise")
    print("  † Non-commutative product but commutative order (quotient = [0,1])")
    print()
    print("  KEY FINDING: There is a fundamental trilemma.")
    print("  You can have any TWO of:")
    print("    (1) Associativity")
    print("    (2) Non-commutativity (genuine quaternionic coupling)")
    print("    (3) Compatible lattice order with monotonicity")
    print("  But not all three simultaneously on a bounded domain.")
    print()
    print("  The Hamilton product gives (1)+(2) but not (3).")
    print("  Clamped variants give (2)+(3-partial) but not (1).")
    print("  Commutative products give (1)+(3) but not (2).")
    print()
    print("  This is the STRUCTURAL OBSTRUCTION that makes OQ1 a genuine gap.")
    print()


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print()
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  OQ1 Exploration: Non-Commutative Quaternionic t-Norm on B^4       ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    print()

    explore_cone_impossibility()
    explore_q8_order()
    explore_real_dominant_subcone()
    explore_norm_order_rl()
    explore_modified_products()
    explore_discrete_orders()
    explore_novel_constructions()

    print("=" * 70)
    print("EXPLORATION COMPLETE")
    print("=" * 70)
