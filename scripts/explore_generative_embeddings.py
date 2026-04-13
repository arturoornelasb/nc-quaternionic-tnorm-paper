#!/usr/bin/env python3
"""
Exploration: Can the 5 generative algebras be formally embedded
into (r,i,j,k) with provable k=0 degeneration?

For each algebra, we ask:
  1. Is there a CANONICAL embedding phi: X -> (r,i,j,k)?
  2. Is the embedding MATHEMATICALLY FORCED or a SEMANTIC CHOICE?
  3. Does k=0 cause PROVABLE degeneration (information loss)?
  4. What is the THEOREM, if any?

Author: J. Arturo Ornelas Brand
Date: 2026-04-12
"""

import numpy as np
np.random.seed(42)

print("=" * 70)
print("GENERATIVE ALGEBRA EMBEDDING EXPLORATION")
print("=" * 70)
print()

# =====================================================================
# ALGEBRA 5: Categorical Duality (ALREADY DONE IN P11)
# =====================================================================
print("ALGEBRA 5: Categorical Duality  [DONE in P11]")
print("-" * 50)
print()
print("  Status: Theorem + proof in P11 (thm:catdual)")
print("  Embedding: (-)^op is k because:")
print("    (a) meta-operation (acts on morphisms, not objects)")
print("    (b) generative (produces new theorems)")
print("    (c) self-inverting ((C^op)^op = C)")
print("  k=0 degeneration: without (-)^op, no free dual theorems")
print("  TYPE: structural correspondence, not algebraic embedding")
print("  STRENGTH: moderate — matches K3, K4 (involution form)")
print()

# =====================================================================
# ALGEBRA 6: Kernel/Image
# =====================================================================
print("ALGEBRA 6: Kernel/Image")
print("-" * 50)
print()

# The rank-nullity theorem: for f: V -> W,
#   dim(ker f) + dim(im f) = dim(V)
# This is a DUALITY: knowing one determines the other.
# The kernel measures what f DESTROYS; the image measures what f PRESERVES.

# Candidate embedding for a linear map f: R^n -> R^n
# phi(f) = (r, i, j, k) where:
#   r = rank(f)/n          (fraction preserved = "crystallized truth")
#   i = ???                (potentiality — what could change?)
#   j = ???                (verifiability)
#   k = nullity(f)/n       (fraction destroyed = "what was selected away")
#
# Problem: what goes in i and j? The assignment is ARBITRARY.
# There's no canonical way to extract "potentiality" from a linear map.

# ALTERNATIVE: forget the semantic labels. Ask a STRUCTURAL question:
# Given f: V -> V, what information do you LOSE if you set k=0?
#
# Claim: if k = nullity(f)/n = 0, then f is injective.
# For injective f: the kernel/image algebra degenerates because
# ker(f) = {0} — there's nothing to "audit."
#
# Is this PROVABLE? Yes, trivially: nullity=0 iff ker={0} iff f injective.
# The "degeneration" is: the kernel functor becomes trivial (constant {0}).
#
# But is this an EMBEDDING or just a RELABELING?

# Test: for random matrices, compute the embedding
n = 10
N = 5000
ranks = []
for _ in range(N):
    A = np.random.randn(n, n)
    # Make some rank-deficient
    if np.random.random() < 0.3:
        A[:, -1] = A[:, 0]  # force rank < n
    if np.random.random() < 0.1:
        A[:, -2] = A[:, 1]
    r = np.linalg.matrix_rank(A)
    ranks.append(r)

ranks = np.array(ranks)
print("  Random 10x10 matrices:")
print(f"  Full rank (k=0): {(ranks == n).sum()}/{N} = {100*(ranks==n).mean():.1f}%")
print(f"  Rank deficient (k>0): {(ranks < n).sum()}/{N} = {100*(ranks<n).mean():.1f}%")
print()

# THE REAL QUESTION: is there a theorem here that goes beyond trivial relabeling?
#
# Candidate Theorem: "The kernel functor Ker: Lin(V,V) -> Sub(V)
# carries information independent of the image functor Im: Lin(V,V) -> Sub(V).
# Specifically, for any proper subspace W < V, there exist maps f,g with
# Im(f) = Im(g) = W but Ker(f) != Ker(g). The k-axis encodes this
# independent information."
#
# Proof: Take W with dim(W) = d < n. Pick any complement W' of W in V.
# Define f: V -> V by f|_W' = 0, f|_W'' = id for some other decomposition.
# Multiple choices of W' give different kernels with same image.

# Verify computationally
n = 5
# Two maps with same image but different kernels
e1, e2, e3, e4, e5 = np.eye(5)
# f: kills e4,e5, maps e1,e2,e3 to themselves
f = np.zeros((5,5))
f[:3, :3] = np.eye(3)  # im(f) = span(e1,e2,e3), ker(f) = span(e4,e5)

# g: kills e3,e5, maps e1,e2,e4 to e1,e2,e3
g = np.zeros((5,5))
g[0,0] = 1; g[1,1] = 1; g[2,3] = 1  # im(g) = span(e1,e2,e3), ker(g) = span(e3,e5)

rank_f = np.linalg.matrix_rank(f)
rank_g = np.linalg.matrix_rank(g)
print(f"  f: rank={rank_f}, nullity={5-rank_f}, im=span(e1,e2,e3), ker=span(e4,e5)")
print(f"  g: rank={rank_g}, nullity={5-rank_g}, im=span(e1,e2,e3), ker=span(e3,e5)")
print(f"  Same image: {rank_f == rank_g}")
print(f"  Same kernel: NO (e4 in ker(f) but not ker(g))")
print()

# ASSESSMENT
print("  ASSESSMENT:")
print("  [+] Rank-nullity IS a duality (proven, classical)")
print("  [+] Kernel carries info independent of image (provable)")
print("  [+] k=0 <=> injective <=> kernel trivial (provable)")
print("  [-] The assignment k=nullity is a SEMANTIC CHOICE, not forced")
print("  [-] Why k and not j? No algebraic reason; it's interpretive")
print("  [-] No connection to K4 (self-inversion): Ker(Ker(f)) doesn't typecheck")
print()
print("  VERDICT: PARTIAL. Can prove a degeneration theorem, but the")
print("  embedding is not canonical — it's a labeled reading of rank-nullity.")
print("  The theorem would be about rank-nullity itself, not about k.")
print()

# =====================================================================
# ALGEBRA 7: Eigenvalue Spectrum
# =====================================================================
print("ALGEBRA 7: Eigenvalue Spectrum")
print("-" * 50)
print()

# For A: V -> V, the eigenvalues {lambda_i} classify dynamics:
#   lambda > 0: expansion
#   lambda < 0: reflection
#   lambda complex: oscillation
#   lambda = 0: annihilation
#
# The spectral decomposition A = sum lambda_i P_i (projectors)
# is a GLOBAL property — you can't read it from A's action on one vector.
#
# Candidate embedding phi(A) = (r, i, j, k):
#   r = spectral radius rho(A) / ||A||  (normalized "strength")
#   i = dim(generalized eigenspace for complex eigenvalues) / n
#   j = number of distinct real eigenvalues / n
#   k = number of DISTINCT eigenvalues / n (spectral complexity)
#
# k=0 degeneration: A = lambda * I (scalar). Trivial spectrum.

# Test: spectral complexity for random matrices
n = 10
N = 2000
spectral_complexities = []
for _ in range(N):
    A = np.random.randn(n, n)
    eigenvalues = np.linalg.eigvals(A)
    # "Distinct" eigenvalues (up to tolerance)
    distinct = []
    for ev in eigenvalues:
        if all(abs(ev - d) > 0.1 for d in distinct):
            distinct.append(ev)
    spectral_complexities.append(len(distinct) / n)

sc = np.array(spectral_complexities)
print(f"  Random 10x10 matrices: mean spectral complexity = {sc.mean():.3f}")
print(f"  Min: {sc.min():.3f}, Max: {sc.max():.3f}")
print(f"  k=0 (all eigenvalues equal): {(sc < 0.15).sum()}/{N}")
print()

# Key mathematical fact: The spectral theorem for normal matrices
# gives a canonical decomposition. But for general matrices, the
# Jordan form is the canonical form, and it's not a "spectrum" in
# the same sense.
#
# The "requires k" claim: classifying the MODES of A requires
# observing A's global structure. A single vector Av tells you
# about one direction; the full spectrum requires all of them.
#
# Is this formalizeable? YES, via:
# Theorem: No linear functional L: Mat(n) -> R determines the spectrum.
# Proof: The spectrum is not a linear function of the matrix entries.
# (trivial: spec(A+B) != spec(A) + spec(B))
#
# But that's a theorem about spectral theory, not about k.

# More relevant: the spectral decomposition SEPARATES modes that are
# COUPLED in the original operator. This is a "selection" operation:
# choosing which modes to look at.

# Check: is the spectral decomposition self-inverting?
# A = P D P^{-1}, so "decompose" sends A to (D, P).
# "Recompose" sends (D, P) to A.
# decompose(decompose(A)) doesn't typecheck (D is already diagonal).
# So: NO self-inversion. K4 doesn't match.

print("  ASSESSMENT:")
print("  [+] Spectrum is a GLOBAL property (not readable from one vector)")
print("  [+] k=0 <=> scalar matrix <=> trivial spectrum (provable)")
print("  [+] Spectral decomposition = 'selecting modes' (interpretive)")
print("  [-] Embedding is ARBITRARY (many reasonable choices for r,i,j,k)")
print("  [-] No self-inversion: decompose(decompose(A)) doesn't typecheck")
print("  [-] The claim 'classifying modes requires k' is SEMANTIC, not algebraic")
print("  [-] WHY k and not j? No formal reason.")
print()
print("  VERDICT: WEAK. The spectral decomposition is genuinely meta-level")
print("  (global, not local), but the connection to k is interpretive.")
print("  No canonical embedding. No K4 match.")
print()

# =====================================================================
# ALGEBRA 8: Brouwer Fixed Point
# =====================================================================
print("ALGEBRA 8: Brouwer Fixed Point")
print("-" * 50)
print()

# Brouwer: every continuous f: D^n -> D^n has a fixed point x*: f(x*)=x*.
# The existence is guaranteed but FINDING x* requires search/iteration.
#
# Candidate embedding: ???
# The fixed point theorem is about EXISTENCE, not about a specific map.
# What would phi map? The function f? The fixed point x*? The iteration?
#
# If phi(f) = (r, i, j, k):
#   r = |f(x*) - x*| = 0 (at fixed point, trivially)
#   k = ... "how hard it is to find x*"? This is COMPUTATIONAL, not algebraic.
#
# The Brouwer theorem itself is NON-CONSTRUCTIVE (in general).
# You can't compute x* from f in general — you need to SEARCH.
# This is the connection to k: search is selection.
#
# But this is an argument about COMPUTABILITY, not about algebra.

# Test: fixed point iteration for random contractions on [0,1]
N = 1000
iterations_needed = []
for _ in range(N):
    # Random Lipschitz contraction on [0,1]
    c = np.random.uniform(0.1, 0.9)  # contraction constant
    b = np.random.uniform(0, 1-c)     # offset
    f = lambda x, c=c, b=b: c*x + b
    # Fixed point: x* = b/(1-c)
    x_star = b / (1 - c)

    # Iterate from random start
    x = np.random.uniform(0, 1)
    for itr in range(1000):
        x = f(x)
        if abs(x - x_star) < 1e-10:
            iterations_needed.append(itr + 1)
            break

iters = np.array(iterations_needed)
print(f"  Contractions on [0,1]: mean iterations to fixed point = {iters.mean():.1f}")
print(f"  Range: [{iters.min()}, {iters.max()}]")
print()

# Is there a THEOREM?
# "Finding the fixed point requires iteration" — yes, but that's
# about computation, not about the quaternionic framework.
#
# The claim "Brouwer fixed point requires k" really means:
# "The EXISTENCE of a fixed point is provable (r,i,j suffice),
#  but LOCATING it requires selection (k > 0)."
#
# This is a claim about CONSTRUCTIVE vs NON-CONSTRUCTIVE proof,
# not about embeddings.

print("  ASSESSMENT:")
print("  [+] Fixed point existence is non-constructive (classical result)")
print("  [+] Finding x* requires iteration = 'search' = 'selection'")
print("  [-] NO embedding: what does phi map? f? x*? The iteration?")
print("  [-] The connection to k is about COMPUTATION, not algebra")
print("  [-] No self-inversion: iterating the iteration is just more iteration")
print("  [-] 'Search requires selection' is TAUTOLOGICAL, not a theorem")
print()
print("  VERDICT: NOT FORMALIZABLE as an embedding.")
print("  The connection to k is philosophical (selection = search).")
print("  No mathematical content beyond Brouwer's theorem itself.")
print()

# =====================================================================
# ALGEBRA 9: Homology / Betti Numbers
# =====================================================================
print("ALGEBRA 9: Homology / Betti Numbers")
print("-" * 50)
print()

# Homology: H_n(X) = Ker(d_n) / Im(d_{n+1})  (chain complex)
# Betti numbers: beta_n = rank(H_n(X))
# These measure "holes" of dimension n.
#
# KEY OBSERVATION: Homology IS ALREADY a kernel/image construction!
# H_n = Ker/Im. So it's a special case of Algebra 6.
#
# But with extra structure: it's a FUNCTOR from Top to Ab (or GrAb).
# And it measures ABSENCE (holes = what's missing).
#
# Candidate embedding phi(X) = (r, i, j, k):
#   r = 1 if X is connected, 0 otherwise (beta_0)
#   k = sum of beta_n for n >= 1 (total topological complexity)
#
# k=0 degeneration: all beta_n = 0 for n >= 1.
# This means X is acyclic (no holes). If simply connected, X is
# contractible by Whitehead's theorem (under mild conditions).
#
# Is this a THEOREM? Only in the sense that "contractible spaces
# have trivial homology" (a well-known fact).

# Computational illustration: simplicial complexes
# Triangle: H_0 = Z, H_1 = 0, H_2 = 0. beta = (1, 0, 0). k=0.
# Circle: H_0 = Z, H_1 = Z. beta = (1, 1). k=1.
# Torus: H_0 = Z, H_1 = Z^2, H_2 = Z. beta = (1, 2, 1). k=3.
# Klein bottle: H_0 = Z, H_1 = Z + Z/2. beta_1 depends on coefficients.

spaces = [
    ("Point",       [1, 0, 0],     0),
    ("Interval",    [1, 0, 0],     0),
    ("Circle S^1",  [1, 1, 0],     1),
    ("Sphere S^2",  [1, 0, 1],     1),
    ("Torus T^2",   [1, 2, 1],     3),
    ("Klein btl",   [1, 1, 0],     1),  # over Q
    ("RP^2",        [1, 0, 0],     0),  # over Q
    ("Wedge S1vS1", [1, 2, 0],     2),
]

print("  Space          beta_0  beta_1  beta_2  k=sum(beta_n>0)")
print("  " + "-" * 55)
for name, betti, k in spaces:
    b = betti + [0] * (3 - len(betti))
    print(f"  {name:<15} {b[0]:>5}   {b[1]:>5}   {b[2]:>5}   {k:>5}")
print()

# The mathematical content:
# 1. Homology = Ker(d)/Im(d) — a kernel/image construction (Algebra 6)
# 2. Betti numbers measure topological complexity
# 3. beta = 0 iff contractible (for simply connected CW complexes)
#
# But: the assignment k = sum(beta) is ARBITRARY.
# Why not k = beta_1? Or k = chi (Euler characteristic)?
# No algebraic reason to prefer one over another.
#
# MORE FUNDAMENTALLY: homology is a FUNCTOR, and functors are
# compositional. H_n(f o g) relates to H_n(f) and H_n(g).
# This is algebraic structure, not "selection."
#
# The "requires meta-knowledge" claim: to compute H_n(X), you need
# the GLOBAL structure of X (all simplices and their boundary maps).
# This IS true — homology is not a local invariant.
# But "global information" != "selection."

print("  ASSESSMENT:")
print("  [+] Homology measures ABSENCE (holes) — genuinely meta-level")
print("  [+] Homology IS kernel/image (Algebra 6 special case)")
print("  [+] Betti=0 <=> contractible (under conditions) — provable")
print("  [+] Requires global structure (not locally computable)")
print("  [-] Embedding is ARBITRARY (which Betti number = k?)")
print("  [-] 'Meta-knowledge' != 'selection' — interpretive leap")
print("  [-] No self-inversion: H_n(H_n(X)) doesn't typecheck")
print("  [-] As functor, homology is COMPOSITIONAL, not 'selective'")
print()
print("  VERDICT: PARTIAL, but inherits from Algebra 6 (kernel/image).")
print("  The extra claim 'measures absence' is interpretive.")
print("  Provable content = Algebra 6 + functoriality (both classical).")
print()

# =====================================================================
# SYNTHESIS: What can actually be proved?
# =====================================================================
print("=" * 70)
print("SYNTHESIS")
print("=" * 70)
print()
print("For each algebra, the question is: is there a THEOREM, or just")
print("a SEMANTIC READING?")
print()
print("  Algebra   Canonical     k=0 degen    K4 match   Theorem?")
print("            embedding?    provable?    (self-inv)")
print("  " + "-" * 60)
print("  5 CatDual    ~YES       ~YES         YES(invol) DONE in P11")
print("  6 Ker/Im     NO*        YES(trivial) NO         PARTIAL**")
print("  7 Eigenval   NO         YES(trivial) NO         WEAK")
print("  8 Brouwer    NO         N/A          NO         NO")
print("  9 Homology   NO*        YES(trivial) NO         PARTIAL***")
print()
print("  * There IS a natural quantity (nullity, Betti) but its assignment")
print("    to the k-axis rather than j or i is a semantic choice.")
print()
print("  ** The kernel/image independence theorem IS provable:")
print("     'Ker carries information not contained in Im.'")
print("     But this is a theorem about LINEAR ALGEBRA, not about k.")
print()
print("  *** Homology inherits from Algebra 6 (it IS kernel/image).")
print("     The extra 'measures absence' claim is interpretive.")
print()
print("KEY INSIGHT: The 'generative algebras require k' claim has")
print("TWO layers:")
print()
print("  Layer 1 (mathematical): Each generative algebra involves a")
print("  META-OPERATION — one that acts on global structure rather than")
print("  local content. This is PROVABLE for each algebra individually:")
print("    6: Ker(f) requires knowing ALL of f, not just f(v)")
print("    7: spectrum requires ALL eigenvalues, not just one")
print("    8: fixed point existence is non-constructive")
print("    9: Betti numbers require global chain complex")
print()
print("  Layer 2 (interpretive): The meta-operation IS the k-axis.")
print("  This is a SEMANTIC IDENTIFICATION, not a theorem.")
print("  There is no algebraic reason why 'meta-operation' = k")
print("  rather than 'meta-operation' = j or some other axis.")
print()
print("P11.2 can prove Layer 1 for each algebra. Layer 2 is the")
print("interpretive claim that connects them to the quaternionic framework.")
print()
print("RECOMMENDATION FOR P11.2:")
print("  1. Define 'meta-operation' formally (operates on global structure)")
print("  2. Prove each generative algebra is a meta-operation (Layer 1)")
print("  3. Prove each cancellative algebra is NOT a meta-operation")
print("  4. State Layer 2 as an IDENTIFICATION, not a theorem")
print("  5. The paper's contribution = the formal cancellative/generative")
print("     split, not the quaternionic embedding")
print()
print("ALTERNATIVE (stronger but harder):")
print("  Define k operationally as 'the component that vanishes when")
print("  the system has no access to its own global structure.'")
print("  Then 'requires k' becomes: 'requires global self-access.'")
print("  Each generative algebra provably requires global self-access.")
print("  Each cancellative algebra provably does not.")
print("  This makes the identification less arbitrary — but it changes")
print("  what k MEANS from 'recursive selection' to 'global self-access.'")
