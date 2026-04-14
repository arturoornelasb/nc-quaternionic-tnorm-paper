# Toward a Non-Commutative Residuated Lattice from Quaternion Multiplication

**Author:** J. Arturo Ornelas Brand — arturoornelas62@gmail.com
**Paper DOI:** [10.5281/zenodo.19561407](https://doi.org/10.5281/zenodo.19561407) (v0.1.0, all versions: [10.5281/zenodo.19561406](https://doi.org/10.5281/zenodo.19561406))
**Repository archive DOI:** [10.5281/zenodo.19563496](https://doi.org/10.5281/zenodo.19563496) (v0.1.0, all versions: [10.5281/zenodo.19563495](https://doi.org/10.5281/zenodo.19563495))
**Status:** Published on Zenodo.

## Result

This paper investigates whether the Hamilton product of quaternions can serve as a genuinely **non-commutative conjunction** for the G-lattice introduced in the companion paper (P11). Two types of results are established.

**Positive results on V = [0,1] × [-1,1]^3:**

- Bilateral identity, bilateral annihilation, and associativity
- Left and right implications via quaternionic division, both satisfying modus ponens
- Under every commutative restriction (Boolean, fuzzy, modal), the two implications coincide and the product reduces to the standard commutative t-norm
- Non-commutativity is **emergent**: it arises only in the full four-dimensional space, not in any proper restriction

**Fundamental obstructions:**

- The Hamilton product is not closed on V and is not monotone with respect to component-wise order
- No total order on H is compatible with Hamilton multiplication
- The component-wise lattice order and the Hamilton monoid cannot coexist on V
- On the full ball B^4 = {q ∈ H : |q| ≤ 1}, four independent impossibility results rule out any non-trivial compatible lattice order

## The Quaternionic Trilemma

The results culminate in a **theorem** on V and on B^4:

> For any monoidal operation on a subset of V = [0,1] × [-1,1]^3, at most two of {associativity, non-commutativity, compatible lattice order} can hold simultaneously.

- **Theorem (trilemma-V):** De Moivre's theorem shows that every non-real element's iterates eventually exit V, so the only sub-monoids of V under Hamilton product are commutative.
- **Theorem (trilemma-B4):** The canonical non-commutative sub-monoid on B^4 (equal-norm generators) admits no compatible lattice order.

## Repository structure

```
README.md                              This file
paper/
  nc_quaternionic_tnorm.tex            Paper source
  nc_quaternionic_tnorm.pdf            Compiled paper (16 pp)
scripts/
  verify_trilemma.py                   Main trilemma verification
  b4_trilemma_solver.py                B^4 trilemma proof (equal-norm generators)
  b4_lattice_solver.py                 Lattice-order impossibility on B^4
  b4_proof_verify.py                   Independent verification of the four impossibility results
  explore_generative_embeddings.py     Exploration of non-commutative embeddings
  explore_oq1.py                       Initial exploration of open question OQ1 (from P11)
  gap_analyze.py                       Structural gap characterization
  gap_deep.py                          Deep gap analysis (higher-dimensional restrictions)
  gap_depth4.py                        Full four-dimensional gap structure
  gap_test.py                          Automated gap tests
```

## Dependencies

- Python 3.8+
- NumPy

## Running

```bash
# Main result: trilemma verification
python scripts/verify_trilemma.py

# B^4 proofs (the four independent impossibility results)
python scripts/b4_trilemma_solver.py
python scripts/b4_lattice_solver.py
python scripts/b4_proof_verify.py

# Gap analysis
python scripts/gap_analyze.py
python scripts/gap_deep.py
```

## Companion paper

This paper extends the G-lattice framework established in:

- **P11** — *Quaternionic Logic: A G-Lattice Unifying Boolean and Fuzzy Frameworks*. DOI: [10.5281/zenodo.19562014](https://doi.org/10.5281/zenodo.19562014). The monoidal operation there is component-wise min; the present paper investigates whether Hamilton multiplication can replace it to obtain a genuinely non-commutative conjunction.

## Citation

Paper:

```
Ornelas Brand, J.A. (2026). Toward a Non-Commutative Residuated Lattice
from Quaternion Multiplication (v0.1.0). Zenodo.
https://doi.org/10.5281/zenodo.19561407
```

Concept DOI (all versions): [10.5281/zenodo.19561406](https://doi.org/10.5281/zenodo.19561406)

Repository archive:

```
Ornelas Brand, J.A. (2026). Toward a Non-Commutative Residuated Lattice
from Quaternion Multiplication (repository) (v0.1.0). Zenodo.
https://doi.org/10.5281/zenodo.19563496
```

Concept DOI (all versions): [10.5281/zenodo.19563495](https://doi.org/10.5281/zenodo.19563495)

## License

Business Source License 1.1 (BSL 1.1). Non-production use is permitted. On the Change Date (**2030-04-13**), or the fourth anniversary of the first publicly available distribution of a given version, whichever comes first, the work becomes available under the Change License (**Apache License, Version 2.0**).

For alternative licensing arrangements, contact arturoornelas62@gmail.com. See [`LICENSE`](LICENSE) for full terms.
