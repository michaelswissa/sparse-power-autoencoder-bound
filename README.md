# Sparsity-dependent power-autoencoder bound

This package contains a mathematical research draft by Michael Swissa, exact CPU verification code, and a machine-readable result summary. No affiliation or credentials are asserted. The manuscript has not received external human peer review; independent expert checking is recommended before relying on it. Originality remains provisional.

## Contents

- `manuscript.md` and `manuscript.pdf`: theorem, proof, references, scope and unreviewed-draft status.
- `verify.py`: original finite-case verification program, requiring only Python's standard library.
- `results.json`: recorded exact-check results.

## Reproduce

Use Python 3.10 or later:

```sh
python3 verify.py --output results.json
```

No package installation, network, GPU, paid service or model training is needed. The recorded run took about four CPU seconds and 22 MB peak memory. Timings vary by machine.

## Validation scope

The run passed 576 matrix/distribution/probability/power cases and 311,525 counted checks. All tested expectations and projections use exact rational arithmetic. It covers powers 3, 5 and 7, input dimensions at most 4, three symmetric finite-support laws, and probabilities 0, 1/100, 1/3 and 1.

The 74,986 moment-state visits include repeated enumerations; separate linear enumerations are excluded. The final inequalities take the easy nonpositive-delta branch because the theorem's constant is loose in these dimensions. Intermediate geometric, moment, interpolation, centering and remainder inequalities provide the substantive falsification tests. There are 405 strictly positive near-optimal scalar rescalings. The three controls are direct counterexamples to false strengthenings or violated assumptions, not mutations of the complete harness.

Finite verification is not a proof, formal proof verification, external peer review, or a training experiment. See the manuscript for the full mathematical argument and model restrictions. The capacity/rank ingredient is attributed to Scherlis et al.; no new capacity inequality is claimed.
