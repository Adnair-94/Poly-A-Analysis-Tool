# Version 2.1 update notes

## Problem addressed

The earlier library generator placed optional C/U residues after the A-run and
inherited the terminal state assigned by the in-silico T1 digest. In some
sequences, C/U candidate masses were approximately 79 Da above the observed
masses. This pattern is consistent with a one-phosphate terminal-state mismatch,
not with C/U positional order.

The monoisotopic mass difference between otherwise identical oligonucleotides
with one additional terminal phosphate is:

```text
HPO3 = 79.9663309 Da
```

Changing `CUU-A(n)` to `A(n)-CUU` does not change intact mass when composition
and termini are unchanged. The update therefore separates two independent
model dimensions:

1. C/U placement relative to the A-run.
2. 5′ and 3′ terminal chemistry.

## Implementation decisions

- Automatic T1 terminal assignment remains the default.
- The user can override 5′ and/or 3′ termini for tail-positive fragments.
- Both automatic and effective terminal states are exported.
- C/U placement can be 5′, 3′ or both.
- Both placement candidates are intentionally retained as ambiguous when they
  are isobaric.
- A phosphate-offset diagnostic is added to annotated output and mode-mass QC.
- DNA coding and DNA template strands are handled explicitly rather than by a
  silent T→U substitution alone.

## Validation included

The test suite verifies:

- coding DNA → RNA conversion;
- template/antisense DNA reverse-complement conversion;
- +79.9663309 Da for a 5′ or 3′ terminal phosphate;
- cyclic-phosphate mass convention;
- equality of 5′-CUU and 3′-CUU intact masses;
- manual terminal override propagation into the library;
- exact mass matching;
- detection of a −79.9663309 Da phosphate offset;
- explicit isobaric duplication when both C/U positions are generated.
