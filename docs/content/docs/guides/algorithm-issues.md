---
title: "Algorithm Issues"
description: "Known numerical limits and behavioural quirks of the SymCurve calculation."
summary: ""
date: 2026-08-22T00:00:00+00:00
lastmod: 2026-08-22T00:00:00+00:00
draft: false
menu:
  docs:
    parent: ""
    identifier: "algorithm-issues-9d2f1c7a4b8e6053f19c27ad5e30b4c8"
weight: 815
toc: true
seo:
  title: "" # custom title (optional)
  description: "" # custom description (recommended)
  canonical: "" # custom canonical URL (optional)
  noindex: false # false (default) or true
---

Curvature values are accurate to about one part in \(10^{8}\), which is finer than the
32-bit floats that bigWig and bedGraph store. Nothing below makes the output wrong, but
each is a place where the implementation makes a choice worth knowing about.

## Numerical precision

**The twist-sum is kept bounded.** The accumulated twist \(T_i\) grows by about 0.6
radians per base and only its sine and cosine are used, so SymCurve reduces it into
\([0, 2\pi)\) as it goes. Left unbounded it loses precision over a chromosome; bounded, the
effect on curvature is below \(10^{-8}\).

**The rolling sums are rebuilt periodically.** Running totals accumulate rounding error
from values that have already left the window. SymCurve rebuilds them from the window
contents every 65,536 positions, which bounds the drift at negligible cost.

**Results depend slightly on `--max-memory`.** Sequence is scored in chunks whose size is
set by the memory budget. Chunking is exact in principle, but different chunk sizes round
differently, so runs at different budgets can differ in the last digit of a few values
(largest relative difference measured: \(2\times10^{-7}\)). For byte-identical output,
use the same `--max-memory`.

## Behavioural quirks

**The tilt term is inert.** The supplied tilt matrix is zero, so curvature is determined
entirely by roll and twist. The term is implemented and tested for anyone supplying a
non-zero matrix.

**The symmetry stage has a wide margin.** Symmetry needs `--symcurve-win` curvature values
on each side of a dyad, and each of those needs its own context, so `--stage symmetry`
yields nothing for the first and last \(a+b+1+\mathtt{win}\) bases of a piece (122 at the
defaults). The reference implementation reserves this full margin even though it only
uses half of it; that is reproduced so positions match.

**Nucleosome calls reproduce three reference quirks.**

- A phantom call at position 0 rejects every dyad at or below 177, so none is ever selected.
- GFF coordinates are zero-based (`dyad - 73` printed directly), one base left of where a
  browser will place them. Add one for spec-conformant GFF.
- The reference breaks tied scores by hash order, which varies between its own runs.
  SymCurve breaks ties by ascending dyad, so tied cases may not match.

The selection itself is the same greedy rule as the reference, implemented in
\(O(n \log n)\) rather than \(O(n^2)\).

**`--curve-step-two` cannot be set independently.** Both window bounds derive from
`--curve-step-one`, since the reference's weighting only works when
\(\mathtt{stepone} = \mathtt{steptwo} + 2\). A disagreeing value produces a warning.

**Non-ACGT bases split the sequence.** Rather than deleting unknown bases, SymCurve splits
the sequence at them, so a window never spans a gap. Each piece loses \(a+b+1\) bases of
scores at each end. Lowercase soft-masked bases are scored normally.
