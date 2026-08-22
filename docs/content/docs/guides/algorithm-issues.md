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

The [Algorithm](../algorithm/) page describes what SymCurve computes. This page describes
where that computation meets the limits of floating-point arithmetic, and the places where
the implementation makes a choice that a user should know about. Every figure quoted here
was measured rather than estimated, on random sequence in a release build.

The short version: the curvature values SymCurve produces are accurate to roughly one part
in \(10^{8}\), which is below the precision of the bigWig format they are usually written
to. None of the issues below make the output wrong. They do explain why two runs over the
same sequence can disagree in the last digit.

## The accumulated twist grows without bound

The twist-sum \(T_i\) is defined as a running total over every 3-mer from the start of the
sequence:

\[T_i = \sum_{j=1}^{i}\Omega_j\]

With the standard twist matrix every \(\Omega_j\) is \(0.598647428\), so \(T_i\) grows by
about \(0.6\) radians per base and never comes back down. Over a genome that reaches
values no one intends to represent in a double:

| sequence | \(T_n\) (radians) | full turns | ulp of \(T_n\) |
| --- | --- | --- | --- |
| 1 Mbase | \(5.99\times10^{5}\) | 95 thousand | \(1.2\times10^{-10}\) |
| chr1, 250 Mbase | \(1.50\times10^{8}\) | 23.8 million | \(3.0\times10^{-8}\) |
| genome, 3.1 Gbase | \(1.86\times10^{9}\) | 295 million | \(2.4\times10^{-7}\) |

Only \(\sin(T_i)\) and \(\cos(T_i)\) are ever used, so the turns themselves carry no
information. What they cost is precision: each addition rounds relative to a magnitude
that keeps growing, and the errors compound. Measured against a reference that computes
\(T_i\) with a single rounding instead of \(i\) of them, an unbounded running total reaches
an absolute error of \(1.75\times10^{-2}\) in \(\sin(T_i)\) by 50 Mbase. That is a
percent-level error in a quantity bounded by 1.

SymCurve keeps \(T_i\) reduced into \([0, 2\pi)\) as it goes, which is mathematically
identical and holds the error near \(7\times10^{-9}\) over the same run.

### Why the output is far less affected than that suggests

A percent-level error in \(\sin(T_i)\) sounds fatal, and it is worth being precise about
why it is not. Curvature is a distance between two locally averaged points:

\[\kappa_i = \lambda\sqrt{(\overline{x}_{i+b}-\overline{x}_{i-b})^2 + (\overline{y}_{i+b}-\overline{y}_{i-b})^2}\]

An error in \(T_i\) rotates the step vector \((dx_i, dy_i)\). If that error changes slowly,
neighbouring steps are rotated by nearly the same angle, so the traced path is rotated
almost rigidly over the span of a window — and a rotation does not change the distance
between two points on it. The drift therefore cancels almost entirely at the scale
curvature is measured over.

Measured end to end, comparing a bounded twist-sum against an unbounded one:

| sequence | largest relative change in \(\kappa\) |
| --- | --- |
| 1 Mbase | \(4.5\times10^{-9}\) |
| 20 Mbase | \(1.5\times10^{-8}\) |

So the practical effect is small, but it grows with sequence length, and at chromosome
scale it approaches the \({\sim}6\times10^{-8}\) relative precision of the 32-bit floats
that bigWig stores. Bounding the twist-sum costs nothing and removes the question.

## The rolling sums ratchet

The rolling averages \(\overline{x}_i\) and \(\overline{y}_i\) are computed incrementally:
each step subtracts the coordinate leaving the window and adds the one entering it. This is
what makes the calculation \(O(1)\) per base instead of \(O(a)\), and it is the reason
SymCurve can stream a genome in bounded memory.

The cost is that a running total never sheds the rounding of the values that have already
left it. Subtracting a coordinate does not undo the rounding its addition caused, so the
error only accumulates. Measured over 20 million coordinates with an 11-base window,
against sums recomputed fresh at every position:

- largest absolute drift: \(1.13\times10^{-7}\)
- largest relative drift: \(5.4\times10^{-8}\)

That is, again, right at the precision of a 32-bit float. SymCurve rebuilds the rolling
sums from the window contents every 65,536 items. Rebuilding costs one pass over about a
hundred values, so amortised it is a fraction of a percent of the work, and it bounds the
drift instead of letting it grow with the length of the sequence.

## Results depend slightly on `--max-memory`

SymCurve divides each piece of sequence into chunks that are scored independently. This is
what allows a single chromosome to be spread across threads and to be held in bounded
memory. It is exact in principle: a chunk that starts partway into a piece begins with the
wrong coordinate and the wrong accumulated twist, but the first is a translation of the
traced path and the second a rotation of it, and the distances that curvature is made of
survive both. Each chunk reads \(a+b+1\) bases of context beyond its own range at each end,
so no score is ever computed without full context.

What is not preserved is the exact rounding. A smaller chunk accumulates twist over a
shorter run, so its arithmetic rounds differently. Since `--max-memory` sets the chunk
size, two runs at different budgets can differ in the last digit of some values. Measured
over 2 million bases, comparing `--max-memory 8G` against `--max-memory 1M`:

- positions: identical
- values differing at all in their printed form: 0.15%
- largest absolute difference: \(1\times10^{-6}\)
- largest relative difference: \(2\times10^{-7}\)

If you need byte-identical output across runs, pass the same `--max-memory` value. The
differences are below the precision bigWig stores, but bedGraph is written as text and will
show them.

## The tilt term has a sign discrepancy

The equations on the Algorithm page, and the original Perl implementation, both write the
tilt contribution as

\[
\begin{aligned}
dx_i &= \rho_i \sin(T_i) + \tau_i \sin(T_i - \pi/2) \\
dy_i &= \rho_i \cos(T_i) + \tau_i \cos(T_i - \pi/2)
\end{aligned}
\]

The Rust implementation uses \(T_i + \pi/2\) rather than \(T_i - \pi/2\). Because the
supplied tilt matrix \(\tau\) is uniformly zero, both terms vanish and the two forms
produce identical results today. The discrepancy is dormant rather than harmless: it would
change every value the moment a non-zero tilt matrix were supplied, which is exactly what
the planned `--matrices` option would allow. It is recorded here so that whoever implements
that option resolves the sign first.

## `--curve-step-two` cannot be set independently

The reference implementation takes two rolling-mean parameters, `stepone` and `steptwo`. Its
weighting only works out when they satisfy \(\mathtt{stepone} = \mathtt{steptwo} + 2\),
which the defaults (6 and 4) do. SymCurve derives both ends of the window from a single
parameter taken from `--curve-step-one`, so `--curve-step-two` has no independent effect. A
value that disagrees with `--curve-step-one` produces a warning rather than being silently
ignored.

## Unscoreable bases end a window rather than being skipped

Bases that are not A, C, G or T cannot be looked up in the 3-mer matrices. SymCurve splits
the sequence at them rather than deleting them, because the calculation slides a window
over a running coordinate sum: if unknown bases were simply removed, a window could span
the join and produce a value from bases that are far apart in the real sequence.

Splitting has a cost of its own. Each piece loses \(a+b+1\) bases of scores at each end,
since the windows need context on both sides. A sequence broken into many short pieces
therefore yields fewer scores than its length suggests, and pieces shorter than
\(2(a+b+1)\) bases yield none at all.

Lowercase `acgt` is not affected by this. Soft-masking marks repetitive DNA and is an
annotation rather than missing data, so it is scored exactly as its uppercase form.

## Output precision

Both output formats store values as 32-bit floats, which carry about seven decimal digits,
or a relative precision near \(6\times10^{-8}\). The calculation is carried out in 64-bit
throughout and only narrowed on the way out. Every effect described on this page is at or
below that threshold, which is the sense in which they do not affect the result: they are
smaller than the format can record.
