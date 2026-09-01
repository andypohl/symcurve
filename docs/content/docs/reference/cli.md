---
title: "CLI"
description: "The symcurve command-line interface: arguments, options, output formats and stages."
summary: ""
date: 2026-08-31T00:00:00+00:00
lastmod: 2026-08-31T00:00:00+00:00
draft: false
menu:
  docs:
    parent: ""
    identifier: "cli-7c1e9a3f5b2d4e6081a9c3f7e5d2b4a6"
weight: 905
toc: true
seo:
  title: "" # custom title (optional)
  description: "" # custom description (recommended)
  canonical: "" # custom canonical URL (optional)
  noindex: false # false (default) or true
---

```
symcurve [OPTIONS] <INPUT> <OUTPUT>
```

`INPUT` is a FASTA file. `OUTPUT` is where the result goes; its extension chooses both the
file format and, indirectly, what can be written to it. `symcurve --help` prints every
option with its default.

## Input

Any FASTA file works. Lowercase (soft-masked) bases are scored as their uppercase form.
Bases other than A, C, G and T split the sequence, so no score is computed across them;
see [Algorithm Issues](../../guides/algorithm-issues/) for what that costs at the edges.

If an index exists alongside the input (`input.fa.fai`, as produced by `samtools faidx`),
it is used automatically: records are read a window at a time instead of whole, which
lowers memory for large chromosomes, and bigWig output no longer needs a first pass over
the file to learn the chromosome sizes. A missing index is fine; a malformed one is an
error.

## Output

The format is taken from the output extension. An unrecognised extension is refused
before any work is done.

| extension | format | contents |
| --- | --- | --- |
| `.bw`, `.bigWig` | bigWig | one value per base, indexed and compressed |
| `.bedGraph`, `.bg` | bedGraph | one value per base as text: `chrom start end value`, 0-based half-open |
| `.gff`, `.gff2`, `.gff3` | GFF | one nucleosome call per line, with the called sequence in the attribute column |

Values are written as 32-bit floats in both signal formats. The GFF feature column is
`act_nucleosome` with `--roll simple` and `stat_nucleosome` with `--roll active`, as in
the reference implementation.

## Stages

`--stage` selects how far along the pipeline to go. Each stage is computed from the one
before it.

| `--stage` | writes | formats |
| --- | --- | --- |
| `curvature` (default) | curvature at each base | bigWig, bedGraph |
| `symmetry` | symmetry of curvature around each dyad | bigWig, bedGraph |
| `calls` | every nucleosome call, which may overlap | GFF |
| `final-calls` | non-overlapping calls, chosen greedily by score | GFF |

A signal stage with a `.gff` output, or a call stage with a `.bw` output, is an error.

## Options

**Curvature**

| option | default | meaning |
| --- | --- | --- |
| `--roll simple\|active` | `simple` | roll matrix: `simple` is the reference's DNase state, `active` its nucleosome state |
| `--curve-step` | `15` | distance between the two averaged points |
| `--curve-step-one` | `6` | rolling-average window; `--curve-step-two` is derived from it |
| `--curve-step-two` | `4` | accepted for compatibility; a value other than `--curve-step-one` minus 2 is ignored with a warning |
| `--curve-scale` | `0.33335` | scale factor applied to every curvature value |

**Symmetry** (apply to `--stage symmetry` and the call stages)

| option | default | meaning |
| --- | --- | --- |
| `--symcurve-win` | `101` | curvature values on each side of the dyad |
| `--symcurve-step` | `1` | spacing of the mirrored pairs that are compared |

**Calls** (apply to `--stage calls` and `final-calls`)

| option | default | meaning |
| --- | --- | --- |
| `--min-linker-size` | `30` | minimum bases between neighbouring final calls; the exclusion distance is this plus the 147-base nucleosome footprint |

**Resources and diagnostics**

| option | default | meaning |
| --- | --- | --- |
| `--max-memory` | `8G` | upper bound on memory used to buffer scores; accepts `K`, `M`, `G` suffixes. Different budgets can change the last digit of some values |
| `-v`, `--verbose` | off | print the memory plan, index use and record counts to stderr |
| `-m`, `--matrices` | none | reserved for custom matrices; not implemented yet, and warns if given |

Options that do not apply to the selected stage produce a warning rather than being
silently accepted.

## Exit status

`0` on success. Any error prints `error: …` to stderr and exits `1`; nothing is written to
the output on failure to resolve the format or read the input.

## Examples

```sh
# curvature as bigWig, using chr21.fa.fai if present
symcurve chr21.fa chr21.bw

# symmetry scores as text, with the activated roll matrix
symcurve --stage symmetry --roll active chr21.fa chr21.bg

# non-overlapping nucleosome calls under a 2 GB memory budget
symcurve --stage final-calls --max-memory 2G chr21.fa chr21.gff
```
