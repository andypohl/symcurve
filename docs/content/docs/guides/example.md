---
title: "Example"
description: "Walk through an example of the calculation."
summary: ""
date: 2023-09-07T16:04:48+02:00
lastmod: 2023-09-07T16:04:48+02:00
draft: false
menu:
  docs:
    parent: ""
    identifier: "example-3b2c4d5e6f7890a1b2c3d4e5f67890ab"
weight: 820
toc: true
seo:
  title: "" # custom title (optional)
  description: "" # custom description (recommended)
  canonical: "" # custom canonical URL (optional)
  noindex: false # false (default) or true
---

### Sequence

Consider the 500 nucleotides of sequences below.  It's enough to wrap at least two eukaryotic nucleosomes.

```
CCAACATTTTGACTTTTTGGGAGGGCACTAGCACCTATCTACCCTGAATC
TAAGTGCTAACAGGAAAGGATGCCAGATTGCATGCCTGCTGATAAAGCCA
CAGTTTGGACTGTCACTCAATCACCATCGTTCCTCCTGTGACTCAGTATA
ACAAGATTGGGAGAATACTCTACAGTTCCTGATTCCCCCACAGAAGTGAA
AGAGAATACTGCAACACACATTCAGTGTTCAGACTTCAGGGGAGTTGCCA
AGTGATTGTTTTCTGTCTTGCCTGAATTTAAGCGTTAACAGGAAAGCTTT
CCAGGTTGGGGATATTAAGAATAAATGAGTTGAGGAAGTTTGGGTTAGCA
TACATTCACTTATCATATCATCCTTCCCTGGATCAATATGCAATGAGTGG
GACAAAACCCTCAACTCCTGGCTTCCCCTTGGAGAAGGAAAAAGCTGGAG
TGTGCATCCAGAATTTCAACTTTTCCCAGTCAGCCTGATGGACTGTTTTC
```

This sequence will be \(S\) with \(n=500\).  Now let's consider \(W\), with \(w_1 = \text{CCA}\), \(w_2 = \text{CAA}\), and so on until \(w_{498} = \text{TTC}\).

### Parameters

For this example, we'll use a typical set of parameters:
* \(a = 5\)
* \(b = 15\)
* \(c = 51\)
* \(\lambda = 0.33335\)
* [\(\rho = \rho^\beta\)](/symcurve/docs/reference/matrices/#roll)

(*To be continued...*)
