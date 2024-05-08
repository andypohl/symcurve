---
title: "Matrices"
description: ""
summary: ""
date: 2023-09-07T16:12:37+02:00
lastmod: 2023-09-07T16:12:37+02:00
draft: false
menu:
  docs:
    parent: ""
    identifier: "matrices-a1b2c3d4e5f67890a1b2c3d4e5f67890"
weight: 911
toc: true
seo:
  title: "" # custom title (optional)
  description: "" # custom description (recommended)
  canonical: "" # custom canonical URL (optional)
  noindex: false # false (default) or true
---
The default matrices SymCurve uses for roll (\(\rho\)), twist (\(\Omega\)), and tilt (\(\tau\)) are listed on this page. The matrix sizes are all 4x4x4, where the values are indexed by nucleotide 3-mers as in the matrix \(M\) below:
\[
M = 
\begin{bmatrix}
\begin{bmatrix} M_{AAA} & M_{AAT} & M_{AAG} & M_{AAC} \\ M_{ATA} & M_{ATT} & M_{ATG} & M_{ATC} \\ M_{AGA} & M_{AGT} & M_{AGG} & M_{AGC} \\ M_{ACA} & M_{ACT} & M_{ACG} & M_{ACC} \end{bmatrix}, \\
\begin{bmatrix} M_{TAA} & M_{TAT} & M_{TAG} & M_{TAC} \\ M_{TTA} & M_{TTT} & M_{TTG} & M_{TTC} \\ M_{TGA} & M_{TGT} & M_{TGG} & M_{TGC} \\ M_{TCA} & M_{TCT} & M_{TCG} & M_{TCC} \end{bmatrix}, \\
\begin{bmatrix} M_{GAA} & M_{GAT} & M_{GAG} & M_{GAC} \\ M_{GTA} & M_{GTT} & M_{GTG} & M_{GTC} \\ M_{GGA} & M_{GGT} & M_{GGG} & M_{GGC} \\ M_{GCA} & M_{GCT} & M_{GCG} & M_{GCC} \end{bmatrix}, \\
\begin{bmatrix} M_{CAA} & M_{CAT} & M_{CAG} & M_{CAC} \\ M_{CTA} & M_{CTT} & M_{CTG} & M_{CTC} \\ M_{CGA} & M_{CGT} & M_{CGG} & M_{CGC} \\ M_{CCA} & M_{CCT} & M_{CCG} & M_{CCC} \end{bmatrix}
\end{bmatrix}
\]



### Roll

\[
\rho^\psi = 
\begin{bmatrix}
\begin{bmatrix} 0.1 & 0.0 & 4.2 & 1.6 \\ 9.7 & 0.0 & 8.7 & 3.6 \\ 6.5 & 2.0 & 4.7 & 6.3 \\ 5.8 & 2.0 & 5.2 & 5.2 \end{bmatrix}, \\
\begin{bmatrix} 7.3 & 9.7 & 7.8 & 6.4 \\ 7.3 & 0.1 & 6.2 & 5.1 \\ 10.0 & 5.8 & 0.7 & 7.5 \\ 10.0 & 6.5 & 5.8 & 6.2 \end{bmatrix}, \\
\begin{bmatrix} 5.1 & 3.6 & 6.6 & 5.6 \\ 6.4 & 1.6 & 6.8 & 5.6 \\ 6.2 & 5.2 & 5.7 & 8.2 \\ 7.5 & 6.3 & 4.3 & 8.2 \end{bmatrix}, \\
\begin{bmatrix} 6.2 & 8.7 & 9.6 & 6.8 \\ 7.8 & 4.2 & 9.6 & 6.6 \\ 5.8 & 5.2 & 3.0 & 4.3 \\ 0.7 & 4.7 & 3.0 & 5.7 \end{bmatrix}
\end{bmatrix}
\]

\[
\rho^\alpha = 
\begin{bmatrix}
\begin{bmatrix} 0.06330 & 0.35000 & 4.67090 & 2.64115 \\ 6.27340 & 0.35000 & 7.71710 & 4.44325 \\ 4.88840 & 3.92320 & 5.05230 & 6.88290 \\ 5.49030 & 3.92320 & 5.30550 & 5.30550 \end{bmatrix}, \\
\begin{bmatrix} 4.67090 & 6.27340 & 5.00295 & 5.06730 \\ 4.67090 & 0.06330 & 4.76180 & 4.06330 \\ 7.70000 & 5.49030 & 3.05865 & 6.75525 \\ 7.70000 & 4.88840 & 7.07195 & 4.99070 \end{bmatrix}, \\
\begin{bmatrix} 4.06330 & 4.44325 & 5.98060 & 5.51645 \\ 5.06730 & 2.64115 & 6.62555 & 5.51645 \\ 4.99070 & 5.30550 & 5.89135 & 9.08230 \\ 6.75525 & 6.88290 & 5.89135 & 9.08230 \end{bmatrix}, \\
\begin{bmatrix} 4.76180 & 7.71710 & 6.89960 & 6.62555 \\ 5.00295 & 4.67090 & 6.89960 & 5.98060 \\ 7.07195 & 5.30550 & 3.86900 & 5.90000 \\ 3.05865 & 5.05230 & 3.86900 & 5.82700 \end{bmatrix}
\end{bmatrix}
\]

### Twist
\[
\Omega = 
\begin{bmatrix}
\begin{bmatrix} 0.5986 & 0.5986 & 0.5986 & 0.5986 \\ 0.5986 & 0.5986 & 0.5986 & 0.5986 \\ 0.5986 & 0.5986 & 0.5986 & 0.5986 \\ 0.5986 & 0.5986 & 0.5986 & 0.5986 \end{bmatrix}, \\
\begin{bmatrix} 0.5986 & 0.5986 & 0.5986 & 0.5986 \\ 0.5986 & 0.5986 & 0.5986 & 0.5986 \\ 0.5986 & 0.5986 & 0.5986 & 0.5986 \\ 0.5986 & 0.5986 & 0.5986 & 0.5986 \end{bmatrix}, \\
\begin{bmatrix} 0.5986 & 0.5986 & 0.5986 & 0.5986 \\ 0.5986 & 0.5986 & 0.5986 & 0.5986 \\ 0.5986 & 0.5986 & 0.5986 & 0.5986 \\ 0.5986 & 0.5986 & 0.5986 & 0.5986 \end{bmatrix}, \\
\begin{bmatrix} 0.5986 & 0.5986 & 0.5986 & 0.5986 \\ 0.5986 & 0.5986 & 0.5986 & 0.5986 \\ 0.5986 & 0.5986 & 0.5986 & 0.5986 \\ 0.5986 & 0.5986 & 0.5986 & 0.5986 \end{bmatrix}
\end{bmatrix}
\]

### Tilt
Though SymCurve makes room for \(\tau\) in its calculations, it's not used in the default implementation.

\[
\tau = 
\begin{bmatrix}
\begin{bmatrix} 0.0 & 0.0 & 0.0 & 0.0 \\ 0.0 & 0.0 & 0.0 & 0.0 \\ 0.0 & 0.0 & 0.0 & 0.0 \\ 0.0 & 0.0 & 0.0 & 0.0 \end{bmatrix}, \\
\begin{bmatrix} 0.0 & 0.0 & 0.0 & 0.0 \\ 0.0 & 0.0 & 0.0 & 0.0 \\ 0.0 & 0.0 & 0.0 & 0.0 \\ 0.0 & 0.0 & 0.0 & 0.0 \end{bmatrix}, \\
\begin{bmatrix} 0.0 & 0.0 & 0.0 & 0.0 \\ 0.0 & 0.0 & 0.0 & 0.0 \\ 0.0 & 0.0 & 0.0 & 0.0 \\ 0.0 & 0.0 & 0.0 & 0.0 \end{bmatrix}, \\
\begin{bmatrix} 0.0 & 0.0 & 0.0 & 0.0 \\ 0.0 & 0.0 & 0.0 & 0.0 \\ 0.0 & 0.0 & 0.0 & 0.0 \\ 0.0 & 0.0 & 0.0 & 0.0 \end{bmatrix}
\end{bmatrix}
\]
