---
title: "Algorithm"
description: "Explains what is calculated and how."
summary: ""
date: 2023-09-07T16:04:48+02:00
lastmod: 2023-09-07T16:04:48+02:00
draft: false
menu:
  docs:
    parent: ""
    identifier: "algorithm-6a1a6be4373e933280d78ea53de6158e"
weight: 810
toc: true
seo:
  title: "" # custom title (optional)
  description: "" # custom description (recommended)
  canonical: "" # custom canonical URL (optional)
  noindex: false # false (default) or true
---

### DNA step parameters
The diagram below shows the six main ways that two consecutive base pairs are distanced from each other.

<div id="dna-img-buffer" style="padding-bottom: 20px">
<img src="/symcurve/images/dna-curve-light.png" class="themed" alt="Alt text for the image" title="Title of the image" width="60%" height="60%" />
</div>

Of these six, the three parameters SymCurve uses are \(\Omega\), \(\rho\), and \(\tau\). Also, SymCurve uses estimates of the 3-mer version of these parameters instead of the 2-mer representations in the diagram.  

### 3-mer windowing
In our equations later, we'll refer to the input nucleotide sequence as \(S\), with length \(n\).  Individual nucleotides \(s_i\) form \(S\) as in the notation below:

\[S = \left(s_1,s_2,...,s_{n-1},s_n\right)\]

Then we'll define \(W\) as the set of all sliding-window 3-mer subsequences \(w_i\) of \(S\) such that:

\[
\begin{aligned}
W &= \left(w_1,w_2,...,w_{n-3},w_{n-2}\right) \\
W &= \left[\left(s_1,s_2,s_3\right),\left(s_2,s_3,s_3\right),...,\left(s_{n-3},s_{n-2},s_{n-1}\right),\left(s_{n-2},s_{n-1},s_n\right)\right]
\end{aligned}
\]

As additional shorthand we'll define the the lookup of a given 3-mer \(w_i\) in the \(4\times 4 \times 4\) matrices \(\Omega\), \(\rho\), and \(\tau\) as \(\Omega_i\), \(\rho_i\), and \(\tau_i\) respectively.  The matrix values at each window are used to compute variables \(T_i\) (the cumulative twist-sum), and the deltas \(dx_i\) and \(dy_i\), also for each window:

\[
\begin{aligned}
T_i  &= \sum_{j=1}^{i}\Omega_j \\
dx_i &= \rho_i \sin(T_i) + \tau_i \sin(T_i - \pi/2) \\
dy_i &= \rho_i \cos(T_i) + \tau_i \cos(T_i - \pi/2)
\end{aligned}
\]

### Mapping into 2D space
Across the space of \(n-2\) 3-mers \(w_i\), we'll define \(x_i\) and \(y_i\) coordinates as the sum of the previous coordinates/deltas as:
\[
\begin{aligned}
x_{i+1} &= x_i + dx_i \\
y_{i+1} &= y_i + dy_i
\end{aligned}
\]
where \(x_1 = y_1 = 0\). Note: the range of valid coordinates \(i\) extends to \(n-1\), which is one past the number of 3-mer windows. For this reason,
\(x_1\) and \(y_1\) are ignored in subsequent steps.
### Rolling coordinate averages
We'll now define parameters \(a\) and \(b\), where \(2b+3\) is a sliding window size over the range of coordinates \(a \lt i \lt n-a\). Usually we set \(a=6\) and \(b=4\), but here they're made variable for the sake of illustration and can be set with alternate values as long as \(a-b=2\) and \(2(a+b+1) \ll n\).
The rolling averages \(\overline{x}\) and \(\overline{y}\) are also slightly weighted centrally,
with the values at either end of the window only contributing half what the central values contribute.

\[
\begin{aligned}
\overline{x}_i &= \left(\frac{x_{i-b-1}+x_{i+b+1}}{2} + \sum_{j=i-b}^{i+b}x_j\right)\left(\frac{1}{2b+2}\right) \\
\overline{y}_i &= \left(\frac{y_{i-b-1}+y_{i+b+1}}{2} + \sum_{j=i-b}^{i+b}y_j\right)\left(\frac{1}{2b+2}\right) \\
\end{aligned}
\]

### Curvature
Using the rolling averages, the curvature values \(\kappa_i\) are now possible to calculate over a range of
\(a+c < i < n-a-c\) where \(c\), is another span, usually set to 15.

\[
\kappa_i = \sqrt{(\overline{x}_{i+c}-\overline{x}_{i-c})^2 + (\overline{y}_{i+c}-\overline{y}_{i-c})^2}
\]