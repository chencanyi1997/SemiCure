
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semicure

<!-- badges: start -->

<!-- badges: end -->

The goal of semicure is to estimate the parameters given some i.i.d data
generated from `semicure model`.

## Introduction

The population survival function for a with covariates \(\mathbf{Z}\) is
given by \[
S_{pop}(t|\mathbf{Z})=G_{\gamma}\{e^{\beta^T\mathbf{Z}}F(t)\},
\] where \[G_{\gamma}(x)=\left\{
\begin{aligned}
&(1+\gamma x)^{-1/\gamma}&,\quad \gamma>0 \\
&e^{-x}&, \quad \gamma=0.
\end{aligned}
\right.
\] This class of models include the proportional hazards cure model and
the proportional odds cure model as special cases. When \(\gamma=0\), it
reduces to the proportional hazards cure model; when \(\gamma=1\), it
has the proportional odds structure.

## Installation

You can install the released version of semicure from
[Github](https://github.com/chencanyi1997/semicure) with:

``` r
devtools::install_github("chencanyi1997/semicure")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(semicure)
## basic example code
demo(semicuredemo)
```

To see more through sample, you can check the vigenitte accompanied with
the package. But you should install this package with the following code
to build vigenitte. It will take a little time to build
it.

``` r
devtools::install_github("chencanyi1997/semicure", build_vignettes = TRUE)
vignette(semicure)
```
