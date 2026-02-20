# The `BayesianMCPmod` package

Simulate, analyze, and evaluate Bayesian MCPMod trials with normally and
binary distributed endpoints. Bayesian MCPMod [(Fleischer et al.,
2022)](https://doi.org/10.1002/pst.2193) is an innovative method that
improves the traditional MCPMod by systematically incorporating
historical data, such as previous placebo group data. This package
offers functions for simulating, analyzing, and evaluating Bayesian
MCPMod trials with normally and binary distributed endpoints. It enables
the assessment of trial designs incorporating historical data across
various true dose-response relationships and sample sizes. Robust
mixture prior distributions, such as those derived with the
Meta-Analytic-Predictive approach [(Schmidli et al.,
2014)](https://doi.org/doi:10.1111/biom.12242), can be specified for
each dose group. Resulting mixture posterior distributions are used in
the Bayesian Multiple Comparison Procedure and modeling steps. The
modeling step also includes a weighted model averaging approach
[(Pinheiro et al., 2014)](https://doi.org/doi:10.1002/sim.6052).
Estimated dose-response relationships can be bootstrapped and
visualized.

## Installation

Install the package from CRAN using

`{r} install.packages("BayesianMCPMod")`

The development version can be installed from this repository.

`{r} # install.packages("remotes") remotes::install_github("https://github.com/Boehringer-Ingelheim/BayesianMCPmod")`

## Documentation

The package documentation is hosted
[here](https://boehringer-ingelheim.github.io/BayesianMCPMod/).
