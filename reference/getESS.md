# getESS

This function calculates the effective sample size for every dose group
via
[`RBesT::ess()`](https://opensource.nibr.com/RBesT/reference/ess.html).

## Usage

``` r
getESS(post_list, method = c("elir", "moment", "morita"), n_digits = 1, ...)
```

## Arguments

- post_list:

  A posterior list object, for which the effective sample size for each
  dose group should be calculated

- method:

  A string specifying the method of ESS calculation, see
  `?RBesT::ess()`.

- n_digits:

  An integer for the number of digits the result should be rounded to.

- ...:

  Optional arguments applicable to specific methods, see
  `?RBesT::ess()`.

## Value

A vector of the effective sample sizes for each dose group
