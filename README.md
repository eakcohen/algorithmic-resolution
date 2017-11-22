# algorithmic-resolution
Matlab code for estimating algorithmic resolution limits 

E.A.K. Cohen and N.M. Adams, Dept of Mathematics, Imperial College London

est_alpha.m estimates the algorithmic resolution limit by
searching for a changepoint in the pair correlation function. The change
point is interpreted to be the radius of correlation rho, and the
algorithmic resolution limit is defined as alpha := rho/2.
This is followed by an F-test at 5% level to ensure it is a genuine radius of
correlation/algorithmic resolution limit.

bootstrap_alpha.m performs a bootstrap to provide bootstrap intervals for the algorithmic resolution limit estimate.

The individual files have instructions on how to use them
