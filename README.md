# algorithmic-resolution
Matlab code for estimating algorithmic resolution limits 

E.A.K. Cohen and N.M. Adams, Dept of Mathematics, Imperial College London

This piece of code estimates the algorithmic resolution limit by
searching for a changepoint in the pair correlation function. The change
point is interpreted to be the radius of correlation rho, and the
algorithmic resolution limit is defined as alpha := rho/2.
This is followed by an F-test at 5% level to ensure it is a genuine radius of
correlation/algorithmic resolution limit.

INPUTS:
 
 r              vector of radial distances at which pair correlation is evaluated
 
 g              vector of pair correlation values at radial distances in r. Must be
                same dimension as r 
                buffer         number of points at either end of g which is are not
                considered as a possible value of rho, the radius of correlatiojn


OUTPUTS:

alpha          algorithmic resolution limit
