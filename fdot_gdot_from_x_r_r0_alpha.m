function [fdot, gdot] = fdot_gdot_from_x_r_r0_alpha(x, r, r0, alpha)
% this function calculates [fdot, gdot] from (x, r, r0, alpha)
%{
inputs: 
    x     - the universal anomaly after time t (km^0.5)
    r     - the radial position after time t (km)
    r0    - the radial position at time to (km)
    alpha - reciprocal of the semimajor axis (1/km)
output: 
    fdot - time derivative of the Lagrange f coefficient (1/s)
    gdot - time derivative of the Lagrange g coefficient (dimensionless)
requirements:
    mu - the gravitational parameter (km^3/s^2) 
    function stumpC(z)
    function stumpS(z)
%}
global mu
z = alpha*x^2;
fdot = sqrt(mu)/r/r0*(z*stumpS(z) - 1)*x;
gdot = 1 - x^2/r*stumpC(z);