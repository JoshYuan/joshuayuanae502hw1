function [f, g] = f_g_from_x_t_r0_alpha(x, t, r0, alpha)
% this function calculates [f, g] from (x, t, r0, alpha)
%{
inputs: 
    x     - the universal anomaly after time t (km^0.5)
    t     - the time elapsed since r0 (s)
    r     - the radial position after time t (km)
    r0    - the radial position at time t0 (km)
    alpha - reciprocal of the semimajor axis (1/km)
output: 
    f - the Lagrange f coefficient (dimensionless)
    g - the Lagrange g coefficient (s)
requirements:
    mu - the gravitational parameter (km^3/s^2) 
    function stumpC(z)
    function stumpS(z)
%}
global mu
z = alpha*x^2;
f = 1 - x^2/r0*stumpC(z);
g = t - 1/sqrt(mu)*x^3*stumpS(z);
end
