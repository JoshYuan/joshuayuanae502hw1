function x = keplerU_x_from_dt_r0_vr0_alpha(dt, r0, vr0, alpha)
% this function applies Newtons method in Kepler equation
% to find universal anomaly x from (dt, r0, vr0, alpha)
%{
inputs:
    dt      - time since x = 0 (s)
    r0      - radial position @x=0 (km) 
    vr0     - radial velocity @x=0 (km/s)
    alpha   - reciprocal of semimajor axis (1/km)
outputs:
    x   - the universal anomaly (km^.5)
requirements:
    mu - the gravitational parameter (km^3/s^2) 
    function stumpC(z)
    function stumpS(z)
%}
global mu

%...Set an error tolerance and a limit on the number of iterations:
error = 1.0e-8;
nMax = 1.0e3;
%...Starting value for x:
x = sqrt(mu)*abs(alpha)*dt;

n = 0;
F_dFdx = 1;
while abs(F_dFdx) > error && n <= nMax
    n = n + 1;
    C = stumpC(alpha*x.^2);
    S = stumpS(alpha*x.^2);
    F = r0*vr0/sqrt(mu)*x.^2*C + (1 - alpha*r0)*x.^3*S + r0*x - sqrt(mu)*dt;
    dFdx = r0*vr0/sqrt(mu)*x*(1 - alpha*x.^2*S) + (1 - alpha*r0)*x.^2*C + r0;
    F_dFdx = F/dFdx;
    x = x - F_dFdx;
end
%...Deliver a value for x, but report that nMax was reached:
if n > nMax
    %fprintf('**No. iterations of Keplers equation = %g', n)
    fprintf('Max iterations reached, F/dFdx = %g\n', F_dFdx)
end