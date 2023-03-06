function [R,V] = State_R0_V0_t(R0, V0, t)
% this function propagates the state vectors
% [R, V] from (R0, V0, t)
%{
inputs:
    R0 - initial position vector (km)
    V0 - initial velocity vector (km/s)
    t  - elapsed time (s)
outputs:
    R  - final position vector (km)
    V  - final velocity vector (km/s)
requirements
    global mu - gravitational parameter (km^3/s^2)
    function keplerU_x_from_dt_r0_vr0_alpha
    function f_g_from_x_t_r0_alpha
    function fdot_gdot_from_x_t_r0_alpha
%}

global mu   % 1.327e11 for sun (km3/s2)

r0 = norm(R0);  % R0 magnitude
v0 = norm(V0);  % V0 magnitude
vr0 = dot(R0, V0)./r0;   % initial radial velocity
alpha = 2./r0 - v0^2./mu; % alpha = 1/a (a = semimajor axis)
x = keplerU_x_from_dt_r0_vr0_alpha(t, r0, vr0, alpha);  % universal anomaly
[f, g] = f_g_from_x_t_r0_alpha(x, t, r0, alpha);    % f and g
%disp([f,g])
R = f.*R0 + g.*V0;    % final position vector
r = norm(R);       % r magnitude
[fdot, gdot] = fdot_gdot_from_x_r_r0_alpha(x, r, r0, alpha); % fdot and gdot
V = fdot.*R0 + gdot.*V0;  % final velocity