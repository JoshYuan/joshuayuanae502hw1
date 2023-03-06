function [V1,V2] = LambertU_V1_V2_from_R1_R2_t(R1, R2, t, grade)
% this function solves Lambert's problem
% to find [V1, V2] from (R1, R2, t)
%{
inputs:
    R1, R2 - position vectors 1 and 2 (km)
    t      - time (s)
    
outputs:
    V1, V2 - velocity vectors 1 and 2 (km)
requirements:
    mu - the gravitational parameter (km^3/s^2) 
    function stumpC(z)
    function stumpS(z)
%}
global mu % = 1.327*10^11 km3/s2 for sun

r1 = norm(R1);
r2 = norm(R2);

r1xr2 = cross(R1,R2);
theta = acos(dot(R1,R2)/r1/r2);

%...Determine whether the orbit is prograde or retrograde:
if nargin < 4 || (wstrcmp(grade,'retro') & (wstrcmp(grade,'pro')))
    grade = 'pro';
    %fprintf('\n ** Prograde trajectory assumed.\n')
end
if strcmp(grade,'pro')
    if r1xr2(3) <= 0
        theta = 2*pi - theta;
    end
elseif strcmp(grade,'retro')
    if r1xr2(3) >= 0
        theta = 2*pi - theta;
    end
end

%...Equation 5.35:
A = sin(theta)*sqrt(r1*r2/(1 - cos(theta)));

% using z = 0 as starting guess
z = -100;
while F(z,t) < 0
    z = z + 0.1;
end
%...Set an error tolerance and a limit on the number of iterations:
tol = 1.e-8;
nmax = 5e3;
%...Iterate on Equation 5.45 until z is determined to within the
%...error tolerance:
ratio = 1;
n = 0;
while (abs(ratio) > tol) & (n <= nmax)
    n = n + 1;
    ratio = F(z,t)/dFdz(z);
    z = z - ratio;
end
%...Report if the maximum number of iterations is exceeded:
if n >= nmax
    fprintf('**Number of iterations exceeds %g, z=%g\n ',nmax,z)
end

%...Equation 5.46a:
f = 1 - y(z)/r1;
%...Equation 5.46b:
g = A*sqrt(y(z)/mu);
%...Equation 5.46d:
gdot = 1 - y(z)/r2;
%...Equation 5.28:
V1 = 1/g*(R2 - f*R1);
%...Equation 5.29:
V2 = 1/g*(gdot*R2 - R1);

return

%...Equation 5.38:
    function temp = y(z)
        temp = r1 + r2 + A*(z*S(z) - 1)/sqrt(C(z));
    end
%...Equation 5.40:
    function temp = F(z,t)
        temp = (y(z)/C(z))^1.5*S(z) + A*sqrt(y(z)) - sqrt(mu)*t;
    end
%...Equation 5.43:
    function temp = dFdz(z)
        if z == 0
            temp = sqrt(2)/40*y(0)^1.5 + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0)));
        else
            temp = (y(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z)) ...
                + 3*S(z)^2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(y(z)) ...
                + A*sqrt(C(z)/y(z)));
        end
    end
%...Stumpff functions:
    function temp = C(z)
        temp = stumpC(z);
    end
    function temp = S(z)
        temp = stumpS(z);
    end
end %lambert
