function [a,e,i,Omega,omega,f] = orbital_elements_from_R_V(R,V)
%find orbital elements from r and v vectors
%{
a = semimajor axis
e = eccentricity
i = inclination
Omega = longitude of ascending node
omega = argument of periapse
f = true anomaly
%}

global mu

r = norm(R);
v = norm(V);
a = r/(2-r*v*v/mu);
H = cross(R,V);
h = norm(H);
E = 1/mu* (cross(V,H)-mu*R/r);
e = norm(E);
N = cross([0,0,1],H/h);
n = norm(N);
i = acos(H(3)/h);

Omega = acos(N(1)/n);
if N(2) < 0
    Omega = 2*pi-Omega;
end

omega = dot(N,E)/(n*e);
if E(3) < 0
    omega = 2*pi-omega;
end

f = acos(dot(E,R)/(e*r));
if dot(R,V) < 0
    f = 2*pi - f;
end

