function c = stumpC(z)
% Stumpff function C(z)
% z - input
% c - return
if z > 0
    c = (1 - cos(sqrt(z)))/z;
elseif z < 0
    c = (cosh(sqrt(-z)) - 1)/(-z);
else
    c = 1/2;
end