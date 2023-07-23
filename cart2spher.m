function [r, theta, phi] = cart2spher(x, y, z)

r = sqrt(x^2 + y^2 + z^2);
theta = acos(z/r);
phi = atan(y/x);
