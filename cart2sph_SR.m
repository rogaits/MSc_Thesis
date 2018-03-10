function [theta, phi, r] = cart2sph_SR(x,y,z)
r = sqrt(x.^2+y.^2+z.^2);
theta = atan2(sqrt(x.^2+y.^2),z);
phi = atan2(y,x);