
function [x,y,z] = sph2cart_SR(theta,phi,r)

x = r.*sin(theta).*cos(phi);
y = r.*sin(theta).*sin(phi);
z = r.*cos(theta);