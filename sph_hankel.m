function hn = sph_hankel(n,x)
%spherical Hankel function of second kind and order n
hn = besselh(n+0.5,2,x).*(pi./(2*x)).^0.5;