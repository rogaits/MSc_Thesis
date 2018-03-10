% Analytical Spherical wave propagation
Lx = 10; Ly = 10; Lz = 3;
res = 0.1;
rcv_x = 5; rcv_y = 5; rcv_z = 1.5; %room dimensions
xs = 5; ys = 0; zs = 0; %point source location
[x,y,z] = meshgrid((0:res:Lx)-rcv_x,fliplr(0:res:Ly)-rcv_y,(0:res:Lz)-rcv_z);
f = 250;
c = 344;
k = 2*pi*f/c;
theta_k = pi/2;
phi_k = pi/9;
phi = pi/2;
rho = 1.2;
q = 0.01;
R = sqrt((x - xs).^2 + (y - ys).^2 + (z - zs).^2); % ||r - r_s|| in analytical expression
p =  1i*2*pi*f*rho*q*exp(-1i*k*R)./(4*pi*R);
p = reshape(p,size(x));
%%
figure
surface(x(:,:,1),y(:,:,1),real(p(:,:,1)),'EdgeColor','None')
colormap(jet)
colorbar
caxis([-1 1])
axis equal
%% Plane wave reconstruction using Spherical harmonics
N = 64;
Lx = 10; Ly = 10; Lz = 0;
res = 0.1;
rcv_x = 5; rcv_y = 5; rcv_z = 1.5;
[x,y,z] = meshgrid((0:res:Lx)-rcv_x,fliplr(0:res:Ly)-rcv_y,(0:res:Lz)-rcv_z);
% Angles of point source
[theta_src, phi_src,r_src] = cart2sph_SR(xs,ys,zs);
[theta_rec, phi_rec,r_rec] = cart2sph_SR(x(:),y(:),z(:));
sum_m = 0;
sum_n = 0;

for n = 0:N
    sum_m = 0;
    for m = -n:n
        Ynm_rec = Ynm(n, m, theta_rec, phi_rec);
        Ynm_src = conj(Ynm(n, m, theta_src, phi_src));
        sum_m = sum_m + Ynm_rec.*Ynm_src;
    end
    sum_n = sum_n + 4*pi*(-1i)*k*sph_hankel(n,k*r_src)*sph_bessel(n,k.*r_rec).*sum_m;
end
%%besselj(n,k.*r_rec)
p_sp = reshape(sum_n,size(x));
%%
figure
surface(x(:,:),y(:,:),real(p_sp(:,:)),'EdgeColor','None')
colormap(jet)
colorbar
axis equal
view(0,90)
%% Animation
tf = 0.01;
t=0:0.0001:tf;
fig = figure;
colormap(jet)
M(length(t)) = struct('cdata',[],'colormap',[]);
filename = ['SPHW_N_' num2str(N) '_f_' num2str(f) '.avi'];
v = VideoWriter(filename);
open(v)
for n=1:length(t)
    p_t = real(exp(1i*2*pi*f*t(n))*p_sp);
    surf(x(:,:),y(:,:),p_t,'EdgeColor','None')
    view(0,90)
    colorbar
    caxis([-1 1])
    drawnow
    %pause(0.001)
    M(n) = getframe;
    writeVideo(v,M(n))
end
close(v)