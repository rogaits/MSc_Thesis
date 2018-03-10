%Sound pressure around a rigid sphere
clear all; clc;
N = 8;
Lx = 10; Ly = 10; Lz = 6;
res = 0.1;
rcv_x = 5; rcv_y = 5; rcv_z = 3;
[x,y,z] = meshgrid((0:res:Lx)-rcv_x,fliplr(0:res:Ly)-rcv_y,(0:res:Lz)-rcv_z);
% Angles of arrival
theta_k = pi/2;
phi_k = pi/9;
f = 63;
c = 344;
k = 2*pi*f/c;
%k = 1;
[x_k,y_k,z_k] = sph2cart_SR(theta_k,phi_k,k);
[theta_rec, phi_rec,r_rec] = cart2sph_SR(x(:),y(:),z(:));
sum_m = 0;
sum_n = 0;
ra = 1;
p_i = 0;
p_s = 0;
for n = 0:N
    sum_m = 0;
    for m = -n:n
        Ynm_rec = Ynm(n, m, theta_rec, phi_rec);
        Anm = conj(Ynm(n, m, theta_k, phi_k));
        sum_m = sum_m + Anm .* Ynm_rec;
    end
    jn_r = sph_bessel(n,k.*r_rec);
    hn_r = sph_hankel(n,k.*r_rec);
    %Spherical bessel and hankel functions first derivatives
    jn_ra = (n*sph_bessel(n-1,k*ra) - (n+1)*sph_bessel(n+1,k*ra))/(2*n + 1); 
    hn_ra = (n*sph_hankel(n-1,k*ra) - (n+1)*sph_hankel(n+1,k*ra))/(2*n + 1);
    Bnm = 4.*pi.*(1i)^n.*(jn_r - jn_ra./hn_ra .* hn_r);
    sum_n = sum_n + sum_m .* Bnm;
    p_i = p_i + 4.*pi.*(1i)^n.*sum_m .* jn_r;
    p_s = p_s + 4.*pi.*(1i)^n.*sum_m .* (-jn_ra./hn_ra * hn_r);
end
%%besselj(n,k.*r_rec)
p_sp = reshape(sum_n,size(x));
p_sp(find(r_rec<ra))=NaN;

p_i = reshape(p_i,size(x));
p_i(find(r_rec<ra))=NaN;

p_s = reshape(p_s,size(x));
p_s(find(r_rec<ra))=NaN;
%%
%load('data_pi_ps.mat')
figure
subplot(131)
surface(x(:,:,30),y(:,:,30),real(p_i(:,:,30)),'EdgeColor','None')
surface(x(find(r_rec<ra)),y(find(r_rec<ra)),p_sp(find(r_rec<ra)))
title('Incident Pressure')
xlim([-5 5])
zlim([-5 5])
colormap(jet)
colorbar
caxis([-1 1])
axis equal
view(0,90)

subplot(132)
surface(x(:,:,30),y(:,:,30),real(p_s(:,:,30)),'EdgeColor','None')
surface(x(find(r_rec<ra)),y(find(r_rec<ra)),p_sp(find(r_rec<ra)))
title('Scattered Pressure')
xlim([-5 5])
zlim([-5 5])
colormap(jet)
colorbar
caxis([-1 1])
axis equal
view(0,90)

subplot(133)
surface(x(:,:,30),y(:,:,30),real(p_sp(:,:,30)),'EdgeColor','None')
surface(x(find(r_rec<ra)),y(find(r_rec<ra)),p_sp(find(r_rec<ra)))
title('Total Pressure')
xlim([-5 5])
zlim([-5 5])
colormap(jet)
colorbar
caxis([-1 1])
axis equal
view(0,90)
%% Animation
tf = 0.01;
t=0:0.0001:tf;
fig = figure;
colormap(jet)
M(length(t)) = struct('cdata',[],'colormap',[]);
filename = ['Rigid_sph_N_' num2str(N) '_f_' num2str(f) '_ra_' num2str(ra) '.avi'];
v = VideoWriter(filename);
open(v)
for n=1:length(t)
    p_t = real(exp(1i*2*p_i*f*t(n))*p_sp);
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