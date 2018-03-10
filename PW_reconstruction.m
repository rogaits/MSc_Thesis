clear; close all; clc;
Lx = 10; Ly = 10; Lz = 3;
res = 0.01;
rcv_x = 5; rcv_y = 5; rcv_z = 1.5;
[x,y,z] = meshgrid((0:res:Lx)-rcv_x,fliplr(0:res:Ly)-rcv_y,(0:res:Lz)-rcv_z);
f = 344;
c = 344;
k = 2*pi*f/c;
theta_k = pi/2;
phi_k = pi/2;
[kx,ky,kz] = sph2cart_SR(theta_k,phi_k,k);
p = exp(1i*(kx*x + ky*y + kz*z));
%%
surface(x(:,:,1),y(:,:,1),real(p(:,:,1)),'EdgeColor','None')
view(0,90)
colormap(jet)
colorbar
axis equal

%% Plane wave reconstruction using Spherical harmonics
N = 8;
% Angles of arrival
[theta_rec, phi_rec,r_rec] = cart2sph_SR(x(:),y(:),z(:));
sum_m = 0;
sum_n = 0;

for n = 0:N
    sum_m = 0;
    for m = -n:n
        Ynm_rec = Ynm(n, m, theta_rec, phi_rec);
        Ynm_src = conj(Ynm(n, m, theta_k, phi_k));
        sum_m = sum_m + Ynm_rec.*Ynm_src;
    end
    sum_n = sum_n + 4*pi*(1i)^n*sph_bessel(n,k.*r_rec).*sum_m;
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
%% Error graph
Error = abs(p - p_sp).^2/(sum(abs(p).^2));


%% Animation
tf = 0.01;
t=0:0.0001:tf;
fig = figure;
colormap(jet)
M(length(t)) = struct('cdata',[],'colormap',[]);
filename = ['PW_N_' num2str(N) '_f_' num2str(f) '.avi'];
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
