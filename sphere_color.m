theta = 0:0.01:pi;
phi = 0:0.01:2*pi;
f = sin(theta.').^2*cos(2*phi);
[THETA, PHI] = meshgrid(theta,phi);


F = sin(THETA).^2.*cos(2*PHI);
R = abs(F); %distance from the origin
x = sin(THETA).*cos(PHI);
y = sin(THETA).*sin(PHI);
z = cos(THETA);

subplot(121)
surf(x,y,z,F,'EdgeColor','None')
title('Function over unit sphere','Interpreter','Latex')
xlabel('x', 'Interpreter','Latex')
ylabel('y', 'Interpreter','Latex')
zlabel('z', 'Interpreter','Latex')
%shading interp
colormap(jet)
colorbar('ylim',[-1 1])
axis equal

subplot(122)
surf(R.*x,R.*y,R.*z,F,'EdgeColor','None')
%shading interp
title('Balloon plot','Interpreter','Latex')
xlabel('x', 'Interpreter','Latex')
ylabel('y', 'Interpreter','Latex')
zlabel('z', 'Interpreter','Latex')
colormap(jet)
colorbar('ylim',[-1 1])
axis equal

figure
contourf(phi,theta,f)
title('Contour map on $\theta \phi$ plane','Interpreter','Latex')
xlabel('$\phi$ (radians)', 'Interpreter','Latex')
ylabel('$\theta$ (radians)', 'Interpreter','Latex')
colormap(jet)
colorbar('ylim',[-1 1])


