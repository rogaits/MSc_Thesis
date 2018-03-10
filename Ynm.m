function Y = Ynm(varargin)
% if 4 input arguments provided then
% n: degree of Associated Legendre Function
% m: order of Associated Legendre Function
% theta: measured from the polar axis z 0 - pi
% phi: Measured in the xy plane from the x-axis (azimuthal direction) 0 - 2pi
% Output
% Y: matrix with values of spherical harmonics for m = -n:n every row
% corresponds to a value of m.

% if 4 input arguments provided then
% n: degree of Associated Legendre Function
% m: order of Associated Legendre Function
% theta: measured from the polar axis z 0 - pi
% phi: Measured in the xy plane from the x-axis (azimuthal direction) 0 - 2pi
% Output
% Y: scalar with value of spherical harmonic
%disp('Number of provided inputs: ' + length(varargin))
% Associated Legendre Functions

if nargin == 3
    n = varargin{1,1};
    theta = varargin{1,2};
    phi = varargin{1,3};
    indx = 1;
    Y = zeros(2*n+1,length(theta));
    for m = -n:n
        if m < 0
            M = -1*m;
            Pmn = (-1)^M*factorial(n - M)/factorial(n + M)*legendre(n,cos(theta));
            pos = M;
        else
            Pmn = legendre(n,cos(theta));
            pos = m;
        end
        Y(indx,:) = sqrt((2*n + 1)/(4*pi)*factorial(n - m)/factorial(n + m))*Pmn(pos+1,:).*exp(1i*m*phi);
        indx = indx + 1;
    end
    
elseif nargin == 4
    n = varargin{1,1};
    m = varargin{1,2};
    theta = varargin{1,3};
    phi = varargin{1,4};
    indx = 1;
    Y = zeros(2*n+1,length(theta));
    if m < 0
        M = -1*m;
        Pmn = (-1)^M*factorial(n - M)/factorial(n + M)*legendre(n,cos(theta));
        pos = M;
    else
        Pmn = legendre(n,cos(theta));
        pos = m;
    end
    if m == 0 && n == 0
        Y = sqrt((2*n + 1)/(4*pi)*factorial(n - m)/factorial(n + m))*Pmn.'.*exp(1i*m*phi);
    else
        Y = sqrt((2*n + 1)/(4*pi)*factorial(n - m)/factorial(n + m))*squeeze(Pmn(pos+1,:,:,:)).'.*exp(1i*m*phi);
    end
else
    disp('Not enough arguments provided!')
end


