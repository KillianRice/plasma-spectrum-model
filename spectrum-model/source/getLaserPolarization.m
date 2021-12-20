function [eps,epsM,eps0,epsP] = getLaserPolarization(pol,theta,phi,b,B)
% pol (string): LIF-laser polarization - 'lin','left','right'
% theta (double): angle (radians) between local magnetic field and y-axis of lab frame
% phi (double): angle (radians) between local field and x-y plane
% b (struct): contains information about basis states - see 'defineBasisStates.m'
% B (double): magnitude of local magnetic field in Tesla

%% Error Checking
% ensure 'pol' is specified correctly
if max(strcmp(pol,{'left','right','lin'})) ~= 1
    error('<pol> specified incorrectly. String options are ''lin'', ''left'', ''right''')
end

%% Function
% the following equations will be written up in Grant Gorman's thesis. For the phi = 0 equations, see the magnetized LIF
% paper.

% define polarization in spherical tensor notation for laser polarization vector in local field frame
i = sqrt(-1);
e0 = @(eta,xi,theta,phi) eta*cos(theta)*cos(phi)+exp(i*xi)*sqrt(1-eta^2)*sin(phi); % dm = 0
eP = @(eta,xi,theta,phi) 1/sqrt(2)*(-exp(i*xi)*sqrt(1-eta^2)*cos(phi)-i*eta*sin(theta)+eta*cos(theta)*sin(phi)); % dm = +1
eM = @(eta,xi,theta,phi) 1/sqrt(2)*(+exp(i*xi)*sqrt(1-eta^2)*cos(phi)-i*eta*sin(theta)-eta*cos(theta)*sin(phi)); % dm = -1

% based on user input of 'pol' and quantization axis, choose 'eta' and 'xi' 
if strcmp(pol,'right') && B~= 0   % for right-handed circular polarization, y = -pi/2
    eta = 1/sqrt(2);
    xi = -pi/2;
    eps0 = e0(eta,xi,theta,phi);
    epsP = eP(eta,xi,theta,phi);
    epsM = eM(eta,xi,theta,phi);
elseif strcmp(pol,'left')  && B~= 0   % for left-handed circular polarization, y = pi/2
    eta = 1/sqrt(2);
    xi = +pi/2;
    eps0 = e0(eta,xi,theta,phi);
    epsP = eP(eta,xi,theta,phi);
    epsM = eM(eta,xi,theta,phi);
elseif strcmp(pol,'lin')  && B~= 0 % linear laser polarization
    eta = 1;
    xi = 0;
    eps0 = e0(eta,xi,theta,phi);
    epsP = eP(eta,xi,theta,phi);
    epsM = eM(eta,xi,theta,phi);
elseif strcmp(pol,'lin') && B == 0 % laser sets quantization axis in absense of local fields
    epsP = 0;
    epsM = 0;
    eps0 = 1;
elseif strcmp(pol,'left') && B == 0 % laser sets quantization axis in absense of local fields: theta = pi/2
    eps0 = 0;
    epsM = 1;
    epsP = 0;
elseif strcmp(pol,'right') && B == 0 % laser sets quantization axis in absense of local fields: theta = pi/2
    eps0 = 0;
    epsM = 0;
    epsP = 1;
end

if isnan(eps0*epsP*epsM)
    error('Spherical tensor components of LIF-laser polarization are not valid. Check user input.')
end

% define eps(i,j) the component of laser polarization aligned with transition |i> to |j>
eps = zeros(length(b.s));
for i = b.exc % final state
    for j = b.gnd % initial state
        dm = b.m(i) - b.m(j);
        if i ~= j
            if dm == 0
                eps(i,j) = eps0;
            elseif dm == 1
                eps(i,j) = epsP;
            elseif dm == -1
                eps(i,j) = epsM;
            end
            eps(j,i) = eps(i,j);
        end
    end
end


end