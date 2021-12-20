function [intRadPat] = integrateRadiationPattern(phi)
% phi (double): angle (radians) subtended between local quantization axis and xy-plane of lab coordinate system

% This function computes the fraction of spontaneously emitted photons our 1:1 optical relay imaging 
% system collects during laser-induced sheet fluorescence (LIF) of Sr+ ions in a UCNP. The fraction 
% depends on the angular momentum (q) of spontaneously emitted photon. We use the D1 line of Sr+ for
% LIF. In this case, the emitted photons can have angular momentum q = m_g - m_e = 0, +/-1, where 
% m_g and m_e are the projection quantum numbers for the ground and excited states, respectively,
% corresponding to photon decay.

% if 'r' not specified, assume ion located away from origin in xy-plane of lab coordinate system
if nargin == 0
    phi = 0;
end

%% Integrate Spontaneous Emission Angular Distribution

% The angular distribution of spontaneous emission is defined in terms of spherical
% coordinates. In quantum mechanics, the z axis is typically taken to be the
% quantization axis, which lies in the imaging plane. Thus, the inherent symmetry of 
% spontaneous emission requires that the angular distribution is independent of the 
% azimuthal angle (in this case 'y' variable in the code). It is additionally independent 
% of the radius (as required by conservation of energy). The 'x' variable in the code represents the 
% angle subtended from imaging plane.

% Note that the angular distribution is the power radiated per solid angle. The angular 
% distributions are integrated over the relevant solid angle as discussed in the 11/07/19 OneNote 
% entry. The area that we image the plasma is much smaller than the diameter of the lenses used in
% the optical relay system, so we treat the plasma as a point source (i.e., the solid angle we 
% integrate over is independent of position within the imaging plane).

intRadPat.dm = [0 -1 1];    % dm = m_g - m_e, change in projection quantum number due to spontaneous decay
intRadPat.fun{1} = @(x,y) sin(x-phi).^2;        % angular distribution for dm = 0
intRadPat.fun{2} = @(x,y) (1+cos(x-phi).^2)./2; % angular distribution for dm = -1
intRadPat.fun{3} = @(x,y) (1+cos(x-phi).^2)./2; % angular distribution for dm = 1

for i = 1:length(intRadPat.fun)    % for each distribution function
    
    % define an anonymous function for the distribution function. note that the solid angle differential is dOm = sin(x)dxdy, 
    % so that's why there is an extra factor of sin(x) below.
    fun = intRadPat.fun{i};
    funForInt = @(x,y) fun(x,y).*sin(x);
    
    numPoints = 301;    % number of points to use in integration w.r.t. x
    
    % define integration range and vector for x - should be +/- 8 degrees, as indicated on 11/07/19
    xmin = pi/2-8*pi/180; % x = pi/2 is the angle that points directly towards the optical relay system
    xmax = pi/2+8*pi/180;
    xForInt = linspace(xmin,xmax,numPoints);
    
    % integration range for y, which is integrated over all y
    ymin = 0;
    ymax = 2*pi;
    yForInt = linspace(ymin,ymax,numPoints);

    [xMeshForInt,yMeshForInt] = meshgrid(xForInt,yForInt);
        
    % do 2D simpson integration - the factor of 8.36921 is the normalization factor that comes from
    % integrating the angular distributions over all solid angle. In this way, the integrated quantity represents
    % the captured fraction.
    wx = getSimpsonWeights(xForInt);
    wy = getSimpsonWeights(yForInt);
    intRadPat.capture(i) = wy'*funForInt(xMeshForInt,yMeshForInt)*wx/8.36921;
end

end