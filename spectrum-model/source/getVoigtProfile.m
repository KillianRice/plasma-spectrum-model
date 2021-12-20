function [V,V_int] = getVoigtProfile(det,det0,sigG,sigL)
% det (row vector double): vector of frequencies to return Voigt profile for
% det0 (double): center of Gaussian (and Voigt) distribution
% sigG (double): RMS width of Gaussian distribution
% sigL (double): FWHM of Lorentzian distribution
% V (row vector double): Voigt profile (integral normalized to one) corresponding to freq. Note that
                        % V(i) corresponds to det(i).

% Note: det, det0, sigG, sigL must all have the same units, but what that unit is doesn't matter.

if sigG/sigL < .015
    error('sigG/sigL should be > 0.015')
end

% Define normalized Gaussian distribution with center det0
G = @(dets,detc) exp(-(dets-detc).^2./(2*sigG^2))./sqrt(2*pi*sigG^2);

% Define normalized Lorentzian distribution
L = @(dets,detc) sigL/(2*pi)./(sigL^2/4 + (dets-detc).^2);

% Define vector of detunings to use for convolution
pts = 21; % points per linewidth - this was chosen to balance maintaining normalization and cost efficiency
numFWHM = 10; % this was chosen to balance numerical accuracy and cost efficiency 
convpts = round(pts*numFWHM/2)*2+1; % odd number of conv points that scales with width
detConvRange = getVoigtFWHM(sigG,sigL)*numFWHM; 
detConv = linspace(-detConvRange,detConvRange,convpts);

% For each detuning in 'det', convolve G and L to get Voigt profile 'V'
V = zeros(size(det));
norm_fac = trapz(detConv,L(detConv,0));
for i = 1:length(det)
    V(i) = trapz(detConv,G(det(i)-detConv,det0).*L(detConv,0))./norm_fac;
end
V_int = trapz(det,V);

end