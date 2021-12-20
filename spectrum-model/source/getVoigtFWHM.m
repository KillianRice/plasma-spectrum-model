function [sigV] = getVoigtFWHM(sigG,sigL)
% sigG (double): RMS width of Gaussian distribution
% sigL (double): FWHM of Lorentzian
% sigV (double): FWHM of Voigt profile

sigG = sigG*2*sqrt(2*log(2)); % FWHM of Gaussian distribution 
sigV = 0.5346*sigL+sqrt(0.2166*sigL^2+sigG^2); % FWHM of Voigt distribution

end