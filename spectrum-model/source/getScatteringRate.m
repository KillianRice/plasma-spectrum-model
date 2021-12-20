function [R] = getScatteringRate(det,v,dEz,Om,gamL,b)
% det (double): LIF-laser detuning from unperturbed resonance with units \gamma_422
% v (double): hydrodynamic flow velocity with units \gamma_422/k
% dEz (1xn double): Zeeman shifts for each basis state |n> with units \hbar\gamma_422
% Om (nxn complex double): Rabi coupling between basis states with units \gamma_422
% b (struct): basis state information
% R (nxn double): scattering rate with units \gamma_422
% gamL (double): linewidth of Lorentzian laser lineshape with units \gamma_422

% \gamma_422 is the natural linewidth of Sr+ D1 line (i.e., the imaging transition)

%% Define Scattering Rate for Each Transition
% Det(i,j) is the Doppler- and Zeeman-shifted resonance condition with units \gamma_422 for velocity v(k)
% R(i,j) is the scattering rate (with units \gamma_422) from state |j> to state |i> for velocity v(k)
% gam(i,j) is the decay rate for a given transition

% See Castro thesis Appendix A for scattering rate equations.

gam_tot = 1+gamL;
R = complex(zeros(b.numstates,b.numstates),0);
for i = b.exc
    for j = b.gnd
        Det = det+v-(dEz(i)-dEz(j));
        R(i,j) = Om(i,j)*conj(Om(i,j))*gam_tot/(gam_tot^2+4*Det^2);
    end
end

end