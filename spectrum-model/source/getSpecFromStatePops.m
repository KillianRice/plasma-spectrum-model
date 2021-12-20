function [S,dSdt] = getSpecFromStatePops(t,p,radpat,b)
% t (1xn double): excitation time with units gam422^-1
% p (kxn double): ion state population in state |k>

% C(i,j) is the fraction of spontaneously-emitted light that is captured by the optical relay
% system for transition |i> to |j>, note C(i,j) = C(j,i)
% gam(i,j) is the decay rate with units gam422 for the corresponding transition |i> to |j>
C = zeros(b.numstates);
gam = zeros(b.numstates);
for k = 1:length(b.exc)
    exc = b.exc(k);
    gnd = b.gnd(k);
    dm = abs(b.m(exc)-b.m(gnd));
    C(exc,gnd) = radpat.capture(dm == radpat.dm);
    gam(exc,gnd) = b.gam(k);
end

% calculate 'dSdt', the rate of LIF signal collected per unit time due to fluorescence from each
% excited state
dSdt = zeros(size(t));
for k = 1:length(b.exc) % loop over spontaneous emission decay path
    exc = b.exc(k);
    gnd = b.gnd(k);
    dSdt = dSdt + p(exc,:).*gam(exc,gnd).*C(exc,gnd);
end

% integrate the signal w.r.t. time
S = trapz(t,dSdt);


end
