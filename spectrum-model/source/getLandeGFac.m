function [g] = getLandeGFac(j,l,s)
% j (1xn double): total angular momentum quantum number
% l (1xn double): orbital angular momentum quantum number
% s (1xn double): spin angular momentum quantum number

%% Output Parameters
% g (1xn double): Lande g-factor

%% Calculate Lande g-factor from Quantum Numbers
g = 1+(j.*(j+1)+s.*(s+1)-l.*(l+1))./(2.*j.*(j+1));

end