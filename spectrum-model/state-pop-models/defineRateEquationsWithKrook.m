function [dfdt] = defineRateEquationsWithKrook(~,p,R,g,G,w,mu,gamd)
% t (1x1 double): time in units of gam, an unused (but required) input parameter
% p (4nx1 double): column vector of the state population of state |i> for each velocity
% R (4x4xn double): R(i,j,n) defines scattering rate for states |i> and |j> for velocity v(n) in units gam^-1
% g (4x4 double): g(i,j) defines decay rate from state |i> to |j> in units gam^-1
% G (1xn double): Maxwell-Boltzmann distribution for velocity v(n)
% w (nx1 double): w(n) contains the weight for Simpson integration for velocity v(n) with units k/gam
% mu (1x1 double): collision frequency for Krook term with units \gamma^-1

% Unpackage state populations, p, s.t. we obtian p(i,n) the population for state |i> for velocity v(n)
p = reshape(p,size(g,1),[]);

% Calculate current ensemble average of each state population
pEns = zeros(4,1);
for i = 1:4 % iterate through each state
    pEns(i) =  p(i,:)*w;
end

% For each velocity, calculate the rate equations
dfdt = zeros(4,length(G));
for vInd = 1:length(G)
    dfdt(1,vInd) = -R(3,1,vInd)*(p(1,vInd)-p(3,vInd))-R(4,1,vInd)*(p(1,vInd)-p(4,vInd))+g(3,1)*p(3,vInd)+g(4,1)*p(4,vInd)-mu*(p(1,vInd)-pEns(1)*G(vInd));
    dfdt(2,vInd) = -R(4,2,vInd)*(p(2,vInd)-p(4,vInd))-R(3,2,vInd)*(p(2,vInd)-p(3,vInd))+g(3,2)*p(3,vInd)+g(4,2)*p(4,vInd)-mu*(p(2,vInd)-pEns(2)*G(vInd));
    dfdt(3,vInd) = +R(3,1,vInd)*(p(1,vInd)-p(3,vInd))+R(3,2,vInd)*(p(2,vInd)-p(3,vInd))-g(3,1)*p(3,vInd)-g(3,2)*p(3,vInd)-mu*(p(3,vInd)-pEns(3)*G(vInd)) - p(3,vInd)*gamd;
    dfdt(4,vInd) = +R(4,2,vInd)*(p(2,vInd)-p(4,vInd))+R(4,1,vInd)*(p(1,vInd)-p(4,vInd))-g(4,2)*p(4,vInd)-g(4,1)*p(4,vInd)-mu*(p(4,vInd)-pEns(4)*G(vInd)) - p(4,vInd)*gamd;
end

dfdt = dfdt(:);

%% Output parameters
% dfdt (4nx1 double): rate of change of state population, time in units gam^-1

end