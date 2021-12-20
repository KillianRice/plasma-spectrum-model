function [b] = defineBasisStates()
% This function defines the relevant quantum state information for the 2S_{1/2} - 2P_{1/2} transition
% of Sr+ used for LIF imaging.

% define quantum numbers of basis states coupled by LIF laser
b.l = [0 0 1 1];                    % orbital angular momentum
b.s = [1/2 1/2 1/2 1/2];            % spin angular momentum
b.j = [1/2 1/2 1/2 1/2];            % total angular momentum
b.m = [-1/2 1/2 -1/2 1/2];          % projection of total angular momentum

b.upInd = [3 4];
b.lowInd = [1 2];

% assign column vectors for each basis state
b.numstates = 4;
b.psi = cell(1,b.numstates);
I = eye(b.numstates);
for i = 1:b.numstates
    b.psi{i} = I(:,i);
end

% calculate 3-j symbols - note that J(i,j) = J(j,i)
% calculate psi^{dagger}*psi for each state
b.J = zeros(b.numstates);
for i = 1:b.numstates % final state
    for j = 1:b.numstates % initial state
        b.J(i,j) = Wigner3j([b.j(i) 1 b.j(j)],[-b.m(i) b.m(i)-b.m(j) +b.m(j)]);
        b.bbd{i,j} = complex(b.psi{i}*b.psi{j}',zeros(4,4));
    end
end

% define possible spontaneous decay paths 
b.exc = [3 3 4 4];  % excited state for each decay path
b.gnd = [1 2 1 2];  % ground state for each decay path
b.gamMan = [1 1 1 1];  % manifold decay rate for each path with units gam422

% define decay rates (with units \gam_422) for each transition)
% calculate other useful quantities for each transition, such as the jump operators (c) and
% multiplications of the decay operators (c^{dagger}*c)
b.gam = zeros(size(b.gamMan));
for i = 1:length(b.gamMan)
    b.gam(i) = (2*b.j(b.exc(i))+1)*b.J(b.exc(i),b.gnd(i))^2*b.gamMan(i);
    b.c{i} = b.psi{b.gnd(i)}*b.psi{b.exc(i)}';
    b.cdc{i} = b.c{i}'*b.c{i};
end

end