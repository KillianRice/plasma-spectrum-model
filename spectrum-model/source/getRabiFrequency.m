function [Om] = getRabiFrequency(Om0,eps,b)
% Om0 (double): sublevel-independent Rabi freq with units gam422
% eps (4x4 complex double): eps(i,j) laser polarization relative to local magnetic field, expressed in
%                           spherical tensor notation that drives |i> to |j>.  see 'source/getLaserPolarization.m'
% b (struct): basis state information, see 'source/defineBasisStates.m'

% Om (4x4 complex double): Om(i,j) is the Rabi coupling betweten state |i> and |j> with units gam422

zeromat = zeros(b.numstates);
Om = complex(zeromat,zeromat);
for i = b.exc
    for j = b.gnd
        % get the largest total angular momentum quantum number between initial and final state
        if b.j(i) > b.j(j)
            j_gt = b.j(i);
        else
            j_gt = b.j(j);
        end
        Om(i,j) = Om0*(-1)^(b.j(i)+b.j(j)+j_gt-b.m(i))*sqrt(2*b.j(i)+1)*conj(eps(i,j))*b.J(i,j);
        Om(j,i) = conj(Om(i,j));
    end
end

end