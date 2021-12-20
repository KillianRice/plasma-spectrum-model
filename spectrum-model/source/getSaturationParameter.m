function [s0] = getSaturationParameter(Om)
% Om (nxn double): Rabi coupling between LIF basis states with units gam422
    % (i.e., Om(i,j) is the Rabi coupling between states |i> and |j>
% s0 (nxn double): s0(i,j) is the saturation parameter for transition |i> to |j>

s0 = 2.*Om.^2; % note that Om is already dimensionless

end