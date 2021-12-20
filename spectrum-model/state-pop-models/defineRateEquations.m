function [dfdt] = defineRateEquations(~,p,R,g,zerovec,gamd)
% t (double): an unused (but required) input parameter
% p (4x1 double): p(i) is the state population of state |i> with a particular velocity, n is the
%               number of basis states
% R (4x4 double): R(i,j) defines scattering rate for states |i> and |j> in units gam422
% g (4x4 double): g(i,j) defines decay rate from state |i> to |j> in units gam422
% zerovec (4x1 double): vector of zeros for each basis state
% gamd (double): decay rate (units: gam422) from 2P_{1/2} - 2D_{3/2} manifold
% dfdt (4x1 double): dfdt(k) is rate of change of population of state |k>, time in units gam422^-1

% This function defines rate equations for the 2S_{1/2}-2P_{1/2} LIF manifold, including decay from
% 2P_{1/2} to 2D_{3/2}, which is grafted on to the rate equations as a sync of the ion population.

dfdt = zerovec;
dfdt(1) = -R(3,1)*(p(1)-p(3))-R(4,1)*(p(1)-p(4))+g(3,1)*p(3)+g(4,1)*p(4);
dfdt(2) = -R(4,2)*(p(2)-p(4))-R(3,2)*(p(2)-p(3))+g(3,2)*p(3)+g(4,2)*p(4);
dfdt(3) = R(3,1)*(p(1)-p(3))+R(3,2)*(p(2)-p(3))-(g(3,1)+g(3,2))*p(3) - p(3)*gamd;
dfdt(4) = R(4,2)*(p(2)-p(4))+R(4,1)*(p(1)-p(4))-(g(4,2)+g(4,1))*p(4) - p(4)*gamd;

end
