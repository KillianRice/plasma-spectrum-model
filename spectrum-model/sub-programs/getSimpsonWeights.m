function [w] = getSimpsonWeights(x)
% x (1xn double): vector containing independent variable
% w (nx1 double): weights for Simpson integration

%% Obtain Simpson Weights for Integration
% Make sure 'x' has linearly-spaced values
dx = unique(round(diff(x),5,'significant'));
if length(dx) ~= 1
    error('Simpson weights only defined for vectors with uniformly-spaced elements.')
end

% Get Simpson Weights
w = ones(size(x'));
w(2:2:end) = 2;
w(3:2:end) = 4;
w(end) = 1;
w = w.*dx/3;

end