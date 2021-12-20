function [w] = getSimpsonWeights(x)
% x (1xn double): vector containing independent variable
% w (nx1 double): weights for Simpson integration

% make sure 'x' has linearly-spaced values
dx = diff(x); % get difference between adjacent elements
if sum(diff(dx)) > 1e-5
    error('Simpson weights only defined for vectors with uniformly-spaced elements.')
end

% get Simpson Weights
w = ones(size(x')).*2;
w(1) = 1;
w(3:2:end) = 4;
w(end) = 1;
w = w.*dx(1)/3;

end