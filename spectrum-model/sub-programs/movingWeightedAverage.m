function [sig_filt] = movingWeightedAverage(sig_in,weights,numpts,pow)
% sig_in (1xn double): signal to be filtered
% sig_se (1xn double): standard error corresponding to sig_in, should have same units as sig_in
% numpts (1x1 double): 2*numpts+1 is the number of points to use in moving average
% pow (double): power to raise weights to
% sig_filt (1xn double): filtered signal

% This function filters 'sig_in' using a local weighted average where the weights are given by the
% inverse standard error.

if nargin < 4
    pow = 1;
end

sig_filt = zeros(size(sig_in));
weights = abs(weights).^pow;
for i = 1:length(sig_in)
    ind = i-numpts:i+numpts; % identify indices of local pts to use in weighted average
    ind = ind(ind > 0 & ind < length(sig_in)); % make sure indices are within range of sig_in (matters for edge cases)
    
    temp = sig_in(ind); % points to average
    temp_w = weights(ind)./sum(weights(ind)); % weights for average renormalized so sum(temp_w) = 1
    
    sig_filt(i) = sum(temp.*temp_w); % compute local weighted average
end

end