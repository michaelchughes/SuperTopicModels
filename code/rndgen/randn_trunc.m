function y = randn_trunc( mu, sigma, doA, nSamples)
%  Draw y from 1D truncated normal distribution
%    y ~  N( mu, sigma)   restricted s.t.  y <= 0  if doA = 0
%                                          y  > 0  if doA = 1
% INPUTS
%   mu : scalar mean of the target Normal distribution
%   sigma: scalar standard deviation of target Normal distribution
%   doA  : binary indicator for whether to truncate ABOVE zero (doA=1)
%                                             or  BELOW zero (doA=0)
%   nSamples  : number of samples to draw from this distribution
% OUTPUT
%   y  :  Nx1 vector of indep. draws from truncated normal

if nargin < 4
    nSamples = 1;
end

if ( doA == 1)
     y = randn_trunc_below( mu, sigma, 0, nSamples );    
else    
     y = randn_trunc_above( mu, sigma, 0, nSamples );    
end

end
