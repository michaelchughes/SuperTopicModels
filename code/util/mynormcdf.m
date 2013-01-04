%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   CDF of the Standard Normal Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = mynormcdf(  x )
%  Find probability p s.t.
%     Pr( X <= x ) = p,  where random variable X ~ Normal( 0, 1 )
% Notes:
% - Assume x is already standardized ( e.g. x = (x-Mu)/Sigma  )
% - Use the complementary error function instead of .5*(1+erf(z/sqrt(2))),
%      to produce accurate near-zero results for large negative x.

p = 0.5 * erfc(-x ./ sqrt(2));
end
