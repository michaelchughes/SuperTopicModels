%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%         Truncated from Above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = randn_trunc_above( mu, sigma, b, N)
% Draws y ~ Norm( mu, sigma )  s.t.   -Inf <= y <= b
% Uses numerical methods based on inverting normal CDF
%   unless truncated region is so far (~5*sigma) away from the mean
%   that these methods become unstable 
% In these extreme cases, uses a rejection sampler with assymptotically 
%   perfect acceptance rates as distance from mean increases. 



if ( mu - b > 5*sigma )
    %  ------------------------------------------  Rejection Sampling
    y = mu - sigma*randn_trunc_tail_rejection( (mu-b)/sigma, N);
else
    c = mynormcdf( (b-mu)/sigma );
    u = c*rand(N,1);      % u ~ Unif(  0, c )
    y = mu + sigma*myinvnormcdf( u );
end

end