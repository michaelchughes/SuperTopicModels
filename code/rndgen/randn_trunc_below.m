%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%         Truncated from Below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = randn_trunc_below( mu, sigma, a, N)
% Draws y ~ Norm( mu, sigma )  s.t.   a <= y <= +Inf
% Uses numerical methods based on inverting normal CDF
%   unless truncated region is so far (>5*sigma) away from the mean
%   that these methods become unstable 
% In these extreme cases, uses a rejection sampler with assymptotically 
%   perfect acceptance rates as distance from mean increases. 

if ( a - mu > 5*sigma)
    y = mu + sigma*randn_trunc_tail_rejection( (a-mu)/sigma, N );    
else
    c = mynormcdf( (a-mu)/sigma );
    u = c + (1-c)*rand(N,1);  % u ~ Uniform( c, 1) [eg c <= u <= 1]
    y = mu + sigma*myinvnormcdf( u );
end

end