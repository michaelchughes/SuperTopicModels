%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%         Rejection sampler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = randn_trunc_tail_rejection(a, N)
%  Draw sample from tail of standard normal via rejection sampling
%   based on a Rayleigh proposal distribution
%  This scheme has acceptance probability approach 1
%   as a -> +Inf (e.g. as truncation boundary gets further from mean )

if nargin < 2
    N = 1;
end

x = sqrt( a^2 - 2 *log(rand(N,1) ) );

while (  all(x.*rand(N,1) > a)  )
    x = sqrt( a^2 - 2 *log(rand(N,1) ) );
end

end