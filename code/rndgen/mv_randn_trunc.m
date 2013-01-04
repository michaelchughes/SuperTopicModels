function Y = mv_randn_trunc( mu, d, nSamples)
%  Draw vector y from truncated normal distribution
%    y ~  N( mu, eye(D) )  restricted s.t.  y(d) >= y(d') for all dims d'
% INPUTS
%   mu : Dx1 col vector indicating mean of the target Normal distribution
%   d  : integer index of dimension (1 <= d <= D) that should be biggest
%   nSamples  : number of samples to draw from this distribution
% OUTPUT
%   Y  :  DxN vector of indep. draws from truncated normal

if ~exist('nSamples','var')
    nSamples = 1;
end
nBurn = 3;

D = length(mu);
y = mu;

notids = [ 1:d-1  d+1:D ];

Y = zeros( D, nSamples );
rc = 0;
for rr = 1:(nBurn + nSamples)
   maxnotd = max( y(notids) ); 
   
   y(d) = randn_trunc_below( mu(d), 1, maxnotd, 1 );
   assert( y(d)> maxnotd, 'BAD!');

   for dd = shuffleVector( notids )
      y(dd) = randn_trunc_above( mu(dd), 1, y(d), 1 ); 
   end
   
   assert( y(d)> maxnotd, 'BAD!');
   
   if rr > nBurn
      rc = rc+1;
      Y(:, rc) = y; 
   end
end


end
