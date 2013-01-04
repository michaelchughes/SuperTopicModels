function [Psi] = initLDA_FreshUniform( Data, model, initP )
% Initialize topics to completely random assignments,
%   no regard for the observed data

Psi.alpha = model.DocTopicM.alpha;
Psi.beta  = model.TopicObsM.beta;

D = length( Data );
K = model.nTopics;
V = max( horzcat( Data(:).words ) );

Psi.DTSuffStats.Ndk = zeros( D, K );
Psi.TWSuffStats.Nkw = zeros( K, V );

for dd = 1:D
   Psi.Topics{dd} = randi( [1 K],  length(Data(dd).words), 1 ); 
   Psi.DTSuffStats.Ndk(dd,:) = histc( Psi.Topics{dd}, 1:K );
   for nn = 1:length( Data(dd).words )
      Psi.TWSuffStats.Nkw( Psi.Topics{dd}(nn), Data(dd).words(nn) ) ...
        = Psi.TWSuffStats.Nkw( Psi.Topics{dd}(nn), Data(dd).words(nn) ) +1;
   end
end

Psi.TWSuffStats.Nk = sum( Psi.TWSuffStats.Nkw, 2 );
