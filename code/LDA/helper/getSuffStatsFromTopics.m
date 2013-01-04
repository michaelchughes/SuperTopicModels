function [Ndk, Nkw, Nk] =getSuffStatsFromTopics( Topics, Data )

D = length( Data );
zAll = [];
wAll = [];
for dd = 1:D
    zs = Topics{dd};
    zAll = [zAll zs];
    wAll = [wAll Data(dd).words];
end

K = max( zAll );
V = max( wAll );

Nkw = zeros( K, V );

for dd = 1:D
   zs = Topics{dd};
   Ndk(dd,:) = histc( zs, 1:K );
   for nn = 1:length(zs)
      Nkw(zs(nn),Data(dd).words(nn) ) =  Nkw(zs(nn),Data(dd).words(nn) )+1;  
   end
end

Nk = sum(Nkw,2);
Nk2 = sum(Ndk,1)';

assert( all(Nk==Nk2), 'Problem in calculation');