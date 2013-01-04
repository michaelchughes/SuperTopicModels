function [zs, Ndk, Nkw, Nk] = sampleTopicsForDoc_sLDA_DirMult( ws, y, zs, ...
                                 eta, lambda, Ndk, Nkw, Nk, ALPH, BETA, ...
                                 permIDs, randChoices, doDebug)
%USAGE
%INPUT
%   ws :
%    y : real-valued target for regression in current document
%              observed in regression probs
%              auxiliary var for classification (via probit Gibbs sampler)
%        1 x 1 in regression/binary classification
%        1 x C in multiclass classification (largest index = class id)
%   zs :
%   eta: regression coeficient weight vector
%        K x 1 in regr/binary classification 
%        K x C

K = size( Nkw, 1);
V = size( Nkw, 2);
Ndk = Ndk';

Nd = length(zs);

if isvector(eta)
    mu = eta'*Ndk;
else
    mu = Ndk'*eta;  % enforce mu = 1xC
end

for n = 1:length(permIDs)
    nn = permIDs(n);
   
    zcur = zs(nn);
    
    if zcur > 0
        Ndk( zcur ) = Ndk( zcur )-1;
        Nkw( zcur, ws(nn) ) = Nkw( zcur, ws(nn) ) -1;
        Nk( zcur ) = Nk( zcur ) -1;
        
        if isvector(eta)
          mu = mu - eta( zcur );
        else
          mu = mu - eta( zcur, :);
        end
    end
    
    priorProbs = Ndk + ALPH;
    likProbs = ( Nkw( :, ws(nn) ) + BETA ) ./ ( Nk + V*BETA );
    
    
    if isvector(eta)
        % regression/binary:
        %   muHat(kk) = eta'*Ndk  if  z(nn) is set to topic kk
        %             = eta'*Ndk (excluding nn)    + eta(kk)
        %             =     mu                     + eta(kk)
        yDiff  = y - (mu+eta)/Nd;
        yProbs = exp( -0.5*lambda*yDiff.^2 );
    else
        % multiclass
        %  muHat(kk,cc) = eta(:,cc)'*Ndk ... if z(nn) is set to topic kk
        muHat = bsxfun(@plus, eta, mu );
        diff = sum( bsxfun(@minus, muHat/Nd, y), 2);
        yProbs   =  exp( -0.5*diff.*diff  );
        
%         yCheck = zeros(K,1);
%         for kk = 1:K
%             barZd = Ndk;
%             barZd(kk) = barZd(kk)+1;
%             mud = barZd'*eta;  % 1 x C
%             SSD = sum( (mud/Nd - y).^2);
%             yCheck(kk) = exp( -0.5*SSD );
%         end
%         assert(  allEq(yCheck, yProbs)  );
    end
       
    
    znew = multinomial_single_draw_given_rand( ...
                    priorProbs .* likProbs.*yProbs, randChoices(nn) );
    zs(nn) = znew;
    
    
    if doDebug && ( n <= 3 || n >= length(zs)-2 )
       ps = priorProbs .* likProbs .* yProbs;
       ps = ps/sum(ps);
       %fprintf( '%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f \n', yProbs);
       fprintf( ' %2d %.3f | %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f --> %d\n', ...
                nn-1, randChoices(nn), ps(1), ps(2), ps(3), ps(4),ps(5),ps(6),ps(7),ps(8),ps(9), ps(10), znew );
    end
    
    Ndk( znew ) = Ndk( znew ) +1;
    Nkw( znew, ws(nn) ) = Nkw( znew, ws(nn) ) +1;
    Nk( znew ) = Nk( znew ) +1;
    if isvector( eta )
        mu = mu + eta( znew );
    else
        mu = mu + eta(znew,:);
    end
    
    %assert( all( abs(mu - Ndk'*eta) < 1e-8 )  );
end