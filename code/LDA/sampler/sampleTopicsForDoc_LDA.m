function [zs, Ndk, Nkw, Nk] = sampleTopicsForDoc_LDA( ...
                     ws, zs, Ndk, Nkw, Nk,  ALPH, BETA, permIDs, myrands)
% Sample the topic assignments for a single document,
%  given observed words and fixed assignments for rest of corpus.
%SPEED
% Note: when possible, use the MEX version of this file.
% Same name (.cpp ext), same interface, and MUCH faster (>>10x).
%USAGE
%  sampleTopicsForDoc_LDA_DirMult(
%      ws, zs, Ndk, Nkw, Nk, ALPH, BETA,  seed* )
%INPUT
%  ws  : 1 x Nd vector of vocabword ids for each token in the document
%  zs  : 1 x Nd vector of topic assignments for each token in document
%  Ndk : 1 x K vector of counts of asgnments in current document
%  Nkw : K x V matrix of global topic-word asgn counts (for entire corpus)
%  Nk  : K x 1 vector of counts for each topic's prevalence corpus-wide
% ALPH : scalar symmetric Dir hyperparameter on Doc-Topic distribution
% BETA : scalar symmetric Dir hyperparameter on Topic-Word distribution

V = size( Nkw, 2);
Ndk = Ndk';
for n = 1:length(zs)
    nn = permIDs(n);
    
    zcur = zs(nn);
    
    if zcur > 0
        Ndk( zcur ) = Ndk( zcur )-1;
        Nkw( zcur, ws(nn) ) = Nkw( zcur, ws(nn) ) -1;
        Nk( zcur ) = Nk( zcur ) -1;
    end
    
    priorProbs = Ndk + ALPH;
    likProbs = ( Nkw( :, ws(nn) ) + BETA ) ./ ( Nk + V*BETA );
   
    znew = multinomial_single_draw_given_rand( priorProbs .* likProbs, myrands(nn) );
    
    zs(nn) = znew;
    
    Ndk( znew ) = Ndk( znew ) +1;
    Nkw( znew, ws(nn) ) = Nkw( znew, ws(nn) ) +1;
    Nk( znew ) = Nk( znew ) +1;
end