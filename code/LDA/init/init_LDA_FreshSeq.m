function [Psi] = init_LDA_FreshSeq( Data, model, initP)
% Initialize topics sequentially, one document at a time

Psi.alpha = model.DocTopicM.alpha;
Psi.beta  = model.TopicObsM.beta;

D = length( Data );
K = model.nTopics;
V = max( horzcat( Data(:).words ) );

Ndk = zeros( D, K );
Nkw = zeros( K, V );
Nk = zeros( K, 1);
for dd = randperm( length( Data ) )
    ws = Data(dd).words;
    initzs = zeros( size(ws) );
    permIDs = randperm(length(ws));
    randChoices = rand(1, length(ws) );
    [zs, Ndk(dd,:),Nkw,Nk] = sampleTopicsForDoc_LDA_MEX(...
                 ws, initzs, ...
                 Ndk(dd,:), Nkw, Nk,  ...
                 Psi.alpha, Psi.beta, ...
                 permIDs, randChoices);  
    Psi.Topics{dd} = zs;
end

Psi.DTSuffStats.Ndk = Ndk;
Psi.TWSuffStats.Nkw = Nkw;
Psi.TWSuffStats.Nk = Nk;
