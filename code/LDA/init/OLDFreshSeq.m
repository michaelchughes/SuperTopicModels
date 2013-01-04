function [Psi] = initLDA_FreshSeq( Data, model, initP)
% Initialize topics sequentially, one document at a time

Psi.alpha = model.DocTopicM.alpha;
Psi.beta  = model.TopicObsM.beta;

D = length( Data );
K = model.nTopics;
V = max( horzcat( Data(:).words ) );

eta = zeros( K, 1);
lambda = model.RegM.lambda.a /  model.RegM.lambda.b; 

Ndk = zeros( D, K );
Nkw = zeros( K, V );
Nk = zeros( K, 1);
for dd = randperm( length( Data ) )
    ws = Data(dd).words;
    initzs = zeros( size(ws) );
       [Psi.Topics{dd}, Ndk(dd,:),Nkw,Nk] = sampleTopicsForDoc_sLDA_DirMult(...
                 ws, Data(dd).y, initzs, ...
                 eta, lambda, ...
                 Ndk(dd,:), Nkw, Nk,  ...
                 Psi.alpha, Psi.beta );   
             
end

Psi.DTSuffStats.Ndk = Ndk;
Psi.TWSuffStats.Nkw = Nkw;
Psi.TWSuffStats.Nk = Nk;
Psi.eta = eta;
Psi.lambda = lambda;
