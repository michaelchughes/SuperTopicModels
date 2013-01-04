function [Psi] = init_sLDA_FreshSeq( Data, model, initP)
% Initialize topics sequentially, one document at a time

Psi.alpha = model.DocTopicM.alpha;
Psi.beta  = model.TopicObsM.beta;

D = length( Data );
K = model.nTopics;
V = max( horzcat( Data(:).words ) );

if model.RegM.nClass > 0
    eta = zeros( K, model.RegM.nClass );
else
    eta = zeros( K, 1);
end
lambda = model.RegM.lambda.a /  model.RegM.lambda.b;

Ndk = zeros( D, K );
Nkw = zeros( K, V );
Nk = zeros( K, 1);

barZ = bsxfun(@rdivide, Ndk+eps, sum(Ndk+eps,2) );
ys = getRealValuedRegressionOutcomes( Data, barZ, eta, model );

for dd = randperm( D )
    ws = Data(dd).words;
    Nd = length(ws);
    initzs = zeros(1, Nd );
    permIDs = randperm( Nd );
    randChoices = rand( 1, Nd );
    [Psi.Topics{dd}, Ndk(dd,:),Nkw,Nk] = sampleTopicsForDoc_sLDA_MEX(...
                 ws, ys(dd,:), initzs, ...
                 eta, lambda, ...
                 Ndk(dd,:), Nkw, Nk,  ...
                 Psi.alpha, Psi.beta, ...
                 permIDs, randChoices, 0);   
             
end

Psi.DTSuffStats.Ndk = Ndk;
Psi.TWSuffStats.Nkw = Nkw;
Psi.TWSuffStats.Nk = Nk;
Psi.eta = eta;
Psi.lambda = lambda;
