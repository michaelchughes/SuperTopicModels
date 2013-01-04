function [Psi] = init_sLDA_DebugMex( Data, model, initP)
% Easy script to verify that the MEXified sampler routine is functional
% Simply tries to initialize all documents to ground truth
%   and then run both the MEX sampler and native MATLAB sampler
% Verifies that their output topic samples are EXACTLY the same
%USAGE
%   runSuperLDA( {'EasyBarsRegressLDA', 'D', 500, 'Nd', 100},{'nTopics',10},{1,1}, {'Niter', 0}, {'InitFunc', @init_sLDA_DebugMex} );

Psi.alpha = model.DocTopicM.alpha;
Psi.beta  = model.TopicObsM.beta;

D = length( Data );
K = model.nTopics;
V = max( horzcat( Data(:).words ) );

INFO = loadSamplerInfo(1,1);

if model.RegM.nClass > 0
    eta = INFO.Truth.eta;
else
    eta = INFO.Truth.eta;
    %eta = zeros( K, 1);
end
lambda = 1; 
%Ndk = zeros( D, K );
%Nkw = zeros( K, V );
%Nk = zeros( K, 1);

[Ndk,Nkw,Nk] = getSuffStatsFromTopics( INFO.Truth.Topics, Data );

barZ = bsxfun(@rdivide, Ndk, sum(Ndk,2) );
ys = getRealValuedRegressionOutcomes( Data, barZ, eta, model );

for dd = 1:length(Data)
ws = Data(dd).words;
Nd = length(ws);

initzs = zeros(1, Nd );

permIDs = randperm( Nd );
randChoices = rand( 1, Nd );

zsMATLAB = sampleTopicsForDoc_sLDA(...
    ws, ys(dd,:), initzs, ...
    eta, lambda, ...
    Ndk(dd,:), Nkw, Nk,  ...
    Psi.alpha, Psi.beta, ...
    permIDs, randChoices, 0);

zsMEX = sampleTopicsForDoc_sLDA_MEX(...
    ws, ys(dd,:), initzs, ...
    eta, lambda, ...
    Ndk(dd,:), Nkw, Nk,  ...
    Psi.alpha, Psi.beta, ...
    permIDs, randChoices, 0);

acc =  sum(zsMATLAB==zsMEX)/length(zsMEX);
fprintf('accuracy = %.3f\n', acc);

assert( acc > 0.9999, 'Badness');

end
Psi.DTSuffStats.Ndk = Ndk;
Psi.TWSuffStats.Nkw = Nkw;
Psi.TWSuffStats.Nk = Nk;
Psi.eta = eta;
Psi.lambda = lambda;
return;