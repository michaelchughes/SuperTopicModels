function [Psi] = init_LDA_DebugMex( Data, model, initP)
% Easy script to verify that the MEXified sampler routine is functional
% Simply tries to initialize all documents to ground truth
%   and then run both the MEX sampler and native MATLAB sampler
% Verifies that their output topic samples are EXACTLY the same
%USAGE
%   runSuperLDA( {'EasyBarsLDA', 'D', 500, 'Nd', 100},{'nTopics',10},{1,1}, {'Niter', 0}, {'InitFunc', @init_LDA_DebugMex} );

Psi.alpha = model.DocTopicM.alpha;
Psi.beta  = model.TopicObsM.beta;

D = length( Data );
K = model.nTopics;
V = max( horzcat( Data(:).words ) );

INFO = loadSamplerInfo(1,1);
[Ndk,Nkw,Nk] = getSuffStatsFromTopics( INFO.Truth.Topics, Data );

for dd = 1:D
ws = Data(dd).words;
Nd = length(ws);

initzs = zeros(1, Nd );

permIDs = randperm( Nd );
randChoices = rand( 1, Nd );

zsMATLAB = sampleTopicsForDoc_LDA(...
    ws, initzs, ...
    Ndk(dd,:), Nkw, Nk,  ...
    Psi.alpha, Psi.beta, ...
    permIDs, randChoices);

fprintf('\n');

zsMEX = sampleTopicsForDoc_LDA_MEX(...
    ws, initzs, ...
    Ndk(dd,:), Nkw, Nk,  ...
    Psi.alpha, Psi.beta, ...
    permIDs, randChoices);

acc =  sum(zsMATLAB==zsMEX)/length(zsMEX);
fprintf('accuracy = %.3f\n', acc);

assert( acc > 0.9999, 'Badness');

end
Psi.DTSuffStats.Ndk = Ndk;
Psi.TWSuffStats.Nkw = Nkw;
Psi.TWSuffStats.Nk = Nk;
return;