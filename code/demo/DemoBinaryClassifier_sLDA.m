% This demo shows how this toolbox can be used to do prediction
%  for binary classification given a set of documents (as bags of words).
% Specifically, we build a binary classification task on top of 
%  the toy bars data seen in "EasyDemo.m".
% Here, we'll use SUPERVISED LDA (aka sLDA)
%  which has a probit regression model built into it.
% Thus, when learning the latent topics for each document,
%  we actually use the binary label to help inform the topic distribution.
% See RunMCMCSimForSuperLDA.m for details of how this works.

close all;
clear variables;

jobID = 101;
taskID = 1;
dataP = {'EasyBarsBinary', 'D', 500, 'Nd', 50};
modelP = {'nTopics',10};
outP = {jobID, taskID};
algP = {'Niter', 250};
testP = {'D', 500, 'Nd', 50};
runSuperLDA( dataP, modelP, outP, algP, testP ); 


fprintf( 'Plotting predictions...\n');
plotPredictions( jobID, taskID );

fprintf( 'Remember, predictions are made using MCMC *samples* from the posterior.\n');
fprintf( '  We SHOULD expect a general upward trend toward a plateau, as the sampler improves on the initial configuration.\n' );
fprintf( '  However, we SHOULD NOT expect these samples to converge tightly... there may be lots of variance in the posterior.\n');


% To facilitate inspection...
INFO = loadSamplerInfo( jobID, taskID );
TrainData = INFO.Data;
TestData  = INFO.TestData;