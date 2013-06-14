% This demo shows how this toolbox can be used to do prediction
%  for regression problems (predict real-valued target variable y_d)
%  for each doc d in a set of documents (as bags of words).
% Specifically, we build a regression task on top of 
%  the toy bars data seen in "EasyDemo.m".
% Here, we'll use UNSUPERVISED LDA (aka vanilla Latent dirichlet alloc.)
%   and just fit regression model post-hoc onto it.
% So, this means the topics "z_d" for each document are learned 
%   *without* using the binary labels attached to the documents. 
% Then, once topic assignments are known, we learn a standard regression
%   model for mapping the topics to a real-valued score.
%   Observed regression score:
%     y_d ~ Normal( eta'*z_d , 1)
%           where z_d = [z_{d1} ... z_{dK}]  (a length-K vector)
%             and z_{dk} = frac of words assigned to topic k in doc d.
%   Basically, we just need to sample the global oefficient vector \eta 
%    from it's posterior given scores y_d and topic features z_d
% See RunMCMCSimForLDA.m for details of how each piece of the puzzle
%   unsupervised topic sampling and probit regression sampling fit
%   together.
close all;
clear variables;

jobID = 201;
taskID = 1;
dataP = {'EasyBarsRegress', 'D', 500, 'Nd', 50};
modelP = {'nTopics',10}; %cheating, giving the "correct" number of topics
outP = {jobID, taskID};
algP = {'Niter', 250};
testP = {'D', 500, 'Nd', 50};
runSuperLDA( dataP, modelP, outP, algP, testP ); 

fprintf( 'Plotting joint probability vs iter (on train data)...\n');
plotLogPr(jobID, taskID);

fprintf( 'Plotting predictions...\n');
plotPredictions( jobID, taskID );

fprintf( 'Remember, predictions are made using MCMC *samples* from the posterior.\n');
fprintf( '  We SHOULD expect a general upward trend toward a plateau, as the sampler improves on the initial configuration.\n' );
fprintf( '  However, we SHOULD NOT expect these samples to converge tightly... there may be lots of variance in the posterior.\n');


% To facilitate inspection...
INFO = loadSamplerInfo( jobID, taskID );
TrainData = INFO.Data;
TestData  = INFO.TestData;