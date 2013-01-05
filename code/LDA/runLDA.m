function [ChainHist] = runLDA( dataParams, modelParams, outParams, algParams, initParams, testP )
% User-facing entry function for configuring and executing MCMC simulation 
%   for posterior inference of standard Latent Dirichlet Allocation (LDA)
% Specifically intended for tasks involving *predictive* side information
%   Given, a document corpus where each doc has a held out label,
%    we first train UNSUPERVISED LDA model using just the bag-of-words info
%    then post-hoc fit a regression/classification model 
%      that discovers how to map between hidden topics and observed labels
%USAGE:
%  To run MCMC using mostly defaults on the classic "toy bars" dataset,
%  >> runLDA( {'EasyBarsLDA'},{'nTopics',10},{1,1}, {'Niter', 100} );
%  See DemoLDAEasy.m in demos/ for more details.
%INPUT:
%  Takes 6 arguments, each a cell array that specifies parameters
%   as Name/Value pairs, overriding default values (see defaults/ dir)
%  dataParams : 
%      params for data preprocessing (# of documents, etc.)
%  modelParams :
%      params that define posterior of interest, like:
%        the number of topics in the model, or prior hyperparameter values%           
%  outputParams :
%      specifies how often to save, how often to display.
%      key params: saveEvery, logPrEvery, printEvery, statsEvery, etc.
%  algParams :
%      MCMC alg behavior (# of iterations, parameters of proposal distr.)
%  initParams :
%   how to construct initial Markov state
%   default (recommended) setting: {'InitFunc', @initLDA_FreshSeq'}
%     *sequentially* visits each doc and assign topics given prev. asgns.
%  testParams :
addpath( genpath('.'));

if ~exist( 'algParams','var')  algParams = {}; end;
if ~exist( 'initParams','var') initParams = {}; end;
if ~exist( 'testP', 'var') testP = {}; end;

% ================================================ SANITIZE INPUT
% Converts strings to doubles when possible, etc. to allow command line input
dataParams  = sanitizeUserInput( dataParams );
modelParams = sanitizeUserInput( modelParams );
outParams   = sanitizeUserInput( outParams );
algParams   = sanitizeUserInput( algParams );
initParams  = sanitizeUserInput( initParams );

% ================================================ INTERPRET INPUT
algDefs = defaultMCMCParams_LDA();
algParams = updateParamsWithUserInput(  algDefs, algParams );

outDefs = defaultOutputParams_LDA( outParams, algParams );
outParams = updateParamsWithUserInput( outDefs, outParams(3:end) );

initDefs = defaultInitParams_LDA();
initParams = updateParamsWithUserInput( initDefs, initParams );

testDefs = defaultTestParams_sLDA();
testParams = updateParamsWithUserInput( testDefs, testP);

% ================================================= LOAD DATA
[Data, TestData, Truth] = loadBagOfWordsData( dataParams, testP );

model = defaultModelParams_LDA( Data );
model = updateParamsWithUserInput( model, modelParams );

info_fname = fullfile( outParams.saveDir, 'Info.mat');
save( info_fname, 'Data', 'TestData', 'Truth', 'model', 'initParams', 'algParams', 'outParams' );

% ================================================= SET MCMC RAND NUM SEED
jobStr = num2str( outParams.jobID );
taskStr = num2str( outParams.taskID );
SEED = force2double( [jobStr taskStr] );

% MATLAB's seed (takes one integer)
SEED = mod( SEED, 2^32);
RandStream.setGlobalStream( RandStream( 'twister', 'Seed', SEED )   );

% ================================================= INITIALIZE MODEL
% Note: InitFunc will often use own random seed (reset internally only)
%   so that different sampling algs can be compared on *same* init state
[Psi] = initParams.InitFunc( Data, model, initParams );

% ================================================= RUN INFERENCE
ChainHist = RunMCMCSimForLDA(Data, Psi, algParams, outParams, model, TestData, testParams);