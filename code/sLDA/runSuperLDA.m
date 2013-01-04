function [ChainHist] = runSuperLDA( dataParams, modelParams, outParams, algParams, initParams, testP )
% runBPHMM
% User-facing entry function for configuring and executing MCMC simulation 
%   for posterior inference of a Beta Process HMM (BP-HMM) model.
% INPUT:
%  Takes five arguments, each a cell array that specifies parameters
%   as Name/Value pairs, overriding default values (see defaults/ dir)
%  dataParams : 
%      params for data preprocessing (# of sequences, block-averaging, etc.)
%  modelParams :
%      params that define posterior of interest, like prior hyperparameters
%  outputParams :
%      specifies how often to save, how often to display, etc.
%      saveEvery, logPrEvery, printEvery, statsEvery, etc.
%  algParams :
%      MCMC alg behavior (# of iterations, parameters of proposal distr.)
%  initParams :
%   initial Markov state (# of initial features, initial state seq., etc.)
%   {'InitFunc', @initBPHMMfromGroundTruth'} initializes to known stateSeq

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
testP  = sanitizeUserInput( testP );

% ================================================ INTERPRET INPUT
algDefs = defaultMCMCParams_LDA();
algParams = updateParamsWithUserInput(  algDefs, algParams );

outDefs = defaultOutputParams_LDA( outParams, algParams );
outParams = updateParamsWithUserInput( outDefs, outParams(3:end) );

initDefs = defaultInitParams_sLDA();
initParams = updateParamsWithUserInput( initDefs, initParams );

testDefs = defaultTestParams_sLDA();
testParams = updateParamsWithUserInput( testDefs, testP);

% ================================================= LOAD DATA
[Data, TestData, Truth] = loadBagOfWordsData( dataParams, testP);

model = defaultModelParams_sLDA( Data );
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
ChainHist = RunMCMCSimForSuperLDA(Data, Psi, algParams, outParams, model, TestData, testParams);