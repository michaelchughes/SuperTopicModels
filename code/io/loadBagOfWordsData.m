function [Data, TestData, Truth] = loadBagOfWordsData( dataParams, testParams )
% Helper for loading datasets for unsupervised/supervised analysis
% NB: For synthetic datasets, we *always* use the same PRNG seed
%   when generating data, so things are consistent across runs

% ------------------------------- Remember old state to use again afterward
curStream = RandStream.getGlobalStream();
entryState = curStream.State;

% Reset PRNG state to default value with SEED 0
%       so that we always get same synth data regardless of when called
reset( RandStream.getGlobalStream(), 0);


Truth = [];
TestData = [];
DATA_DIR = 'data/';
dataName = dataParams{1};

matFileList = dir( fullfile(DATA_DIR, dataName, '*.mat') );

dataN = lower(dataName);
if ~isempty( strfind( dataN, 'toy' ) ) ||  ~isempty( strfind(dataN, 'bars' ))
    genFunction = eval( ['@genSynthData_' dataName] );
    [Data, Truth] = genFunction( dataParams(2:end) );
    
    if exist( 'testParams', 'var')
        [TestData, TestTruth] = genFunction( testParams );
    end
    
elseif ~isempty( matFileList )
    Data = load(     fullfile(DATA_DIR, dataName, matFileList(1).name ) );
    Data = Data.Data;
else
    error(['Dataset ' dataName ' not found.'] );
end


% ---------------------------------------------------------  Reset stream
curStream = RandStream.getDefaultStream();
curStream.State = entryState;