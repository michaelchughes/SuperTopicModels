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

dataName = dataParams{1};

if dataName(1) ~= filesep
    DATA_DIR = 'data/';
    matFilePath = fullfile(DATA_DIR, dataName);
else
    matFilePath = dataName;
end

dataN = lower(dataName);
if ~isempty( strfind(dataN, 'bars' ))
    genFunction = eval( ['@genSynthData_' dataName] );
    [Data, Truth] = genFunction( dataParams(2:end) );
    
    if exist( 'testParams', 'var')
        [TestData, TestTruth] = genFunction( testParams );
    end
    
elseif ~isempty( matFilePath )
    Data = load(  matFilePath );
    Data = Data.Data;
    
    if isnumeric( Data )
       Dstruct = struct();
       D = size( Data,1);
       V = size( Data,2);
       for dd = 1:D
          ws = [];
          for vv = 1:V
            if Data(dd,vv) > 0  
              ws = [ws repmat(vv, 1, Data(dd,vv))];
            end            
          end
          Dstruct(dd).words = ws;
       end
       
       Data = Dstruct;
        
    end
    
else
    error(['Dataset ' dataName ' not found.'] );
end


% ---------------------------------------------------------  Reset stream
curStream = RandStream.getDefaultStream();
curStream.State = entryState;

end