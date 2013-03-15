% This demo shows how this toolbox can be used to do prediction
%  for multi-way classification given a set of documents (as bags of words).
% We'll use a simple 3-class 

close all;
clear variables;

fprintf( 'Generating toy data with C=3 possible class labels\n' );
[Data,Truth] = genSynthData_EasyBarsMultiClass( {'D', 500, 'Nd', 100} );

fprintf( 'Each doc in the Data array has "y" field that gives the integer label (1...C) of the class it belongs to\n');
fprintf( 'Data(1).y = %d\n', Data(1).y );
fprintf( 'Data(15).y = %d, and\n', Data(15).y );
fprintf( 'Data(77).y = %d, etc.\n', Data(77).y );

fprintf( 'Visualizing example documents...\n' );
fprintf( 'Class 1 favors vertical bars\n' );
fprintf( 'Class 2 favors horizontal bars\n' );
fprintf( 'Class 3 favors a mixture\n' );
ys = vertcat( Data(:).y );
for cc = 1:3
    figure( 'Name', sprintf('Class %d: Example Docs',cc)  );
    set( gcf, 'Units', 'normalized', 'Position', [0.3*(cc-1)  0.5 0.3 0.3] );
    exampleIDs = find(  ys == cc );
    exampleIDs = exampleIDs( 1:min(25,length(exampleIDs) ) )';
    L = length( exampleIDs );
    
    for rr = 1:L
        subplot(5, 5, rr);
        wCounts = histc( Data( exampleIDs(rr) ).words, 1:25 );
        wCounts = wCounts/sum(wCounts);
        sqIm = reshape( wCounts, 5,5);
        imagesc( sqIm, [0 .25] );
        axis image;
    end
    colormap hot;
    drawnow;
end


fprintf('Running MCMC on this dataset...\n' );
fprintf('... this training process is quite a bit slower (for now)\n' );
fprintf('    due to complications involved in multiclass probit regression\n');
jobID = 301;
taskID = 1;
dataP = {'EasyBarsMultiClass', 'D', 500, 'Nd', 100};
modelP = {'nTopics',10};
outP = {jobID, taskID};
algP = {'Niter', 250};
testP = {'D', 100, 'Nd', 100};
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