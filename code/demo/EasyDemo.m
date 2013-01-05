% Welcome!
% This easy demo first shows a toy bars dataset [Griffiths+Steyvers 04]
%  and then runs fast MCMC inference for LDA and visualizes the results!
% Make sure you've done these simple things to run this script:
%   -- Make code/ the wording directory
%       >> cd <SuperTopicsRoot>/code/
%   -- Create a local path preference file where results will be saved
%       This is just a one-line plain text file that contains
%         a valid absolute file path to where you'd like to store results
%       >> echo "path/to/my/result/directory/" > SimulationResults.path
% See QuickStartGuide.pdf in doc/ for details on configuring the toolbox

clear variables;
close all;

% This script automatically compiles necessary MEX functions
%     and adds all necessary paths to the MATLAB path
ConfigToolbox;

% -------------------------------------------------   CREATE TOY DATA!
fprintf( 'Creating some toy data...\n' );
% First, we'll create some toy data
%   D=500 documents, each with Nd=100 words
% Data : struct array where
%   Data(dd).words gives word list for the dd-th document 

dataParams = {'D', 500, 'Nd', 100, 'K',10};
[Data,Truth] = genSynthData_EasyBarsLDA( dataParams );


fprintf( 'Visualizing true topic-word distributions and example docs...\n');
% Visualize the 10 true topics
figure( 'Name', 'True Topic-Word Distributions');
set( gcf, 'Units', 'normalized', 'Position', [0.1 0.6 0.5 0.3] );
for kk = 1:10
    subplot( 2, 5, kk);
    sqIm = reshape( Truth.TopicObsPr(kk,:), 5,5);
    imagesc( sqIm, [0 0.5] );
    axis image;
end
colormap hot;

% Visualize some example documents
figure( 'Name', 'Example Documents');
set( gcf, 'Units', 'normalized', 'Position', [0.1 0.2 0.5 0.3] );
for dd = 1:10
    subplot( 2, 5, dd);
    wCounts = histc( Data(dd).words, 1:25 );
    sqIm = reshape( wCounts, 5,5);
    imagesc( sqIm, [0 25] );
    axis image;
end
colormap hot;
drawnow;
pause(1);

fprintf('Now running MCMC inference to recover topic-word distributions\n');
jobID = 1;
taskID = 1;
dataP = {'EasyBarsLDA'};
dataP(2:length(dataParams)+1) = dataParams(:);
modelP = {'nTopics',10};
outP = {jobID, taskID};
algP = {'Niter', 500};
runLDA( dataP, modelP, outP, algP ); 

fprintf( 'Plotting recovered \n' );
X = loadSamplerOutput( jobID, taskID );
Psi = X.Psi(end);
Nkw = Psi.TWSuffStats.Nkw;

EstTopicObsPr = bsxfun( @rdivide, Nkw, sum(Nkw,2) );
figure( 'Name', 'Recovered Topic-Word Distributions');
set( gcf, 'Units', 'normalized', 'Position', [0.1 0.6 0.5 0.3] );
for kk = 1:10
    subplot( 2, 5, kk);
    sqIm = reshape( EstTopicObsPr(kk,:), 5,5);
    imagesc( sqIm, [0 0.5] );
    axis image;
end
colormap hot;

fprintf( 'Remember: actual label ids for each topic are *irrelevant* from model perspective\n');
fprintf( '  what matters: *aligned* topics consistently assigned to same datapoints as ground truth\n' );


fprintf( 'Some tools in viz/ directory can aid visualization of sim results.\n');
fprintf( ' Here is a movie of the previous simulation.\n' );
plotEmissionParamsOverTime( 1, 1 );