function [] = plotPredictions( jobIDs, taskIDs, queryDataName )

if ~exist( 'varName', 'var' ) || isempty( varName )
    varName = 'all';
end

figure;
set( gcf, 'Units', 'normalized', 'Position', [0.5 0.5 0.5 0.5] );
if exist( 'jobNames', 'var' ) && ~isempty( jobNames )
    plotColors  = get(0,'defaultAxesColorOrder'); %jet( length( jobNames )  );
    hold on;
else
    plotColors = get(0,'defaultAxesColorOrder');
    hold all;
end


for jobID = jobIDs
    taskNames = {};
    
    doFirstTask = 1;
    for taskID = taskIDs
        
        DATA_DIR = getUserSpecifiedPath( 'SimulationResults');
        DATA_DIR = fullfile( DATA_DIR, num2str(jobID), num2str(taskID) );
        P = load( fullfile( DATA_DIR, 'Predictions.mat') );
        if ~isstruct(P)
            continue;
        end
        
        if exist( 'jobNames', 'var' )  && ~isempty( jobNames )
            jj = 1 + mod( find( jobID == jobIDs )-1, size(plotColors,1) );
            curColor = plotColors(jj,:);
            if doFirstTask;
                taskVis = 'on';
                doFirstTask = 0;
            else
                taskVis = 'off';
            end
        else
            if taskID == taskIDs(1)
                taskNames{end+1} = 'train';
                taskNames{end+1} = 'test';
            end
            
            jj = 1 + mod( find( taskID == taskIDs )-1, size(plotColors,1) );
            curColor = plotColors(jj,:);
            taskVis = 'on';
        end
        
        if isfield( P.train, 'acc' )
            plot( P.iters, [P.train(:).acc], '+--', 'MarkerSize', 10, 'LineWidth', 1.5, ...
                'HandleVisibility', taskVis, 'Color', curColor );
            
            plot( P.iters, [P.test(:).acc], '.-', 'MarkerSize', 20, 'LineWidth', 1.5, ...
                'HandleVisibility', taskVis, 'Color', curColor );
            metricName = 'accuracy';
        elseif isfield( P.train, 'pR2' )
            
            plot( P.iters, [P.train(:).pR2], '+--', 'MarkerSize', 10, 'LineWidth', 1.5, ...
                'HandleVisibility', taskVis, 'Color', curColor );
            
            plot( P.iters, [P.test(:).pR2], '.-', 'MarkerSize', 20, 'LineWidth', 1.5, ...
                'HandleVisibility', taskVis, 'Color', curColor );
            metricName = 'predictive R2';
        end
    end
end

if exist( 'jobNames', 'var' ) && ~isempty( jobNames )
    legend( jobNames , 'Location', 'SouthEast' );
else
    legend( taskNames, 'Location', 'SouthEast' );
end

if strcmp( varName, 'all' )
    varName = '';
end

ylabel( metricName, 'FontSize', 18 );

xlabel ('iteration', 'FontSize', 18);

grid on;
set( gca, 'FontSize', 16 );
end % main function