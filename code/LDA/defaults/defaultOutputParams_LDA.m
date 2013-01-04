function outP = defaultOutputParams_LDA( outParams, algP)
% Creates a struct encoding the default outP for MCMC inference


saveDir = getUserSpecifiedPath( 'SimulationResults' );

for aa = 1:length( outParams )
    switch aa
        case 1
            jobID = force2double(  outParams{aa} );
        case 2
            taskID = force2double( outParams{aa} );
    end
end
outP.jobID = jobID;
outP.taskID = taskID;
outP.saveDir = fullfile( saveDir, num2str(jobID), num2str(taskID) );
if ~exist( outP.saveDir, 'dir' )
    [~,~] = mkdir( outP.saveDir );
end

if isfield( algP, 'TimeLimit' ) && ~isempty( algP.TimeLimit )
   TL = algP.TimeLimit;
   Niter = Inf;
else
   TL = Inf;
   Niter = algP.Niter;
end

if TL <= 5*60 || Niter <= 200
    outP.saveEvery = 5;
    outP.printEvery = 5;
    outP.logPrEvery = 1;
    outP.statsEvery = 1;
elseif TL <= 2*3600 || Niter <= 5000    
    outP.saveEvery = 25;
    outP.printEvery = 25;
    outP.logPrEvery = 5;    
    outP.statsEvery = 5;
else
    outP.saveEvery = 50;
    outP.printEvery = 100;
    outP.logPrEvery = 10;    
    outP.statsEvery = 10;
end


