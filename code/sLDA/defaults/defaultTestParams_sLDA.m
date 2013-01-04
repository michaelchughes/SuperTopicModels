function testP = defaultTestParams_sLDA( )

% update predictions every X iters of sampler training
%   even during burnin!
testP.predictEvery = 25; 
testP.Niter = 16;
testP.BURNfrac = 0.25; % discard first fraction of the run