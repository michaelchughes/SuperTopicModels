function algP = defaultMCMCParams_LDA()
% Creates a struct encoding the default settings for MCMC inference

algP.Niter = 50;

algP.ObsM.doSampleHypers = 0;

algP.sampleRegEvery = 100; % every so many documents, resample reg coefs.