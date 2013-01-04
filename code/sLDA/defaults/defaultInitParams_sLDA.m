function initP = defaultInitParams_sLDA()
% Creates a struct encoding the default settings for MCMC inference

initP.InitFunc = @init_sLDA_FreshSeq;

