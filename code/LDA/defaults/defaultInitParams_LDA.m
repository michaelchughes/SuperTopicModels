function initP = defaultInitParams_LDA()
% Creates a struct encoding the default settings for MCMC inference

initP.InitFunc = @init_LDA_FreshSeq;

