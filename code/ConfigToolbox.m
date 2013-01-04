addpath(genpath('.'));

% ---------------------------------  Compile MEX files
cd sLDA/sampler/
mex sampleTopicsForDoc_sLDA_MEX.cpp
cd ../..

cd LDA/sampler/
mex sampleTopicsForDoc_LDA_MEX.cpp
cd ../..

