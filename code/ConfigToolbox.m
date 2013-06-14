addpath(genpath('.'));

% ---------------------------------  Compile MEX files
cd sLDA/sampler/
mex sampleTopicsForDoc_sLDA_MEX.cpp
cd ../..

cd LDA/sampler/
mex sampleTopicsForDoc_LDA_MEX.cpp
cd ../..

cd rndgen
mex randn_trunc_MEX.cpp
cd ..

EIGENPATH = getenv('EIGENPATH');
BOOSTPATH = getenv('BOOSTPATH');
if length( EIGENPATH) > 0 && length(BOOSTPATH) > 0 && exist(EIGENPATH,'dir') && exist(BOOSTPATH,'dir')
    CMD = sprintf( 'mex mv_randn_trunc_MEX.cpp -I%s -I%s', EIGENPATH, BOOSTPATH);
    cd rndgen
    eval(CMD);
    cd ..
    
else
    fprintf('WARNING: unable to build basic routines for multi-class sLDA sampler\n');
    fprintf('  please install Eigen and Boost libraries, link their paths to environment vars as described in README, and try again.\n');
end

fprintf('Done building MEX files.\n Remember to run all demos from current directory: PROJECTROOT/code/\n');