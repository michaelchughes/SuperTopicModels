function logPr = calcJointLogPr_sLDA( Psi, Data, model )

Nkw = Psi.TWSuffStats.Nkw;
Nk  = Psi.TWSuffStats.Nk;
Ndk = Psi.DTSuffStats.Ndk; 

ALPH = Psi.alpha;
BETA = Psi.beta;
K = size( Nkw, 1 );
V = size( Nkw, 2 );
Nd = sum( Ndk, 1);

%  Log Prob for w|z,BETA --------------------------------------------
logPr.w = sum( sum( gammaln( Nkw + BETA ) ) ) - sum( gammaln( Nk + V*BETA ) );

%  Log Prob for z|Alpha --------------------------------------------
logPr.z = sum( sum( gammaln( Ndk + ALPH ) ) ) - sum( gammaln( Nd + K*ALPH));

%  Log Prob for y | z, eta, lambda
ys = vertcat( Data(:).y );
barZ = bsxfun( @rdivide, Ndk, sum(Ndk,2) );
lambda = Psi.lambda;
eta = Psi.eta;
if model.RegM.doRegression
    logPr.y = sum( 0.5*log( lambda ) - 0.5*lambda*( ys - barZ*eta ).^2 );
elseif model.RegM.doBinary
    %  y ~ Bernoulli( Probit(mu) )
    mu = barZ*eta;
    mu = min( mu, 7 );  %Enforce bounds for numerical stability
    mu = max( mu, -7);  % mu outside these bounds ~= 0 or 1
    pON = mynormcdf( mu );
    logPr.y = sum( ys.*log(pON) + (1-ys).*log(1-pON) );

elseif model.RegM.doMultiClass
    %  Not sure how to do this analytically,
    %     so we can just draw auxiliary parameters and get their probs
    yAux = getRealValuedRegressionOutcomes( Data, barZ, eta, model);
    logPr.y = sum( sum( - 0.5*lambda*( yAux - barZ*eta ).^2 ) );
end 

%  Log Prob for eta, lambda | Regression Model prior params
RegM = model.RegM;
if model.RegM.doMultiClass
    logPr.eta = sum( sum( -0.5/RegM.eta.var * (eta).^2) );
else
    logPr.eta = sum( -0.5/RegM.eta.var * (eta).^2 );
end
if model.RegM.doRegression
    logPr.lambda = (RegM.lambda.a-1)*log( lambda ) - RegM.lambda.b*lambda;
end

% Compute joint probability of *all* variables
fNames = fieldnames(logPr);
logPr.all = 0;
for n = 1:length( fNames )
    logPr.all = logPr.all + logPr.(fNames{n});
end

assert( 1==1);