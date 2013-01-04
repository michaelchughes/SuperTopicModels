function Psi = trainRegressionModel_LDA( Psi, Data, model, testP )
%  Given fixed topic definitions in training corpus,
%    train a regression model that predicts the supervised side info
%    using the empirical distribution of topic assignments

barZ = Psi.DTSuffStats.Ndk;
barZ = bsxfun( @rdivide, barZ, sum(barZ,2)  );

K = size(barZ,2);
if model.RegM.doMultiClass
    eta = zeros(K, model.RegM.nClass);
else
    eta = zeros(K,1);
end
lambda = 1;

for rr = 1:10
    ys = getRealValuedRegressionOutcomes( Data, barZ, eta, model);
    [eta, lambda] = sampleRegressionParams_sLDA( ys, barZ, eta, lambda, model.RegM );
end
Psi.eta = eta;
Psi.lambda = lambda;