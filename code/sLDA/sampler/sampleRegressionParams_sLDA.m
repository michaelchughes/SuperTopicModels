function [eta, lambda] = sampleRegressionParams_sLDA( y, X, eta, lambda, RegModel)
% Sample regression parameters (w,lambda) for supervised LDA
%   For each document d:
%     y_d  ~  Normal(  w^T z_d, lambda*eye(K)  )

D = size(y,1);
K = size(X,2);

% ---------------------- Update precision lambda (fixed to one if classify)
if RegModel.doRegression
    SSD = sum( (y - X*eta).^2 );
    ap = D/2 + RegModel.lambda.a;
    bp = SSD/2 + RegModel.lambda.b;
    lambda = gamrnd( ap, 1./bp ); % Matlab parameterizes gammarnd as (a, 1/b)
else
    lambda = 1;
end

%  --------------------- Update regression coef vector Eta
%  Posterior for eta: Normal(  mean=mup, inv covar = invSp )
if RegModel.doRegression || RegModel.doBinary
    invSp = lambda*(X'*X)  + 1./RegModel.eta.var*eye(K);
    mup = invSp \ ( lambda*(X'*y) );
    eta = randmvnormal( mup, invSp );
elseif RegModel.doMultiClass  
    % remember: lambda must = 1 here
    invSp = X'*X  + 1./RegModel.eta.var*eye(K);
    for cc = 1:RegModel.nClass
        mu_cc = invSp \ ( X'*y(:,cc) );
        eta(:,cc) = randmvnormal( mu_cc, invSp );
    end    
end