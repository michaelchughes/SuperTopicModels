function [Data, Truth, Params] = genSynthData_EasyBarsBinaryLDA( dataParams )

defs.ALPH = 0.25;
defs.D = 1000;
defs.Nd = 100;
defs.K = 10;
defs.V = 25;
defs.lambda = 1;

Params = updateParamsWithUserInput( defs, dataParams );
% Add all fields of Params struct directly to the workspace
fnames = fieldnames( Params );
for ff = 1:length( fnames)
    eval( [fnames{ff} '= Params.' fnames{ff} ';'] );
end

eta = linspace( -10, 10, K )';

M = sqrt(V);

ObsPr = ones( K, V );
FAVOR_COUNT = 15;

% Horizontal bars
for rowID = 1:K/2
    favorIDs = ( rowID -1 )*M + (1:M);
    ObsPr( rowID, favorIDs) = FAVOR_COUNT;
end

for colID = 1:K/2
    favorIDs = (colID-1) + (1:M:V);
    ObsPr( K/2+colID, favorIDs) = FAVOR_COUNT;
end
ObsPr = bsxfun( @rdivide, ObsPr, sum(ObsPr,2) );

% ============================================  Generate data
for dd = 1:D
    pi = gamrnd( ALPH, 1, 1, K);
    pi = pi/sum(pi);
    
    zs = multinomial_many_draw( pi, Nd );
    
    curNdk = histc( zs, 1:K);
    ws = [];
    for kk = 1:K
        ws = [ws multinomial_many_draw( ObsPr(kk,:),  curNdk(kk) )];
    end
    
    assert( length(ws)==Nd );
    Data(dd).words = ws;
    
    mu = curNdk*eta  / Nd;
    ybar = randn(1,1) + mu;
    
    Data(dd).y = ybar > 0;
    
    Ndk(dd,:) = curNdk;
    
    Truth.Topics{dd} = zs;
end
Truth.TopicObsPr = ObsPr;
Truth.eta = eta;

% %DOUBLE CHECK: Can we recover the regression params??
% barZ = bsxfun( @rdivide, Ndk, sum(Ndk,2) );
% model = defaultModelParams_sLDA( Data );
% ys = getRealValuedRegressionOutcomes( Data, barZ, eta, model );
% 
% 
% yTrue = vertcat( Data(:).y );
% yHat = zeros(size(yTrue));
% 
% for a = 1:10
%         EE(:,a) = eta;
% 
%     for dd = 1:D
%         yHat(dd) = (barZ(dd,:)*eta ) > 0;
%     end
%     fprintf('y(1)=%.2f y(2)=%.2f | Lam %.2f | Error rate: %.3f\n', ys(1), ys(2), lambda, sum(yHat~=yTrue)/D );
%     ys = getRealValuedRegressionOutcomes( Data, barZ, eta, model );
% 
%     [eta,lambda] = sampleRegressionParams_sLDA( ys, barZ, eta, lambda, model.RegM);    
% end
% plot( EE' );

