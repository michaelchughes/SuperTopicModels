function [Data, Truth, Params] = genSynthData_EasyBarsMultiClass( dataParams )

defs.ALPH = 0.1;
defs.D = 500;
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

% Three classes:  
%  1 favors horizontal bars
eta(:,1) = [ 5*ones( 1, K/2)   -5*ones( 1, K/2) ];
%  2 favors vertical bars
eta(:,2) = [ -5*ones( 1, K/2)   5*ones( 1, K/2) ];
%  3 favors a mix
eta(:,3) = [ 2*ones( 1, K/2)   2*ones( 1, K/2) ];


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
    ybar = randn( size(mu) ) + mu;
    [~,maxid] = max( ybar );
    Data(dd).y = maxid;
    
    Ndk(dd,:) = curNdk;
    
    Truth.Topics{dd} = zs;
end
Truth.TopicObsPr = ObsPr;
Truth.eta = eta;



% % Plot example documents! 
% ys = vertcat( Data(:).y );
% histc( ys, 1:3 );
% close all;
% for cc = 1:size(eta,2)
%    exampleIDs = find(  ys == cc );
%    exampleIDs = shuffleVector( exampleIDs );
%    exampleIDs = exampleIDs( 1:min(30,length(exampleIDs) ) )';
%    L = length( exampleIDs );
%    
%       figure(cc);
%       set(gcf, 'Name', ['Class ' num2str(cc)] );
%       for rr = 1:L
%           ws = histc( Data( exampleIDs(rr) ).words, 1:V );
%           ws = ws/sum(ws);
%           docIm = reshape( ws, M, M);
%           subplot(3, 10, rr);
%           imagesc( docIm, [0 0.3] );
%           axis image;
%       end
%       colormap hot;
% end


% %DOUBLE CHECK: Can we recover the regression params??
% barZ = bsxfun( @rdivide, Ndk, sum(Ndk,2) );
% model = defaultModelParams_sLDA( Data );
% ys = getRealValuedRegressionOutcomes( Data, barZ, eta, model );
% 
% yTrue = vertcat( Data(:).y );
% yHat = zeros(size(yTrue));
% 
% for a = 1:25
%     
%     EE(a,:,:) = eta;
%     
%     for dd = 1:D
%         [~, yHat(dd)] = max( barZ(dd,:)*eta );
%     end
%     fprintf('y(1)=%.2f y(2)=%.2f | Lam %.2f | Error rate: %.3f\n', ys(1), ys(2), lambda, sum(yHat~=yTrue)/D );
%     ys = getRealValuedRegressionOutcomes( Data, barZ, eta, model );
% 
%     [eta,lambda] = sampleRegressionParams_sLDA( ys, barZ, eta, lambda, model.RegM);
% end
% 
% figure;
% plot( EE( :, :, 1) )
% 
% figure;
% plot( EE( :, :, 2) )