function [Result] = predict_LDA( Psi, Data, model, testP )
% Given a trained LDA model Psi
%   which includes a post-hoc trained regression model
%         ala trainRegressionModel_LDA.m
%   fit its topics to each doc in a held out dataset "Data"
%    and then predict each doc's response (regression or classification)

D=length( Data );
yTrue = zeros(D,1);
yHat  = zeros(D,1);
pON   = zeros(D,1);

for dd = 1:D
    Nkw = Psi.TWSuffStats.Nkw;
    Nk  = sum( Nkw,2);
    
    K = size( Nkw,1);
    
    curNdk = zeros( K,1 );
    zs = zeros( size(Data(dd).words)  );
    
    rc = 0;
    
    % Make sure we have odd # of samples so there are no ties in prediction!
    BURNiter = ceil( testP.BURNfrac * testP.Niter );
    if mod(testP.Niter - BURNiter,2)==0
        testP.Niter = testP.Niter+1;
    end
    
    if model.RegM.doMultiClass
        yHats = zeros( testP.Niter - BURNiter, model.RegM.nClass );
    else
        yHats = zeros( testP.Niter - BURNiter, 1);
    end
    for rr = 1:testP.Niter
        [zs, curNdk, Nkw, Nk] = sampleTopicsForDoc_LDA_DirMult(...
            Data(dd).words, zs, ...
            curNdk', Nkw, Nk,  ...
            Psi.alpha, Psi.beta );
        
        if rr > BURNiter
            rc = rc+1;
            yHats( rc, : ) = curNdk'*Psi.eta / length(zs);
        end
    end
    
    yTrue(dd) = Data(dd).y;
    if model.RegM.doRegression
        yHat(dd) = mean( yHats );
    elseif model.RegM.doBinary
        yBin = double( yHats > 0 );
        
        yHat(dd) = mode( yBin );
        pON(dd)  = mean( yBin );
    elseif model.RegM.doMultiClass
        [~,yLabel] = max( yHats, [], 2 );
        yHat(dd) = mode( yLabel );
    end
    
end


Result.yHat = yHat;
Result.yTrue = yTrue;

if model.RegM.doRegression
    Result.pR2 = predictiveR2( yTrue, yHat );
elseif model.RegM.doBinary || model.RegM.doMultiClass
    Result.acc = sum( yHat==yTrue )/length( yTrue );
end



end