function [Psi] = sLDAsample( Data, Psi, model, algP )

Nkw = Psi.TWSuffStats.Nkw;
Nk  = Psi.TWSuffStats.Nk;
Ndk = Psi.DTSuffStats.Ndk;
K = size( Ndk, 2);

barZ = bsxfun( @rdivide, Ndk, sum(Ndk,2) );
ys = getRealValuedRegressionOutcomes( Data, barZ, Psi.eta, model );

D = length(Data);
dcount = 0;
for dd = randperm( D )
    dcount = dcount + 1;
    
    Nd = length(Data(dd).words);
    permIDs = randperm(Nd );
    randChoices = rand(1, Nd);
    
    % ======================================  z | alpha, beta
    [Psi.Topics{dd}, Ndk(dd,:),Nkw, Nk] = sampleTopicsForDoc_sLDA_MEX(...
        Data(dd).words, ys(dd,:), Psi.Topics{dd}, ...
        Psi.eta, Psi.lambda, ...
        Ndk(dd,:), Nkw, Nk,  ...
        Psi.alpha, Psi.beta, ...
        permIDs,  randChoices, 0 );
    
    
    % ======================================  w,lam | y, z
    if mod( dcount, algP.sampleRegEvery )==0
        barZ = bsxfun( @rdivide, Ndk, sum(Ndk,2) );
        
        [Psi.eta, Psi.lambda] = sampleRegressionParams_sLDA( ys, barZ, Psi.eta, Psi.lambda, model.RegM );
        
        % Update the ys so they reflect most recent eta/lambda values
        ys = getRealValuedRegressionOutcomes( Data, barZ, Psi.eta, model );
    end
end
Psi.TWSuffStats.Nkw = Nkw;
Psi.TWSuffStats.Nk = Nk;
Psi.DTSuffStats.Ndk = Ndk;

end

