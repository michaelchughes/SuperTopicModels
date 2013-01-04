function [Psi] = LDAsample( Data, Psi, model, algP )

% ======================================  z(dd) | alpha, beta
D = length( Data );
Nkw = Psi.TWSuffStats.Nkw;
Nk  = Psi.TWSuffStats.Nk;
Ndk = Psi.DTSuffStats.Ndk;
seeds = randi( [1 100000], D, 1);
for dd = randperm(D)
    
    Nd = length(Data(dd).words);
    permIDs = randperm(Nd );
    randChoices = rand(1, Nd);
    [Psi.Topics{dd}, Ndk(dd,:),Nkw,Nk] = sampleTopicsForDoc_LDA_MEX(...
        Data(dd).words, Psi.Topics{dd}, ...
        Ndk(dd,:), Nkw, Nk,  ...
        Psi.alpha, Psi.beta, permIDs, randChoices );
end
Psi.TWSuffStats.Nkw = Nkw;
Psi.TWSuffStats.Nk = Nk;
Psi.DTSuffStats.Ndk = Ndk;

% ====================================== TO DO: sampler hypers

end