function logPr = calcJointLogPr_LDA( Psi )
%  see equations 3 and 4 in Griffiths and Steyvers (2004)

Nkw = Psi.TWSuffStats.Nkw;
Nk  = Psi.TWSuffStats.Nk;
Ndk = Psi.DTSuffStats.Ndk; 

ALPH = Psi.alpha;
BETA = Psi.beta;
V = size( Nkw, 2 );

%  Log Prob for w|z,BETA --------------------------------------------
logPr.w = sum( sum( gammaln( Nkw + BETA ) ) ) - sum( gammaln( Nk + V*BETA ) );

%  Log Prob for z|Alpha --------------------------------------------
logPr.z = sum( sum( gammaln( Ndk + ALPH ) ) );

% Compute joint probability of *all* variables
fNames = fieldnames(logPr);
logPr.all = 0;
for n = 1:length( fNames )
    logPr.all = logPr.all + logPr.(fNames{n});
end