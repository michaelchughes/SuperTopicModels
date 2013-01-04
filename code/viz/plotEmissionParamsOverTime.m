function [] = plotEmissionParamsOverTime( jobID, taskID )


X = loadSamplerOutput(jobID, taskID );

Psi = X.Psi(end);
Nkw = Psi.TWSuffStats.Nkw;

K = size( Nkw, 1);
Phi = Nkw + Psi.beta;
Phi = bsxfun(@rdivide, Phi, sum(Phi,2)  );
INFO = loadSamplerInfo(jobID, taskID );
Truth = INFO.Truth;
if ~isempty( Truth )
    Dist = zeros(K,K);
    for jj = 1:K
        for kk = 1:K
            Dist( jj, kk ) = sum( abs( Truth.TopicObsPr(jj,:) - Phi(kk,:) ) );
        end
    end
    est_labels = assignmentoptimal( Dist' );
end

hh = figure();
for iter = X.iters.Psi
    plotEmissionParams( jobID, taskID, iter, 'SamplerOutput', X, 'est_labels', est_labels, 'figHandle', hh );
    
    pause(.2);
end