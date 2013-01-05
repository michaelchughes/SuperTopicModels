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
set( gcf, 'Units', 'normalized', 'Position', [0.1 0.6 0.5 0.3] );
for iter = X.iters.Psi
    plotEmissionParams( jobID, taskID, iter, 'SamplerOutput', X, 'est_labels', est_labels, 'figHandle', hh );
    
    annotation('textbox', [0, 0.5, 0, 0], 'string', sprintf( 'iter %3d', iter), 'FontSize', 25)
    pause(.5);
end