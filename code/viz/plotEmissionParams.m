function [] = plotEmissionParams( jobID, taskID, queryIter, varargin )

Prefs.doAlign = true;
Prefs.est_labels = [];
Prefs.figHandle = [];
Prefs.SamplerOutput = [];
Prefs = updateParamsWithUserInput( Prefs, varargin );

if isempty( Prefs.SamplerOutput )
    X = loadSamplerOutput(jobID, taskID );
else
    X = Prefs.SamplerOutput;
end

if exist( 'queryIter', 'var')
    [~,bestID] = min( abs( X.iters.Psi - queryIter ) );
    Psi = X.Psi(bestID);
    queryIter = X.iters.Psi(bestID);
else
    Psi = X.Psi(end);
    queryIter = X.iters.Psi(end);
end

% ============================== Plot square topic images
Nkw = Psi.TWSuffStats.Nkw;

K = size( Nkw,1);
V = size( Nkw,2);
R = 2;
C = 5;
if R*C ~= K
    R = K;
    C = 1;
end


for kk = 1:K
    subplot( R, C, kk );
    phi = Nkw( kk, :) + Psi.beta;
    Phi(kk,:) = phi / sum(phi);
end

% Align to true topics
if Prefs.doAlign
    if isempty( Prefs.est_labels )
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
    else
        est_labels = Prefs.est_labels;
    end
end

if isempty( Prefs.figHandle )
    hh=figure();
    set( hh, 'Visible', 'off' );
else
    hh = figure( Prefs.figHandle );
    clf;
end

for kk = 1:K
    subplot( R, C, est_labels(kk) );
    sqIm = reshape( Phi(kk,:), sqrt(V), sqrt(V) );
    imagesc( sqIm, [0 .5] );
    
    if isfield( Psi, 'eta')
        title(  sprintf( '% .1f', Psi.eta(kk) ), 'FontSize', 25 );
    end
    axis image;
end
colormap hot;

drawnow;
set( hh, 'Visible', 'on' );
hold off;