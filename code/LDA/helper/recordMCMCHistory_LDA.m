function ChainHist = recordMCMCHistory_LDA( n, outParams, ChainHist, Psi, logPr )
% Save current state of sampler

% -------------------------------------------------- update logPr trace
if n == 1 || rem( n, outParams.logPrEvery ) == 0
    if isfield( ChainHist, 'logPr' )
        dC = length( ChainHist.logPr ) + 1;
        ChainHist.logPr(dC) = logPr;
    else
        ChainHist = struct();
        ChainHist.logPr = logPr;
        dC = 1;
    end
    
    ChainHist.iters.logPr(dC) = n;
    ChainHist.times.logPr(dC) = toc;
end

% -------------------------------------------------- save current config
if ( n==1 || rem( n, outParams.saveEvery)==0 )
    if isfield( ChainHist, 'Psi' )
        storeCount = length( ChainHist.Psi ) + 1;
        ChainHist.Psi( storeCount ) = Psi;
    else
        storeCount = 1;
        ChainHist.Psi = Psi;
    end
    
    ChainHist.iters.Psi( storeCount) = n;
    ChainHist.times.Psi( storeCount) = toc;
    
    ChainHist.RandSeed(storeCount).matlabPRNGState   = RandStream.getGlobalStream.State;    
end

end % MAIN FUNCTION