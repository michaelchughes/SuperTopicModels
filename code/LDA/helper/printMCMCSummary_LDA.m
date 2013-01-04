function [] = printMCMCSummary_LDA( iter, Psi, logPr, algParams )

fprintf( '\t % 5d/%d after %6.0f sec | logPr % .2e \n', ...
     iter, algParams.Niter, toc, logPr.all );

end