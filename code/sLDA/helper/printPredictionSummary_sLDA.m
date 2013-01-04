function [] = printPredictionSummary_sLDA( iter, Predict, model, algP )

if model.RegM.doRegression
fprintf( '\t % 5d/%d \t\t\t\t\t TRAIN pR2 %.3f | TEST pR2 %.3f  \n', ...
     iter, algP.Niter, Predict.train.pR2, Predict.test.pR2 );
elseif model.RegM.doBinary || model.RegM.doMultiClass
fprintf( '\t % 5d/%d \t\t\t\t\t TRAIN acc %.3f | TEST acc %.3f  \n', ...
     iter, algP.Niter, Predict.train.acc, Predict.test.acc );
end
    

end