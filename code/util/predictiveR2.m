function pR2 = predictiveR2( truey, esty )
%  Compute Predictive R^2 measure for estimator
%  Syntax:  pR2 = predictiveR2( truey, esty )
%  Inputs:
%      truey: ground truth scalar values
%              Nx1 column vector
%      esty:  estimates of y
%              Nx1 column vector

varEst = sum( (truey - esty).^2 );
varTrue = sum( (truey - mean(truey) ).^2 );
pR2 = 1 - varEst/varTrue;