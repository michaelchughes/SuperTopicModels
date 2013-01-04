function [y] = randmvnormal(  mu,  invSigma )
%MVNRAND - generate multivariate normal random values 
% mu :  K x 1 column vector
% invSigma : K x K inverse covariance matrix 

K = length(mu);
L = chol(invSigma);
y = L\randn(K,1) + mu;
end
