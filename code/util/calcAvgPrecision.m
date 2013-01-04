function [avgPrecision] = calcAvgPrecision( prTrue, trueLabels, doVerbose)
%  Computes average precision for single binary classifier
% Input:
%     prTrue: Predictive output of binary classifier for entire corpus
%             D-length column vector, where D is num documents in corpus
%             prTrue(d) = scalar prob. document d is positive instance
%                 does not have to be true probability e.g. in (0,1)
%                  can just be a scalar score, more positive = higher rank
%     trueLabels:  Ground truth binary labels for entire corpus
%             D-length column vector
%             trueLabels(d) = 1 if doc d is truly a positive instance
%                             0 otherwise
% Output:
%    avgPrecision = scalar valued average precision for this classifier
%      computed following formula at this URL:
%      http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-ranked-retrieval-results-1.html
% Toy Example:
%    prTrue     = [.5 .9 .2 .7      .3 .1 .6 .05];
%    trueLabels = [ 1  1  1  1       0  0  0   0];
%    AP = calcAvgPrecision( prTrue, trueLabels )
%    should yield AP = .85...
if ~exist( 'doVerbose', 'var' )
    doVerbose = 0;
end


[rankedProbs, idxs] = sort( prTrue, 'descend' );
 rankedLabels = trueLabels(idxs);
%  rankedTrueSum(k) gives # true examples found in top k ranked results
rankedTrueSum = cumsum( rankedLabels );
Precision = zeros( 1, sum(trueLabels ) );
for k = 1:sum( trueLabels )
   kthIdx = find( rankedTrueSum == k, 1, 'first');
   Precision(k) =  rankedTrueSum(kthIdx) / kthIdx;   
   if doVerbose
       disp( [ 'Retrieving the ' num2str(k) '-th relevant doc.    Precision=' num2str(Precision(k)) ] );
   end
end

avgPrecision = mean( Precision );