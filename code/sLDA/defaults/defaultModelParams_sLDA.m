function model = defaultModelParams_sLDA( Data )
K = 10;
model.nTopics = K;

V = max( horzcat( Data(:).words ) );
model.TopicObsM.nVocab = V;
model.TopicObsM.beta = min(1, 10/V );

model.DocTopicM.alpha = 1;

model.RegM.eta.type = 'Gaussian';
model.RegM.eta.mean = zeros(K,1);
model.RegM.eta.var  = 100;

model.RegM.lambda.type = 'Gamma';
model.RegM.lambda.a = 1/10;
model.RegM.lambda.b = 1/10;

ys = vertcat( Data(:).y );
us = unique(ys);
doBinary = false;
doMultiClass = false;
doRegression = false;
nClass= 0;
if length( us)==2 && us(1)==0 && us(2)==1
    doBinary = true;
elseif length(us) < length(ys) && us(1)==int32(us(1)) && us(2)==int32(us(2))
    doMultiClass = true;
    nClass = length(us);
else
    doRegression = true;
end
model.RegM.doBinary = doBinary;
model.RegM.doMultiClass =doMultiClass;
model.RegM.doRegression = doRegression;
model.RegM.nClass = nClass;