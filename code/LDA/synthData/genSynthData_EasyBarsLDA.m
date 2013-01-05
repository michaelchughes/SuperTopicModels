function [Data, Truth, Params] = genSynthData_EasyBarsLDA( dataParams )

defs.alph = 0.1;
defs.D = 100;
defs.Nd = 100;
defs.K = 10;
defs.V = 25;

Params = updateParamsWithUserInput( defs, dataParams );
% Add all fields of Params struct directly to the workspace
fnames = fieldnames( Params );
for ff = 1:length( fnames)
    eval( [fnames{ff} '= Params.' fnames{ff} ';'] );
end

M = sqrt(V);

ObsPr = ones( K, V );
FAVOR_COUNT = 15;

% Horizontal bars
for rowID = 1:K/2
    favorIDs = ( rowID -1 )*M + (1:M);
    ObsPr( rowID, favorIDs) = FAVOR_COUNT;
end

for colID = 1:K/2
    favorIDs = (colID-1) + (1:M:V);
    ObsPr( K/2+colID, favorIDs) = FAVOR_COUNT;
end


% ============================================  Generate data
for dd = 1:D
    pi = gamrnd( alph, 1, 1, K);
    pi = pi/sum(pi);
    
    zs = multinomial_many_draw( pi, Nd );
    
    ws = [];
    for kk = 1:K
       ws = [ws multinomial_many_draw( ObsPr(kk,:),  sum(zs==kk) )];
       
    end
    assert( length(ws)==Nd );
    Data(dd).words = ws;
    
    Truth.Topics{dd} = zs;
end
Truth.TopicObsPr = bsxfun(@rdivide,ObsPr, sum(ObsPr,2)  );
Truth.alpha = alph;