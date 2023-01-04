%% Main Driver for Estimation
clear all;
close all;
clc;

disp('------------');
disp(' ESTIMATION ')
disp('------------');


global alg;
initAlg('Baseline');

[params,~] = readParam(alg.paramFile);


alg.pScale  = (cell2mat(struct2cell(params)))';
alg.nParams = length(alg.pScale);
                   
alg.fixParams  = [2 4 6 7 10 11 13 15]; 
alg.estParams  = setdiff(1:alg.nParams,alg.fixParams);
alg.nFixParams = length(alg.fixParams);
alg.nEstParams = length(alg.estParams);

alg.plb    = zeros(1,alg.nParams) + 0.000001;  
alg.plb(2) = 1.00001;
alg.plb(4) = 1.00001;
alg.pub    = Inf*ones(1,alg.nParams); 
alg.pub([6 7 8 9 10 12]) = 1.0; 


alg.simann = 0;
init_scale = 0.25;
start = ones(1,length(alg.estParams));
alg.bestval = Inf;

    
if (alg.simann == 1)
    [parfin,~] = anneal0(@estimFunc,start,init_scale);
 
else
    mopts = optimset('Display','iter','MaxFunEvals',20000,'MaxIter',10000);
    [parfin,~] = fminsearch(@estimFunc,start,mopts);
end



params = alg.pScale;
params(alg.estParams) = parfin.*alg.pScale(alg.estParams); 
