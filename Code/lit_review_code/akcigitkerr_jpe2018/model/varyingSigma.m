%% Generalized model with different values of sigma
function [] = varyingSigma()

fprintf('Solving the model under different sigmas...')
global alg
initAlg('Baseline');
gridSigma = [0.000000000001,0.2,.4,0.5];

for i = 1:length(gridSigma)

    pS = readParam(alg.paramFile);   
    pS.sigma = gridSigma(i); 
    p  = putParams(pS);
    [A(i),Bn(:,i),eq,flag] = solveEq(p);     % solve model
    innovation(:,i) = eq.xn;
end

save(['Mat files' filesep 'resultBaselineVarySigma.mat'],'A','Bn','innovation');
fprintf('Done!\n')

end