function [pS] = putParamsStructure(params)

global alg

[~,names]      = readParam(alg.paramFile);

pS = struct();

for i = 1:length(names)
    pS.(char(names(i)))=params(i);
end


end
    