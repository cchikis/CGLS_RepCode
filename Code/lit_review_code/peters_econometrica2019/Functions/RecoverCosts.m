% This function recovers the innovation and entry costs given the
% equilibrium innovation flow intensities x, I and z 



function param = RecoverCosts(EqObjects,param)

l = param.lambda;
rho = param.rho;
zeta = param.zeta;


tau = EqObjects.tau;
I = EqObjects.I;
x = EqObjects.x;

% Innovation Costs
Lambda = tau/( l*tau + (l-1)*I );
phiI = param.L^(-1)*(zeta * Lambda * ( l/(l-1) * I^(zeta-1) * (rho + tau) + I^zeta ) + ...
         (tau - (zeta-1)/zeta*x) * zeta * (I^(zeta-1) * (rho + tau) + (zeta-1)/zeta * I^zeta) / (rho + tau - (zeta-1)/zeta*x));
     
% Entry costs
phiz = (rho + tau - (zeta-1)/zeta*x) / (I^(zeta-1) * (rho + tau) + (zeta-1)/zeta * I^zeta) * (1/zeta) * phiI;

% Expansion barriers
phix = x^(zeta-1) * zeta * phiz;

% Save output
param.phiz = phiz;
param.phix = phix;
param.phiI = phiI;


end

