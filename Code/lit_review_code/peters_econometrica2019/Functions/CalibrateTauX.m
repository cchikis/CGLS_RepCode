% This function calibrates x and tau to match the entry rate and the
% employment life-cycle


function [y_sol, fval] = CalibrateTauX(Moments)

Minobject = @(y) TauXObjectiveFunction(y,Moments);
y_init = [0.2, 0.25]; 

%Constraint x < tau: A*X  <= B
A = [1, -1]; B = 0;
% Lower bound on tau and x for non-negative 
LB = [0, 0];
options = optimoptions('fmincon','Display','off');
[y_sol, fval] = fmincon(Minobject, y_init, A, B, [], [], LB, [], [], options);
end



function diff = TauXObjectiveFunction(y,Moments)

x = y(1);
tau = y(2);

% (1) Calculate Entry Rate in model:
ER = x / (log(tau/(tau-x)));

% (2) Calculate employment lifecycle in model:

lnnOLD = lnNbyAge(7.5, x, tau);     % log Sales for for 7.5 year old firms
lnnYoung = lnNbyAge(0.5, x, tau);   % log Sales for for 0.5 year old firms

% Employment lifecycle = sales growth - markup growth
LClnEmpl.Model = lnnOLD - lnnYoung - Moments.LCMU;


% Set up objective function
diff = (Moments.LCEmpl - LClnEmpl.Model).^(2) + (Moments.EntryRate - ER).^(2);
 
end


