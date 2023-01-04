% This function simulates the cross sectional distribution of the economy.
% We focus on products and keep track of the (1) the gap, (2) the "name" of
% the owning firm and (3) the age of the owning firm. This function only
% simulates the cross-sectional distribution so we do not have to keep
% track of the whole evolution

% The stochastic process for innovations is modelled as realizations of a
% arrival times for the respective Poisson process. 

    
function [Gap, Names, Age] = SimulateCrossSection(I, z, tau, NumProducts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Set seed to replicate
rng(2204)

%% Set up objects for simulation
IntervalsWithinPeriods = 50;
NumberOfPeriods = 200;
time=linspace(0,NumberOfPeriods,NumberOfPeriods*IntervalsWithinPeriods-1);


%%%% State Variables
Gap = ones(NumProducts,1);        % quality gap
Names = [1:NumProducts]';         % Names of firms ("-1" refers to a firm not in our sample)
Age = 1/IntervalsWithinPeriods * ones(NumProducts,1);
ProductIndex = [1:NumProducts]';

% Initial arrival times

NextArrival_own_innovation = log(1-rand(NumProducts,1))/(-I);
NextArrival_cd = log(1-rand(NumProducts,1))/(-tau);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Simulate evolution of products  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:length(time)
        
    % Which innovation happens
    cd = abs(time(i)-NextArrival_cd) <= 1/IntervalsWithinPeriods;
    N_cd = sum(cd);
    
    own_innovation = abs(time(i)-NextArrival_own_innovation) <= 1/IntervalsWithinPeriods;
    own_innovation(cd == 1) = 0;
    N_own = sum(own_innovation);
        
    % Update Gap
    
    if N_own > 0
    Gap(own_innovation) = Gap(own_innovation) + 1;    
    end
    if N_cd > 0
    Gap(cd) = 1;    
    end
    
    % Update arrival times if innovation occurred
    NextArrival_own_innovation(own_innovation) = time(i) + log(1-rand(sum(N_own),1))/(-I);
    NextArrival_cd(cd) = time(i) + log(1-rand(sum(N_cd),1))/(-tau);
    
    % Update age
    Age = Age + 1/IntervalsWithinPeriods;

    % Update firm names and age for each product where creative destruction occurs   
    if N_cd > 0
    newNames = ones(N_cd,1)*(-1); %Names of produts, which get creatively destroyed
    newAge = ones(N_cd,1)*1/IntervalsWithinPeriods;
    N_cd_inc = sum(rand(N_cd,1) <= (1-z/tau)); %Share x/tau of creative destruction events are by incumbents
    N_cd_Entr = N_cd - N_cd_inc;
    if N_cd_inc > 0
    %Draw with replacement from existing firm name distribution 
    cd_inc_index = randsample(ProductIndex, N_cd_inc,'true'); % Find 
    newNames(1:N_cd_inc) = Names(cd_inc_index); 
    newAge(1:N_cd_inc) = Age(cd_inc_index); 
    end
    if N_cd_Entr > 0
    newNames(N_cd_inc+1:N_cd) = [max(Names)+1:max(Names)+1+N_cd_Entr-1];    
    end  
    Names(cd) = newNames;
    Age(cd) = newAge;
    end
    
end    

end

