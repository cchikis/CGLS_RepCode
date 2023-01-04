%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: emp_share_sim.m
% Author: Craig A. Chikis
% Date: 09/15/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = emp_share_sim(m) 

    lt10 = sum(m.panel_save.emp(m.panel_save.Age < 11));
    lt5 = sum(m.panel_save.emp(m.panel_save.Age < 6));
    tot = sum(m.panel_save.emp); 


    m.share10 = lt10/tot;
    m.share5 = lt5/tot;


end 