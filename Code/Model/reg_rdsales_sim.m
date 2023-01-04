%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: reg_rdsales_sim.m
% Author: Craig A. Chikis
% Date: 08/27/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = reg_rdsales_sim(m)

    panel = m.panel_out; 

    panel.bucket = nominal(panel.bucket);
    panel.bucket_tm1 = nominal(panel.bucket_tm1);
    panel.Time = nominal(panel.Time); 
    reg_rdsales_profit_tm1_bucket = fitlm(panel, 'rdsales_w ~ bucket_tm1 + Time'); 


    m.reg_rdsales_profit_tm1_bucket = reg_rdsales_profit_tm1_bucket; 

end