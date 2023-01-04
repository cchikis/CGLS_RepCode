%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: run_paper_new.m
% Author: Craig A. Chikis
% Date: 11/22/2022
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []  = run_paper_new(nlopt_dir) 

    data1(); 
    close all 
    run_benchmark(); 
    close all 
    run_multistep(); 
    close all

    [retcode, ~] = system("Rscript Code/Data/outregs_model_data_v4.R");

    fig_code_1();
    close all 
    moment_sensitivity();
    close all  
    kap_sig();
    close all 

    robustness_code();
    close all

    vanishing_QCU_fig();
    close all

    transition_solve();
    close all

    SCUexplore();
    close all 

    error_analysis(); 
    close all

    wrapper_lit_review(); 
    close all

    run_LMS();
    close all  
    lms_compare(); 
    close all

    financial_frictions_fig(); 
    close all

    if (~strcmp(nlopt_dir, ""))

        cd Code/lit_review_code/LMS_optimaleta/
        LMS_opt_eta(nlopt_dir);
        close all 

        cd ../../../

    end

end