#################################################
# Filename: outregs_model_data_v4.R
# Author: Craig A. Chikis
# Date: 10/19/2022
# Note(s):
#################################################
rm(list = ls(all.names = TRUE))
gc()


set.seed(20220719)

packages = c("lfe", "tidyverse")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    install.packages(setdiff(packages, rownames(installed.packages())))
}

library(lfe)
library(tidyverse)


cstat_out = readRDS("Output/Store_Data/cstat_out_three_dropm1.rds")
ds_iter_read_in = function(ds_name) {


    patent_out_CGLS = map(.x = c(patent_data_out_p90 = paste0("patent_data_out_p90"),
                                 patent_data_out_p80 = paste0("patent_data_out_p80")),
                          .f = ~readRDS(paste0("Output/Store_Data/patent_out_", .x, "_", .x, "_", ds_name, ".rds")))

}
data_in = list(CGLS = ds_iter_read_in("kpst"))
data_in_nber = list(CGLS = ds_iter_read_in("nber_nber"))
kpss1 = read_csv("Output/Store_Data/kpss1.csv")

params = read_csv("Output/Store_Data/parameters130.csv")

panel_out_wrap = list(panel_out = read_csv("Output/Store_Data/panel_out_130.csv")) %>%
                      map(.f = ~mutate(.x, rdsales_perc_w = rdsales_w*100,
                                           profitgrowth_perc_w = profitgrowth_w*100))

patent_model_wrap = list(patent_model = list(p80_patent = read_csv("Output/Store_Data/qualitysale_tm1_130_80.csv"),
                                             p90_patent = read_csv("Output/Store_Data/qualitysale_tm1_130_90.csv"))) %>%
                                                       map(.f = ~map(.x, .f = ~(.x %>% mutate(quality_perc = quality*100))))
citation_wrap = list(citation_model = read_csv("Output/Store_Data/citation_track2_130.csv"))


regression_CSTAT_model = function(LHS, data, patent = FALSE, fe = "Time") {

    data = data %>%
        mutate(Time = as.factor(Time),
               Firm = as.factor(Firm))
    if (!patent) { 
        data = data %>%
            mutate(Industry = as.factor(Industry))
    }

    if (!patent) {
        data = data %>%
            mutate(bucket_tm1 = as.factor(bucket_tm1),
                   bucket_tm1_2 = as.factor(bucket_tm1_2))
    } else {
        data = data %>%
            mutate(bucket_tm1_patent = as.factor(bucket_tm1_patent),
                   bucket_tm1_patent_2 = as.factor(bucket_tm1_patent_2),
                   CPC = as.factor(CPC))
    }

    if (!patent) {
        bucket_tm1 = "bucket_tm1"
        bucket_tm1_2 = "bucket_tm1_2"
    } else {
        bucket_tm1 = "bucket_tm1_patent"
        bucket_tm1_2 = "bucket_tm1_patent_2"
    }
    dummy5 = felm(formula(paste0(LHS, " ~ ", bucket_tm1, " | ", fe, " | 0 | Firm")), data)
    dummy2 = felm(formula(paste0(LHS," ~ ", bucket_tm1_2, " | ", fe, " | 0 | Firm")), data)
    
    return(list(dummy5 = dummy5, dummy2 = dummy2))
}
regs_cstat_model = map(.x = panel_out_wrap, 
                       .f = ~map(.x = c(profitgrowth_w = "profitgrowth_perc_w", rdsales_w = "rdsales_perc_w"),
                                 .f = regression_CSTAT_model, data = .x, patent = FALSE, fe = "Time"))



model_patent_run = function(qualitysale_tm1) {
    reg_patent_model_Time_CPC = regression_CSTAT_model("quality_perc", qualitysale_tm1, patent = TRUE, fe = "Time:CPC")

}
regs_patent = map(.x = patent_model_wrap, .f = ~map(.x = .x, .f = ~model_patent_run(.x)))




cdfiter = function(data) {
    data = data %>% 
        filter(!is.na(fcit020)) %>%
        distinct(wku, .keep_all = TRUE) %>%
        summarise(lt0 = sum(fcit020 <= 0)/sum(!is.na(fcit020)),
                  lt10 = sum(fcit020 <= 10)/sum(!is.na(fcit020))) %>%
        ungroup()
}
cdf_model = map(.x = citation_wrap, .f = cdfiter)
cdf_data = cdfiter(kpss1)




gammatable = paste0("\\begin{tabular}{ccc|clcc}  \\hline \\hline  \\multicolumn{2}{c}{Parameter estimate} & & & ",
                    " \\multicolumn{3}{c}{Moments used in the estimation} \\\\ \\cline{1-2} \\cline{5-7}  ", 
                    " \\multicolumn{1}{l}{Parameter} &  \\multicolumn{1}{c}{Value} & & &  ",
                    " \\multicolumn{1}{l}{Description} &  \\multicolumn{1}{c}{Model} & \\multicolumn{1}{c}{Data} \\\\",
                    " \\hline $\\chi$ & ", sprintf('%0.3f', params$values[params$parameter == "zeta"]), " & ",
                    " &  & $\\Pr\\{\\text{Cit.} = 0\\}$ & ", sprintf('%0.2f', 100*cdf_model$citation_model$lt0),
                                "\\% ", " & ", sprintf('%0.2f', 100*cdf_data$lt0), "\\% \\\\  &  &  &  & ", 
                    " $\\Pr\\{\\text{Cit.} \\leq 10 \\}$ & ", sprintf('%0.2f',  100*cdf_model$citation_model$lt10), "\\%",
                            " & ", sprintf('%0.2f', 100*cdf_data$lt10), "\\%", " \\\\ \\hline \\end{tabular}")
cat(gammatable, file = paste0("Output/LaTeX_Output/gammatable", ".txt"))

model_cstat_iter = function(regs_cstat_model, add_on) {

    profitgrowth_coef_model = sprintf('%0.2f', regs_cstat_model$panel_out$profitgrowth_w$dummy5$coefficients)
    profitgrowth_se_model = sprintf('%0.2f', regs_cstat_model$panel_out$profitgrowth_w$dummy5$cse)
    profitgrowth_N_model = sprintf('%0.0f', regs_cstat_model$panel_out$profitgrowth_w$dummy5$N)

    profitgrowth_coef_data = sprintf('%0.2f', 100*cstat_out$regs$reg_profitgrowth_tm1$coefficients)
    profitgrowth_se_data = sprintf('%0.2f', 100*cstat_out$regs$reg_profitgrowth_tm1$cse)
    profitgrowth_N_data = sprintf('%0.0f', cstat_out$regs$reg_profitgrowth_tm1$N)


    rdsales_coef_model = sprintf('%0.2f', regs_cstat_model$panel_out$rdsales_w$dummy5$coefficients)
    rdsales_se_model = sprintf('%0.2f', regs_cstat_model$panel_out$rdsales_w$dummy5$cse)
    rdsales_N_model = sprintf('%0.0f', regs_cstat_model$panel_out$rdsales_w$dummy5$N)

    rdsales_coef_data = sprintf('%0.2f', 100*cstat_out$regs$reg_rdsales_tm1$coefficients)
    rdsales_se_data = sprintf('%0.2f', 100*cstat_out$regs$reg_rdsales_tm1$cse)
    rdsales_N_data = sprintf('%0.0f', cstat_out$regs$reg_rdsales_tm1$N)




    quintile_reg = paste0("\\begin{adjustbox}{center} \\begin{tabular}{lcc} ",
                            " \\hline \\hline  & \\hspace{0.3em} Profit growth \\hspace{0.3em} & \\hspace{0.3em} R\\&D-to-sales \\hspace{0.3em}  \\\\ ",
                            " \\hline & \\multicolumn{2}{c}{A. Model, effects relative to top profit quintile} \\\\ ", 
                " \\cline{2-3} Second quintile \\hspace{1.6em} & ", profitgrowth_coef_model[1], " (", profitgrowth_se_model[1], ")", " & ",
                rdsales_coef_model[1], " (", rdsales_se_model[1], ")", " \\\\ ", 

                " Third quintile & ", profitgrowth_coef_model[2], " (", profitgrowth_se_model[2], ")", " & ",
                rdsales_coef_model[2], " (", rdsales_se_model[2], ")", " \\\\ ", 

                " Fourth quintile & ", profitgrowth_coef_model[3], " (", profitgrowth_se_model[3], ")", " & ",
                rdsales_coef_model[3], " (", rdsales_se_model[3], ")", " \\\\ ", 

                " Smallest quintile & ", profitgrowth_coef_model[4], " (", profitgrowth_se_model[4], ")", " & ",
                rdsales_coef_model[4], " (", rdsales_se_model[4], ")", " \\\\ ", 

                " \\cline{2-3} & \\multicolumn{2}{c}{B. Data, effects relative to top profit quintile} \\\\ ",

                " \\cline{2-3} Second quintile & ", profitgrowth_coef_data[1], " (", profitgrowth_se_data[1], ")", " & ",
                rdsales_coef_data[1], " (", rdsales_se_data[1], ")", " \\\\ ", 

                " Third quintile & ", profitgrowth_coef_data[2], " (", profitgrowth_se_data[2], ")", " & ",
                rdsales_coef_data[2], " (", rdsales_se_data[2], ")", " \\\\ ", 

                " Fourth quintile & ", profitgrowth_coef_data[3], " (", profitgrowth_se_data[3], ")", " & ",
                rdsales_coef_data[3], " (", rdsales_se_data[3], ")", " \\\\ ", 

                " Smallest quintile & ", profitgrowth_coef_data[4], " (", profitgrowth_se_data[4], ")", " & ",
                rdsales_coef_data[4], " (", rdsales_se_data[4], ")", " \\\\ ", 

                "\\hline  \\end{tabular} \\end{adjustbox}")

    cat(quintile_reg, file = paste0("Output/LaTeX_Output/quintile_new_", add_on, ".txt"))

    return(invisible(NULL)) 
}
list(regular = regs_cstat_model) %>%
    walk2(.y = names(.), .f = model_cstat_iter)

iter_patent_data = function(data_in, name) { 

    patent_table = paste0("\\begin{tabular}{lcc} \\hline \\hline & \\hspace{0.3em}  Top 10\\% patent share & \\hspace{0.3em}  Top 20\\% patent share \\\\ ",
                        " \\hline & \\multicolumn{2}{c}{A. Model, effects relative to larger half} \\\\ ",
                        " \\cline{2-3} Smaller half \\hspace{0.3em} & ", sprintf('%0.2f', regs_patent$patent_model$p90_patent$dummy2$coefficients), 
                            " (", sprintf('%0.2f', regs_patent$patent_model$p90_patent$dummy2$cse), ") & ", 
                            sprintf('%0.2f', regs_patent$patent_model$p80_patent$dummy2$coefficients), " (", 
                            sprintf('%0.2f', regs_patent$patent_model$p80_patent$dummy2$cse), ") \\\\ ",
                        " \\cline{2-3} & \\multicolumn{2}{c}{B. Data, effects relative to larger half} \\\\ ",
                        " \\cline{2-3} Smaller half \\hspace{0.3em}  & ", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p90$coefficients$reg2profdummy_tm1_2), " (", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p90$cse$reg2profdummy_tm1_2), ") & ", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p80$coefficients$reg2profdummy_tm1_2), " (", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p80$cse$reg2profdummy_tm1_2), ")", " \\\\ ", 
                        " \\hline \\end{tabular}") 


    patent_table2 = paste0("\\begin{tabular}{lcc} \\hline \\hline & \\hspace{0.3em}  Top 10\\% patent share & \\hspace{0.3em}  Top 20\\% patent share \\\\ ",
                        " \\hline & \\multicolumn{2}{c}{A. Model, effects relative to larger half} \\\\ ",
                        " \\cline{2-3}  Second quintile \\hspace{0.3em} & ", sprintf('%0.2f', regs_patent$patent_model$p90_patent$dummy5$coefficients[1]), 
                            " (", sprintf('%0.2f', regs_patent$patent_model$p90_patent$dummy5$cse[1]), ") & ", 
                            sprintf('%0.2f', regs_patent$patent_model$p80_patent$dummy5$coefficients[1]), " (", 
                            sprintf('%0.2f', regs_patent$patent_model$p80_patent$dummy5$cse[1]), ") \\\\ ",

                        " Third quintile \\hspace{0.3em} & ", sprintf('%0.2f', regs_patent$patent_model$p90_patent$dummy5$coefficients[2]), 
                            " (", sprintf('%0.2f', regs_patent$patent_model$p90_patent$dummy5$cse[2]), ") & ", 
                            sprintf('%0.2f', regs_patent$patent_model$p80_patent$dummy5$coefficients[2]), " (", 
                            sprintf('%0.2f', regs_patent$patent_model$p80_patent$dummy5$cse[2]), ") \\\\ ",

                        " Fourth quintile \\hspace{0.3em} & ", sprintf('%0.2f', regs_patent$patent_model$p90_patent$dummy5$coefficients[3]), 
                            " (", sprintf('%0.2f', regs_patent$patent_model$p90_patent$dummy5$cse[3]), ") & ", 
                            sprintf('%0.2f', regs_patent$patent_model$p80_patent$dummy5$coefficients[3]), " (", 
                            sprintf('%0.2f', regs_patent$patent_model$p80_patent$dummy5$cse[3]), ") \\\\ ",

                        " Smallest quintile \\hspace{0.3em} & ", sprintf('%0.2f', regs_patent$patent_model$p90_patent$dummy5$coefficients[4]), 
                            " (", sprintf('%0.2f', regs_patent$patent_model$p90_patent$dummy5$cse[4]), ") & ", 
                            sprintf('%0.2f', regs_patent$patent_model$p80_patent$dummy5$coefficients[4]), " (", 
                            sprintf('%0.2f', regs_patent$patent_model$p80_patent$dummy5$cse[4]), ") \\\\ ",
                        
                        " \\cline{2-3} & \\multicolumn{2}{c}{B. Data, effects relative to larger half} \\\\ ",
                        " \\cline{2-3}  Second quinitile \\hspace{0.3em}  & ", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p90$coefficients$reg2profdummy_tm1[1]), " (", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p90$cse$reg2profdummy_tm1[1]), ") & ", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p80$coefficients$reg2profdummy_tm1[1]), " (", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p80$cse$reg2profdummy_tm1[1]), ")", " \\\\ ",


                          " Third quintile \\hspace{0.3em}  & ", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p90$coefficients$reg2profdummy_tm1[2]), " (", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p90$cse$reg2profdummy_tm1[2]), ") & ", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p80$coefficients$reg2profdummy_tm1[2]), " (", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p80$cse$reg2profdummy_tm1[2]), ")", " \\\\ ",


                          " Fourth quintile \\hspace{0.3em}  & ", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p90$coefficients$reg2profdummy_tm1[3]), " (", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p90$cse$reg2profdummy_tm1[3]), ") & ", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p80$coefficients$reg2profdummy_tm1[3]), " (", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p80$cse$reg2profdummy_tm1[3]), ")", " \\\\ ",

                          " Smallest quintile \\hspace{0.3em}  & ", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p90$coefficients$reg2profdummy_tm1[4]), " (", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p90$cse$reg2profdummy_tm1[4]), ") & ", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p80$coefficients$reg2profdummy_tm1[4]), " (", 
                        sprintf('%0.2f', 100*data_in$patent_data_out_p80$cse$reg2profdummy_tm1[4]), ")", " \\\\ ",

                        " \\hline \\end{tabular}") 


    cat(patent_table, file = paste0("Output/LaTeX_Output/patent_table", name,".txt"))
    cat(patent_table2, file = paste0("Output/LaTeX_Output/patent_table_2", name,".txt"))



    return(invisible(NULL))

}
walk2(.x = data_in, .y = names(data_in), .f = iter_patent_data)


s_nber = paste0("\\begin{tabular}{lcc} \\hline \\hline & Top 10\\% patent share & Top 20\\% patent share \\\\ ", 
                " \\hline & \\multicolumn{2}{c}{A. Effects relative to the top profit half} \\\\ ",
                " \\cline{2-3} Smaller half & ", 
                sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p90$coefficients$reg2profdummy_tm1_2), 
                " (", sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p90$cse$reg2profdummy_tm1_2), ")", 
                " & ", sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p80$coefficients$reg2profdummy_tm1_2),
                " (",  sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p80$cse$reg2profdummy_tm1_2), ")", 
                " \\\\ ", 
                " \\cline{2-3} & \\multicolumn{2}{c}{B. Effects relative to the top profit quintile} \\\\ ",
                " \\cline{2-3} Second quintile & ", 
                sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p90$coefficients$reg2profdummy_tm1[1]), 
                " (", sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p90$cse$reg2profdummy_tm1[1]), ") ", 
                " & ",  sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p80$coefficients$reg2profdummy_tm1[1]), 
                " (", sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p80$cse$reg2profdummy_tm1[1]), ")  \\\\ ", 

                " Third quintile & ", 
                sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p90$coefficients$reg2profdummy_tm1[2]), 
                " (", sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p90$cse$reg2profdummy_tm1[2]), ") ", 
                " & ",  sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p80$coefficients$reg2profdummy_tm1[2]), 
                " (", sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p80$cse$reg2profdummy_tm1[2]), ")  \\\\ ", 

                " Fourth quintile & ", 
                sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p90$coefficients$reg2profdummy_tm1[3]), 
                " (", sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p90$cse$reg2profdummy_tm1[3]), ") ", 
                " & ",  sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p80$coefficients$reg2profdummy_tm1[3]), 
                " (", sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p80$cse$reg2profdummy_tm1[3]), ")  \\\\ ", 


                " Smallest quintile & ", 
                sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p90$coefficients$reg2profdummy_tm1[4]), 
                " (", sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p90$cse$reg2profdummy_tm1[4]), ") ", 
                " & ",  sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p80$coefficients$reg2profdummy_tm1[4]), 
                " (", sprintf('%0.2f', 100*data_in_nber$CGLS$patent_data_out_p80$cse$reg2profdummy_tm1[4]), ")  \\\\ ", 
                " \\hline \\end{tabular}")

cat(s_nber, file = "Output/LaTeX_Output/nber_tab.txt")