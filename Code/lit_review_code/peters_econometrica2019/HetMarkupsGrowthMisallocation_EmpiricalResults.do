**************************************************************** 
/*
THIS FILE CONTAINS THE EMPIRICAL ANALYSIS FOR THE PAPER 
"HETEROGENEOUS MARKUPS, GROWTH AND ENDOGENOUS MISALLOCATION" 
*/
****************************************************************   


clear
set more off

use "HetMarkupsGrowthMisallocation_FirmlevelData.dta"

******************************************************************************************
****************** Table 1: DESCRIPTIVE STATISTICS
******************************************************************************************

preserve 

keep if Sample == 1
drop if year >= 1998
drop if year == 1990

egen TotL = sum(laborforce), by(year)
egen EntryLtmp = sum(laborforce), by(year Entry)
replace EntryLtmp = 0 if Entry == 0
egen EntryL = max(EntryLtmp), by(year)
gen EntryLShare = EntryL/TotL

egen TotSales = sum(sales), by(year)
egen EntryStmp = sum(sales), by(year Entry)
replace EntryStmp = 0 if Entry == 0
egen EntryS = max(EntryStmp), by(year)
gen EntrySalesShare = EntryS/TotSales

egen ExitLtmp = sum(laborforce), by(year Exit)
replace ExitLtmp = 0 if Exit == 0
egen ExitL = max(ExitLtmp), by(year)
gen ExitLShare = ExitL/TotL

egen ExitStmp = sum(sales), by(year Exit)
replace ExitStmp = 0 if Exit == 0
egen ExitS = max(ExitStmp), by(year)
gen ExitSalesShare = ExitS/TotSales

collapse (mean) lmean = laborforce (p25) l25 = laborforce (p50) l50 = laborforce (p90) l90 = laborforce ///
(mean) Entry  EntryLShare  EntrySalesShare Exit  ExitLShare ExitSalesShare

export excel using "Tables/Table1_DescriptiveTable.xls", firstrow(variables) replace

restore


******************************************************************************************
****************** Table 2: The Life-cycle of markups in Indonesia
******************************************************************************************

preserve

** Generate quantiles of capital-labor ratio
xtile kl_pctile_50 =lnkl, n(50)

eststo clear
** Basic profile
eststo: quietly xi: areg ln_sl_corr cohortage i.year if Sample == 1, absorb(mainproduct) robust
** Controlling for ln kl
eststo: quietly xi: areg ln_sl_corr cohortage lnkl i.year if Sample == 1, absorb(mainproduct) robust
** Controlling for ln kl non-parametrically
eststo: quietly xi: areg ln_sl_corr cohortage i.kl_pctile_50 i.year if Sample == 1, absorb(mainproduct) robust
** Balanced panel
eststo: quietly xi: areg ln_sl_corr cohortage lnkl i.year if NoYears == 2000 - cohort + 1 & Sample == 1, absorb(mainproduct) robust

esttab using "Tables/Table2_Markuplifecycle.tex", replace ///
se nolabel r2 star(* 0.10 ** 0.05 *** 0.01) ///
keep(cohortage lnkl) ///
indicate("Year & Industry FE = _Iyear*" "k/l percentiles = _Ikl*")

restore


******************************************************************************************
****************** Figure 3: Export data for "Life-cycle of markups in Indonesia"
******************************************************************************************


preserve


** Take out year and product fixed effects
quietly: xi: areg ln_sl i.year, absorb(mainproduct)
predict rellnmarkup, residuals

** Calculate mean and standard deviation by age of the cohort
collapse (mean) meanmu = rellnmarkup (semean) semu = rellnmarkup (count) N = psid,  by(cohortage)
drop if missing(cohortage)

** Normalize markups of entrants to zero
gen varyoung = meanmu if cohortage == 0
egen tmp = max(varyoung)
gen meanmu_rel = meanmu - tmp

** Get confidence intervals
gen muH =  meanmu_rel + 1.645*semu
gen muL =  meanmu_rel - 1.645*semu

drop meanmu semu tmp varyoung 


** Export data to plot figure in MATLAB
export delimited using "DataInputs/Data_MarkupLifeCycle.csv", replace

restore


******************************************************************************************
****************** Figure 4: Export data on markups and size for "Non-targeted moments"
******************************************************************************************

preserve 

quietly: xi: areg ln_sl i.year, absorb(mainproduct)
predict rellnmarkup, residuals

quietly: xi: areg lnlaborforce i.year, absorb(mainproduct)
predict rellnl, residuals

** Calculate averages of employment and markups by cohort and size of the cohort
collapse (mean) meanmu = rellnmarkup meanlnl =  rellnl,  by(cohortage)
 
** Normalize markups and size of entrants to zero 

drop if missing(cohortage)

foreach var of varlist meanmu meanlnl {
gen varyoung = `var' if cohortage == 0
egen tmp = max(varyoung)
gen `var'_rel = `var' - tmp
drop tmp varyoung
}

** Export data to plot figure in MATLAB
export delimited cohortage meanmu_rel meanlnl_rel using "DataInputs/Data_IndonesiaPanel.csv", novarnames replace

restore


******************************************************************************************
****************** Figure 4: Export data on survival for "Non-targeted moments"
******************************************************************************************

preserve

keep if missing(cohortage) == 0
collapse (count) N = psid, by(cohortage year)

gen EntrySizeCohort = .
gen Cohort = .

set more off
foreach t of numlist 1991(1)2000 {

* Number of entrants per year
gen tmpvar = N if cohortage == 0 & year == `t'
egen NN = max(tmpvar)

* Label cohort by their entry year
replace Cohort = `t' if cohortage == 0 & year == `t'

* Keep track of cohort label and initial size for each year
foreach i of numlist `t'(1)2000 {
replace EntrySizeCohort = NN if year == `i' & cohortage == `i' - `t'
replace Cohort = `t' if year == `i' & cohortage == `i' - `t'

}
drop tmpvar* NN 

}

** Calculate share of surviving firms
gen Survivalshare = N/EntrySizeCohort

** Export Data to get Figure in MATLAB
export delimited Survivalshare Cohort cohortage year using "DataInputs/Data_IndonesiaSurvival.csv", novarnames replace

restore

******************************************************************************************
****************** Figure 4: Export data on sales concentration for "Non-targeted moments"
******************************************************************************************

preserve

keep if year == 1993

keep va
sort va

gen index = _n

xtile firms_pctile =index, n(100)

collapse (count) index (sum) va, by(firms_pctile)

sort firms_pctile

drop if firms_pctile >= 97
drop if firms_pctile <= 3

gen sales_cdf = sum(va)
gen numfirms_cdf = sum(index)

egen totsales = max(sales_cdf)
egen totfirms = max(numfirms_cdf)

replace sales_cdf = sales_cdf/totsales
replace numfirms_cdf = numfirms_cdf/totfirms

export delimited sales_cdf numfirms_cdf using "DataInputs/Data_SalesDistribution.csv", novarnames replace

restore





    


