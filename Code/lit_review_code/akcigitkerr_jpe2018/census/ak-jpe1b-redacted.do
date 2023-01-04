#delimit;
clear all; set matsize 8000; set more off;
cd [[REDACTED]];

cap n log close; log using ak-jpe1b.log, replace;

***********************************************************;
***********************************************************;
*** ADDTIONAL SURVEY-BASED MATERIAL                     ***;
***********************************************************;
***********************************************************;

* qstata --dofile=ak-jpe1b.do --statatype=mp --memsize=20000;
* Original programs: rd-ak1-fig1a/b and variants;
* Streamlined for results reported in final JPE paper;

* Data ingredients:
	brdis2008: 2008 BRDIS survey
	rd-ak1-prod1.dta: appends together the NSF R&D surveys
	rd-ak1-prod1-lbd.dta: LBD employment in three-year windows to merge to the NSF R&D surveys
	nber_ssel_bridge-merge-comb: earlier patent bridge to match NSF R&D survey;

* qstata --dofile=ak-jpe1b.do --statatype=mp --memsize=20000;

************************************************************;
** 2008 BRDIS CORRELATIONS                                **;
************************************************************;

*** 2008 RD data;
* Partially weakens with zero-inferred shares;
use id wgt smpwgt emp_r wbusp_r ybusp_r wbthp_r ybthp_r wnewp_r ynewp_r
    using brdis2008, clear;
gen lemp=ln(emp);
pwcorr emp* wbusp* ybusp* wbthp* ybthp* wnewp* ynewp*, star(0.05); 
pwcorr lemp* wbusp* ybusp* wbthp* ybthp* wnewp* ynewp*, star(0.05); 
pwcorr lemp* wbusp* ybusp* wbthp* ybthp* wnewp* ynewp* [aw=smpwgt], star(0.05); 
pwcorr lemp* wbusp* ybusp* wbthp* ybthp* wnewp* ynewp* [aw=wgt], star(0.05); 
ren id firm; keep firm emp lemp ybthp_r; save brdis, replace;

************************************************************;
** NSF PRODUCT-PROCESS                                    **;
************************************************************;
* Abbreviated analsis from original NBER working paper;

*** Patent prep work (matches R&D survey, origin: rd-full-prep4);
use nber_ssel_bridge-merge-comb, clear;
gen yr2=.; for any 1979 1981 1985 1987 1989 1991: replace yr2=X if (yr>=X-1 & yr<=X+1);
drop if yr2==.; drop yr; ren yr2 yr;
collapse (sum) dct_self00s dct_self01s dct_self25s, by(firm yr) fast; 
drop if (dct_self00s==0 | dct_self00s==.); drop if yr==1991; compress;
sort firm yr; save nber_ssel_bridge-merge-brdc1-firm-temp, replace;

*** Summarize restricted dataset;
use id yr wt rdtot rdfed rdcomp proc* prod* dne dns sic using rd-ak1-prod1.dta, clear;
gen temp1=1 if prodtot>0 & prodtot!=.; egen temp2=sum(temp1), by(yr); drop if temp2==0; drop temp*;
drop if (rdtot==0 | rdtot==.); drop if yr==1991;

*** Compare and clean-up firm records;
replace id=substr(id,1,6)+"0000" if substr(id,1,1)!="0"; egen records=count(rdtot), by(id); 
egen temp1=sum(rdtot), by(id yr); egen temp2=max(rdtot), by(id yr); gen temp3=temp1-temp2;
drop if records>5 & temp3==temp2 & rdtot!=temp2; 
collapse (sum) rdtot proctot prodtot, by(id yr) fast;

*** Merge in LBD employment data and patent data - longitudinal;
sort id; merge id using rd-ak1-prod1-lbd; drop if _m==2; drop _m;
ren id firm; sort firm yr; merge firm yr using nber_ssel_bridge-merge-brdc1-firm-temp, nok;
drop _m; ren firm id; erase nber_ssel_bridge-merge-brdc1-firm-temp.dta;

*** Prepare LBD measures;
for num 1/3: gen empX=.;
replace emp1=emp78 if yr==1979; replace emp2=emp79 if yr==1979; replace emp3=emp80 if yr==1979;
replace emp1=emp80 if yr==1981; replace emp2=emp81 if yr==1981; replace emp3=emp82 if yr==1981;
replace emp1=emp84 if yr==1985; replace emp2=emp85 if yr==1985; replace emp3=emp86 if yr==1985;
replace emp1=emp86 if yr==1987; replace emp2=emp87 if yr==1987; replace emp3=emp88 if yr==1987;
replace emp1=emp88 if yr==1989; replace emp2=emp89 if yr==1989; replace emp3=emp90 if yr==1989;
replace emp1=emp90 if yr==1991; replace emp2=emp91 if yr==1991; replace emp3=emp92 if yr==1991;
egen empav=rmean(emp1 emp2 emp3);

*** Prepare metrics;
gen ratio=prodtot/(prodtot+proctot); 
replace ratio=1 if prodtot>0 & prodtot!=. & proctot==.;
replace ratio=0 if proctot>0 & proctot!=. & prodtot==.;
gen ratioz=1-ratio; gen doproc=(ratio<1 & ratio!=.);
for var empav: gen lX=ln(X); keep if empav!=.; 
for any 01: gen patratX=dct_selfXs/dct_self00s;

*** Correlations;
pwcorr lempav ratioz doproc patrat01, star(0.05);
keep if lempav!=.; ren id firm; ren empav emp;
keep firm emp lempav ratioz doproc patrat01; save rdprod, replace;

***********************************************************;
*** Patent & R&D Intensity                              ***;
***********************************************************;

/* Conversion background: cd [[REDACTED]];
use rd-2ind4-prep2a2, clear; 
sort firm yr; merge firm yr using rd-2ind4-prep2a3, nok; tab _m yr; drop _m;
sort firm yr; merge firm yr using rd-2ind-prep3b1a, nok; tab _m yr; drop _m;
keep if rdtot>0 & CNF_tvs>0 & CNF_emp>0; drop if rdtot==. | zdcts==. | CNF_tvs==. | CNF_emp==.; 
sum; tab yr; keep if yr>=1982 & yr<=1992;
gen z=1; collapse (sum) rdtot zdcts CNF_tvs CNF_emp LBF_emp, by(z) fast;
gen rdsal=rdtot/CNF_tvs; gen patemp=zdcts/CNF_emp; gen conv=rdsal/patemp; sum; */
use if yr>=1982 & yr<=1992 using ak-jpe1a-core, clear;
keep tctr emp firm; save patint, replace;
gen z=1; collapse (sum) tctr emp, by(z) fast; gen patemp=tctr/emp; gen rdsalconv=patemp*[[REDACTED]]; sum rdsalconv, mean; display r(mean);

*** End of Program;
cap n log close;

***********************************************************;
*** Saved Materials                                     ***;
***********************************************************;