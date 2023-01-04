#delimit;
cap n log close; 
cd C:\kerr\docs-ces\akcigit\growth1\data-jpe\model;
log using AK-comparison.log, replace; 

* William Kerr;
* To God's Glory;

clear all; set matsize 5000; set more off;

**************************************************************;
**************************************************************;
**************************************************************;

use dataGrowth15, clear;
sort id; drop if id==id[_n-1]; sort id; save temp, replace;
use dataInnovation15, clear;
sort id; merge id using temp; tab _m; keep if _m==3; drop _m; erase temp.dta;
gen szbin=0;
for num 20 40 60 80: egen tempX=pctile(firmSize), p(X) \ replace szbin=X if firmSize>tempX \ drop tempX;
drop if totalPatent==0;
for var radical internal patentq*: replace X=X/totalPatent \ replace X=0 if X==.;
gen selfL=(internal>0 & internal<=.2);
gen selfH=(internal>0.2 & internal!=.);
gen patemp=totalPatent/firmSize; egen temp2=std(patemp); replace patemp=temp2; drop temp2;
xi i.szbin;
for var growth patemp internal radical patentq1 patentq2 patentq3 patentq4: regress X _I*, r;
for var growth: regress X firmSize totalPatent patentq3 patentq4 selfL selfH, r;

*** End of program;
log close;