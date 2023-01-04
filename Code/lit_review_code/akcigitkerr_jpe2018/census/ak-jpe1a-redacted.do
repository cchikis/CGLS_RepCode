#delimit;
clear all; set matsize 8000; set more off;
cd [[REDACTED]];

cap n log close; log using ak-jpe1a.log, replace;

***********************************************************;
***********************************************************;
*** Patent Data Preparation                             ***;
***********************************************************;
***********************************************************;

* qstata --dofile=ak-jpe1a.do --statatype=mp --memsize=20000;
* Original programs: rd-ak1-revision4a and variants;
* Streamlined for results reported in final JPE paper;

* Data ingredients:
	p-work-pat_2008_CES: core patent data granted thru 2008
	BRDC01-MI-comb: the concordance file of patent assignees to firm identifiers
	citations: citations across patents, citations-expanded working file has assignee data already mapped to citing and cited patents
	lbdc-raw2-distr-reall1-tl: provides aggregated annual employment and establishment data for firm-year observations
	lbdc-raw2-distr-firm-sz: provides mean annual employment, industry counts, and establishment counts in 5-year blocks for firm
	lbdc-raw2-distr-firm-sic2: provides most-important sic2 code for firm across full time period in terms of employment
	lbdest: establishment-level employment data;

***********************************************************;
*** Prework on Citations Data for Self Cites            ***;
***********************************************************;

*** Prepare interim self citation file by patent;
* Citations file expanded in 3a1 to include assignees;
use if Sayear>=1977 using citations-expanded, clear;
gen cite=1; gen citeself=(Sassignee==Dassignee);
gen lag=(Sayear-Dayear); 
gen cite10=cite if lag>=0 & lag<10; 
gen citeself10=citeself if lag>=0 & lag<10;
ren Sciting patent;
collapse (sum) cite*, by(patent) fast;
compress; sort patent;
save temp-self, replace;

*** Save interim external forward citations file by patent;
use citations-expanded, clear;
gen fcitetl=1; gen fciteself=(Sassignee==Dassignee);
gen lag=(Sayear-Dayear);
gen fcitetl10=fcitetl if lag>=0 & lag<10; 
gen fciteself10=fciteself if lag>=0 & lag<10;
ren Dcited patent;
collapse (sum) fcite*, by(patent) fast;
compress; sort patent;
save temp-fwd, replace;

***********************************************************;
*** Patent Citation and Claim Prep for Own Assignees    ***;
***********************************************************;

*** Collapse 2008 data (chosen to match assignee-firmid bridge);
use patent ayear ctryn cat scat uspc invfor claims ind assignee
    if (ctryn=="US" & invfor==0 & uspc!=. & ind==1 & assignee!=.)
    using p-work-pat_2008_CES, clear; 
sort patent; drop if patent==patent[_n-1];
keep if (ayear>=1973 & ayear<=2002);
gen int yr2=.; for any 1977 1982 1987 1992 1997 2002: replace yr2=X if (ayear>=X-4 & ayear<=X);
sort patent; merge patent using temp-self, nok; erase temp-self.dta; drop _m; 
sort patent; merge patent using temp-fwd, nok; erase temp-fwd.dta; drop _m; 
for var cite* claims fcite*: replace X=0 if X==.;

*** Basic preparation - citations;
gen tct=1; 
gen self10=citeself10/cite10;
gen self10_49=(self10>=1/2 & self10!=.);
gen ecrec=fcitetl10-fciteself10;
for any 25 50 75 90: egen tempX=pctile(ecrec), p(X) by(scat yr2);
gen ec01=(ecrec<temp25 & ecrec!=.);
gen ec25=(ecrec>=temp25 & ecrec<temp50 & ecrec!=.);
gen ec50=(ecrec>=temp50 & ecrec<temp75 & ecrec!=.);
gen ec75=(ecrec>=temp75 & ecrec!=.);
gen ec90=(ecrec>=temp90 & ecrec!=.);
drop temp*; 

*** Basic preparation - claims;
for any 25 50 75 90: egen tempX=pctile(claims), p(X) by(scat yr2);
gen cl01=(claims<temp25 & claims!=.);
gen cl25=(claims>=temp25 & claims<temp50 & claims!=.);
gen cl50=(claims>=temp50 & claims<temp75 & claims!=.);
gen cl75=(claims>=temp75 & claims!=.);
gen cl90=(claims>=temp90 & claims!=.);
drop temp*;
save temp-pathold, replace;

*** Table A2;
sort assignee; compress;
joinby assignee using BRDC01-MI-comb.dta, unmatched(both);
keep if _m==3; drop _m; keep if ord==1 & qual==1;
egen temp1=min(ayear), by(firm); gen entrant=(ayear-temp1<=2); keep if temp1>=1977 & temp1<=1994; drop temp1;
egen temp1=count(ec01), by(assignee); gen wt1=1/temp1; drop temp1;
egen temp1=count(ec01), by(firm); gen wt2=1/temp1; drop temp1;
egen scatyr=group(scat yr2); xi i.scatyr;
for var ecrec claims ec01 ec25 ec50 ec75: qui areg X entrant _I* [aw=wt2], a(firm) cl(firm) \ lincom entrant;
preserve; keep if e(sample)==1; keep firm tct ecrec claims ec01 ec25 ec50 ec75 entrant; save patquality, replace; restore;

*** Collapse on assignee-yr2;
use temp-pathold, clear;
collapse (sum) tct* (mean) ec* cl* self*, by(assignee yr2) fast;

** Merge in new patent bridge; 
* Trailing r designates restricted matching on order and quality;
sort assignee; joinby assignee using BRDC01-MI-comb.dta, unmatched(both);
drop if assignee==.; keep if _m==3; drop _m;
for var tct* ec* cl* self*: gen Xr=X if ord==1 & qual==1;

*** Collapse to firm-yr2 level;
save temp1, replace; keep if tct!=. & tct!=0;
collapse (rawsum) tct
         (mean) cl25 cl50 cl75
                ecrec ec01 ec25 ec50 ec75 ec90
                self10 self10_49
		[aw=tct], by(firm yr2) fast;
sort firm yr2; save temp2, replace;
use temp1, clear; keep if tctr!=. & tctr!=0;
collapse (rawsum) tctr
         (mean) cl25r cl50r cl75r
                ecrecr ec01r ec25r ec50r ec75r ec90r
                self10r self10_49r
		[aw=tctr], by(firm yr2) fast;
sort firm yr2; merge firm yr2 using temp2;
drop _m; erase temp1.dta; erase temp2.dta;

*** Merge LBD firm and LBD-Operations SIC code;
sort firm yr2; merge firm yr2 using lbdc-raw2-distr-firm-sz; drop if _m==1; drop _m; 
sort firm; merge firm using lbdc-raw2-distr-firm-sic2; drop if _m==2; drop _m;

*** Shorten to patenting firms and save file;
gen temp1=(tct==. | tct==0); egen temp2=min(temp1), by(firm); drop if temp2==1; drop temp*;
sort firm yr2; compress; save ak-jpe1a1a, replace;

***********************************************************;
*** Prepare Highly-Cited Patent Designations            ***;
***********************************************************;

*** Create 2nd generation designations;
use temp-pathold, clear;
for var assignee ayear ec* cl*: ren X CX;
keep patent Cayear C*; ren patent citation;
sort citation; save temp1s-cited.dta, replace;
erase temp-pathold.dta;

*** Collapse 2008 data (chosen to match assignee-firmid bridge);
use patent ayear ctryn cat scat uspc invfor claims ind assignee
    if (ctryn=="US" & invfor==0 & uspc!=. & ind==1 & assignee!=.)
    using p-work-pat_2008_CES, clear; 
sort patent; drop if patent==patent[_n-1];
keep if (ayear>=1973 & ayear<=2002);
gen int yr2=.; for any 1977 1982 1987 1992 1997 2002: replace yr2=X if (ayear>=X-4 & ayear<=X);
keep patent assignee ayear yr2; 
sort patent; drop if patent==patent[_n-1]; 
sort patent; save temp1s-citing, replace;

*** Load and merge citations data;
use citations, clear;
sort patent; merge patent using temp1s-citing;
keep if _m==3; drop _m; erase temp1s-citing.dta;
sort citation; merge citation using temp1s-cited, nok; 
keep if _m==3; drop _m; erase temp1s-cited.dta;
drop if assignee==Cassignee; drop Cassignee;

*** Basic preparation and collapse;
gen age=ayear-Cayear; drop if age<0; drop Cayear;
gen Cct=1; for var C*: gen X_5=X if age<=5;
collapse (sum) Cct* (mean) Cec* Ccl*, by(assignee yr2) fast;

** Merge in new patent bridge; 
sort assignee; joinby assignee using BRDC01-MI-comb.dta, unmatched(both);
drop if assignee==.; keep if _m==3; drop _m;
for var C*: gen Xr=X if ord==1 & qual==1;

*** Save output file;
collapse (rawsum) Cct* (mean) Cec* Ccl*, by(firm yr2) fast;
sort firm yr2; compress;
save ak-jpe1a1b, replace;

***********************************************************;
*** Unite Data and Build Growth Variables               ***;
***********************************************************;

*** Prepare innovative sample;
use if yr2>1977 & yr2<2002 using ak-jpe1a1a, clear;
sort firm yr2; merge firm yr2 using ak-jpe1a1b, nok; drop _m; compress;
for any a b: erase ak-jpe1a1X.dta;

*** Group preps;
egen fgrp=group(firm); 
egen xyr2=group(xsic yr2);
for var emp tct*: gen lX=ln(X); 

*** Prepare employment growth variables;
sort fgrp yr2;  
gen empf1=emp[_n+1] if fgrp==fgrp[_n+1] & yr2==yr2[_n+1]-5;
gen empshgr=(empf1-emp)/emp;
gen empshgr2=empshgr; replace empshgr2=-1 if empshgr2==. & yr2!=1997;
for any empshgr empshgr2: replace X=. if yr2==1997;
replace empshgr=10 if empshgr>10 & empshgr!=.;
replace empshgr2=10 if empshgr2>10 & empshgr2!=.;

*** Narrow to restricted sample;
gen temp1=(tctr==. | tctr==0); egen temp2=min(temp1), by(firm); drop if temp2==1; drop temp*;

*** Scaling panel definitions;
* min year still allows zero values;
egen minyr2=min(yr2), by(fgrp);

*** Patent intensity prep;
gen tctremp=tctr/emp if empshgr2!=.;
for var tctremp: egen temp1=mean(X), by(yr2) \ egen temp2=sd(X), by(yr2) \ gen Z1X=(X-temp1)/temp2 \ drop temp*;

*** Define samples;
gen patgrp=(tctr!=0 & tctr!=.) & empshgr2!=.;
gen temp1=tctr; replace temp1=0 if temp1==.; replace temp1=. if yr2==1997;
egen temp2=min(temp1), by(firm); gen contfirm=(temp2>0 & temp2!=.); drop temp*;

*** Save working file;
compress; save ak-jpe1a-core, replace;

***********************************************************;
***********************************************************;
*** ANALYSIS USING FIVE-YEAR INTERVALS                  ***;
***********************************************************;
***********************************************************;

*** Table 5 - entry share;
use ak-jpe1a-core, clear;
keep if yr>=1982 & yr<=1992; gen ct=1;
gen entry=(yr2==minyr2); keep if yr2>=1987;
egen temp1=sum(emp), by(entry); tab entry, s(temp1);
gen rate=[[REDACTED]]; sum rate;

*** Fact 2.4-2.5 and Table 1;
use ak-jpe1a-core, clear;
areg Z1tctremp lemp, a(xyr2) cl(fgrp);
keep if contfirm==1;
areg Z1tctremp lemp, a(xyr2) cl(fgrp);
predict temp1, resid; drop if temp1==.; drop temp1;
for var Z1tctremp ec90r ec01r ec25r ec50r ec75r: areg X lemp, a(xyr2) cl(fgrp);
xtile ctile=emp, nq(20); for any emp tctremp: tab ctile, s(X); drop ctile;
xtile ctile=est, nq(20); for any est tctremp: tab ctile, s(X); drop ctile;
keep if e(sample)==1; keep firm emp est tctremp Z1tctremp ec90r ec01r ec25r ec50r ec75r; save contfirm1, replace;

*** Table 8-9;
use if contfirm==1 & yr2!=1997 using ak-jpe1a-core, clear;
gen patemp=tctr/emp; egen temp1=std(patemp); replace patemp=temp1; drop temp1;
gen sszbin=0; for num 20 40 60 80: egen tempX=pctile(lemp), p(X) \ replace sszbin=X if lemp>tempX \ drop tempX;
xi i.sszbin;
for var empshgr2 patemp self10_49r ec90r ec01r ec25r ec50r ec75r: areg X _Issz*, a(xyr2) cl(fgrp);

*** Table 10 (cut at 20%);
gen selfL=(self10r>0 & self10r<=.2);
gen selfH=(self10r>0.2 & self10r<=1);
for any empshgr2:
\ areg X lemp ltctr ec50r ec75r selfL selfH, a(xyr2) cl(fgrp)
\ areg X lemp ltctr cl50r cl75r selfL selfH, a(xyr2) cl(fgrp);
keep if e(sample)==1; keep firm empshgr2 emp ltctr ec50r ec75r selfL selfH _Issz*; save contfirm2, replace;

***********************************************************;
***********************************************************;
*** ANALYSIS USING SINGLE-YEAR INTERVALS                ***;
***********************************************************;
***********************************************************;

*** Develop annual panel of our firms;
use if patgrp==1 using ak-jpe1a-core, clear;
gen yrmin=yr2; gen yrmax=yr2; gen empmin=emp; gen empmax=emp; 
collapse (min) yrmin empmin (max) yrmax empmax (mean) xsic2, by(firm) fast;
for var empmin empmax: gen temp1=(X>500) \ tab temp1 \ drop temp1 X;
sort firm; save rd-ak1-firmlist, replace;
use lbdc-raw2-distr-reall1-tl.dta, clear;
gen temp1=emp>500 & emp!=.; sum temp1; drop temp1; 
sort firm; merge firm using rd-ak1-firmlist; keep if _m==3; drop _m;

*** Calculate growth properties;
sort firm yr;
gen empf1=emp[_n+1] if firm==firm[_n+1] & yr==yr[_n+1]-1;
gen empshgr2a=(empf1-emp)/emp; 
replace empshgr2a=-1 if empshgr2a==. & yr<2001;
replace empshgr2a=10 if empshgr2a>10 & empshgr2a!=.;
replace empshgr2a=. if yr<=1978 | yr>=2001;
gen lemp=ln(emp); 

*** Fact 2.3 Gibrat's Law;
egen fgrp=group(firm); egen xyr=group(xsic2 yr);
keep if empshgr2a!=.;
areg empshgr2a lemp if (yr>=yrmin-2 & yr<=yrmax+2), a(xyr) cl(fgrp);
predict temp1, resid; drop if temp1==.; drop temp1;
sum empshgr2a, mean; display r(mean);
egen fct=group(firm); egen firmct=max(fct); sum firmct, mean; display r(mean); drop firmct fct;
for num 25 50 75: egen tempX=pctile(emp), p(X) \ sum emp if emp>=tempX*0.90 & emp<=tempX*1.10, mean \ display r(mean) \ drop tempX;
xtile ctile=emp, nq(20); for any emp empshgr2a: tab ctile, s(X); drop ctile;
xtile ctile=est, nq(20); for any est empshgr2a: tab ctile, s(X); drop ctile;
keep if e(sample)==1; keep fgrp emp est empshgr2a lemp; save gibrat1, replace;

*** Extension beyond core sample;
use firm yr emp bestsic if yr>=1987 & yr<=1997 using lbdest.dta, clear;
egen fgrp=group(firm); drop firm;
gen sic4=substr(bestsic,1,4); destring sic4, gen(sic) force; drop bestsic sic4;
gen empM=emp if sic>=2000 & sic<4000; ren emp empT;
collapse (sum) emp*, by(fgrp yr) fast; compress;
for any M:
\ sort fgrp yr 
\ gen empXshgr=(empX[_n+1]-empX)/empX
\ replace empXshgr=-1 if empXshgr==. & yr!=1997 & empX!=.
\ replace empXshgr=10 if empXshgr>10 & empXshgr!=. & empX!=.
\ gen lempX=ln(empX)
\ areg empXshgr lempX, a(yr) cl(fgrp);
keep if e(sample)==1; keep fgrp empM empMshgr lemp; save gibrat2, replace;

*** End of Program;
cap n log close;

***********************************************************;
*** Saved Materials                                     ***;
***********************************************************;