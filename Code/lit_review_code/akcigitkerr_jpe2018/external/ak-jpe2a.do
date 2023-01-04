#delimit;
cd /export/projects/wkerr_h1b_project/kerr/mainwork-2009/programs/akcigit/replication;
cap n log close; log using ak-jpe2a.log, replace; 

* William Kerr;
* To God's Glory;
* Provides external patent calculations;
* stata-mp8-10g -b do ak-jpe2a.do;

clear all; set matsize 5000; set more off;

************************;
** CITATION DISTR    ***;
************************;

* Data folder contains file of US industrial patents (number, assignee, app year)
  and also citations file (including non-US and non-industrial) for patents granted thru 2009;

for any citations-2009 patents-2009 patents-nber ak-cite-compare1: cap n !gunzip ./data/X.dta.gz;

*** Merge in data for citing patent;
use ./data/citations-2009, clear;
sort patent; merge patent using ./data/patents-2009, nok; keep if _m==3; drop _m; 
for var ayear assignee: ren X SX; ren patent Sciting;

*** Merge in data for cited patent;
ren citation patent; sort patent; merge patent using ./data/patents-2009, nok; keep if _m==3; drop _m; 
for var ayear assignee: ren X DX; ren patent Dcited;
gen lag=Sayear-Dayear; keep if lag>=0 & lag<=10; drop lag;

*** Identify external patents;
preserve;
gen ct=1; gen ect=Sassignee!=Dassignee;
collapse (sum) ct ect, by(Sciting) fast; replace ect=0 if ect==.;
gen temp1=(ct-ect)/ct; keep if temp1>=0.5;
keep Sciting; ren Sciting patent; sort patent; save temp1-int, replace;
restore;

*** Collapse citation counts;
gen ct=1; gen ect=Sassignee!=Dassignee;
collapse (sum) ct ect (mean) Dayear, by(Dcited) fast;

*** Merge into patents;
ren Dcited patent; sort patent; merge patent using ./data/patents-2009, keep(ayear); 
replace Dayear=ayear if _m==2; drop _m;
for var ct ect: replace X=0 if X==.;
sort patent; merge patent using temp1-int, nok; erase temp1-int.dta; 
gen ext=(_m!=3); drop _m; 

*** Table 5 Moments and 1985-1994 distribution;
keep if Dayear>=1985 & Dayear<1994;
tab ext; tab ext, s(ect); keep if ext==1; gen obs=1;
collapse (sum) obs, by(ect);
compress; egen temp1=sum(obs); gen obs_per=obs/temp1; drop temp1; 
format obs %8.0f; format obs_per %6.5f;
sort ect; list;

************************;
** MISC PIECES       ***;
************************;

*** Table A3;
use ./data/ak-cite-compare1.dta, clear;
egen Zsyr=group(Zscat Zyr2); xi i.Zsyr;
for var Zclaims Zec*: qui areg X external _I*, a(cited) cl(cited) \ lincom external;

*** Figure 3 (uses original NBER data to match 2010 NBER working paper);
* set is US industrial patents with no foreign inventors applied for 1975-1984;
use ./data/patents-nber, clear;
gen nocite=crec==0; replace nocite=1 if secd==1; 
replace crec=crec*(1-secd); replace crec=0 if crec==.; replace crec=int(crec);
gen ct=1; gen ctS0=self==0; gen ctS1=self>0.00 & self<=0.50; gen ctS2=self>0.50 & self<=1.00;
replace crec=50 if crec>50 & crec!=.;
collapse (sum) ct ctS0 ctS1 ctS2, by(crec) fast;
for var ct*: egen temp1=sum(X) \ gen pX=X/temp1 \ drop temp1;
format ct* %5.0f; format p* %4.3f;
sort crec; list, noobs clean;

for any citations-2009 patents-2009 patents-nber ak-cite-compare1: cap n !gzip ./data/X.dta;

*** End of program;
log close;