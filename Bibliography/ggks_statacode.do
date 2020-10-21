*************************************************************************************************************************
* STATA CODE FOR:													*
* "Marry Your Like: Assortative Mating and Income Inequality"					                        *
* by Jeremy Greenwood, Nezih Guner, Georgi Kocharkov and Cezar Santos                                                   *
* published in the American Economic Review Papers & Proceedings, May 2014		                                *
*                                                                                                                       *
* Description of performed tasks:											*
* 1. This code cleans the IPUMS sample (for a description of the data, see the online appendix). 			*
*    Note: the data for each individual year should be downloaded as a dta file called sampleYEAR.dta for 		*
*    YEAR=1960,1970,1980,1990,2000,2005											*
* 2. The code computes the Kendall's tau statistics that measure the correlation of education between spouses		*
*    over the years.													*
* 3. The code runs the regression between years of schooling for husband vs wife (plus controls)			*
* 4. The data is used to derive  the fractions of each household type out of the total population of households.        *
*    The average income relative to mean household income is derived for each household type as well.                   *
*    These statistics are recorded and are used to compute inequality statistics and perform counterfactual experiments *
*    in the MATLAB environment (see ggks_matlabcode.m)                                                                  *
* 5. Some more statistics are computed:                                                                                 *
*    (i)  Share of wife's labor income out of total household labor income in a contingency table                       *
*    (ii) Wife's share income and LFP by centiles	                                                                *
*    (iii) Contingency tables					                                                        *
*************************************************************************************************************************

* Preliminaries
clear
set more off
pause on

*Years used in the analysis
local year_codes 1960 1970 1980 1990 2000 2005
local years_ex60 1970 1980 1990 2000 2005

*Stata directory
local stata_dir `c(pwd)'

*Change to input directory
cd "Input"

*Cleaning the Data
******************

local counter_tloop = 0 /*loop counter*/

foreach t of local year_codes {

    use sample`t', clear

	*Update counter
	local counter_tloop = `counter_tloop' + 1 
        
	*Define the universe of households
	keep if age >= 25 & age <= 54 /*age restriction*/
	keep if eldch <= 19 | eldch == 99 /*keep household with children younger than 19 or N/A children*/
	keep if yngch <= 19 | yngch == 99 /*keep household with children younger than 19 or N/A children*/
	drop if marst == 2 /*2 is spouse missing, so we cannot match these families with sploc and these will automatically have to go*/
	drop if marst == 5 /*5 is widowed, shall we drop them or treat them as singles?*/
	drop if empstat == 0 /*employment status N/A*/

        *Income levels restrictions 
	drop if ftotinc < 0 /*remove negative incomes*/
    	drop if ftotinc >= 9999999 /*remove N/A values*/
        			
	*Education, marital status, fertility, sex, participation
	gen lhs = educ<6 /*less than highschool*/
	gen hs = educ==6 /*highschool*/
	gen smcol = (educ>6&educ<10) /*some college*/ 	
	gen col = educ==10 /*college*/
	gen mcol = educ==11 /*college  plus*/ 
	
	gen mar = marst == 1 /*2 is spouse missing, so we cannot match these families with sploc and these will automatically have to go*/
	 
	gen nmar = marst == 6 /*never married*/
	
	gen div = (marst == 4 | marst == 3)  /*separated included*/
	
	gen fem = sex == 2

	gen mal = sex == 1
	
	gen kd0 = nchild==0 /* 0 kids*/
	gen kd1 = nchild==1 /* 1 kids*/
	gen kd2 = nchild==2 /* 2 kids*/
	gen kdm2 = nchild>2 /* more than 2 kids*/

	gen wrk = empstat == 1 /*working*/
	gen nwrk = 1 - wrk /*not working*/
	
	*Characteristics of the marital partner
	bysort serial:gen peduc =  educ[sploc]
	bysort serial:gen pempstat =  empstat[sploc]
	bysort serial:gen pincwage = incwage[sploc]

	/*drop if characteristics are missing for married people*/
	drop if mar == 1 & (peduc == . | pempstat ==. | pincwage == .)

	/*check once again whether all married have their spouses with characteristics*/
	gen plhs = peduc<6 /*less than highschool*/
	gen phs = peduc==6 /*highschool*/
	gen psmcol = (peduc>6 & peduc<10) /*some college*/ 	
	gen pcol = peduc==10 /*college*/
	gen pmcol = peduc==11 /*college  plus*/ 
	
	gen pwrk = pempstat == 1
	gen pnwrk = 1 - pwrk

	/*drop married females cause we are going to work only with the married males*/
	drop if mar==1 & fem==1
	
	// Keep only married people:
	keep if mar==1

	/*more restrictions: remove households with subfamilies*/
	gen xfamsize = famsize - nchild
	drop if mar == 1 & xfamsize > 2 /*only the two parents should live with their children*/
	drop if nmar == 1 & xfamsize > 1 /*only the single parent should live with her children*/
	drop if div == 1 & xfamsize > 1 /*only the divorced parent should live with her children*/

	*Equivalence scale
	/*OECD equivalence scale: first adult 1, second adult 0.7, each child 0.5
	  OECD-modified scale:    first adult 1, second adult 0.5, each child 0.3
	  Square root scale:      income / sqrt(family size)                   */

	*Family income: OECD equivalence scales
	gen ftotinc_adj = .
	replace ftotinc_adj = ftotinc / (1 + 0.5 + 0.3 * nchild) if mar == 1 /*married*/
 	replace ftotinc_adj = ftotinc / (1 + 0.3 * nchild) if (nmar == 1 | div == 1) /*single or divorced*/

	*Labor income
	gen labinc =.
	replace labinc = incwage if nmar == 1 | div == 1
	replace labinc = incwage + pincwage if mar == 1

	*Labor income: OECD equivalence scales
	gen labinc_adj =.
	replace labinc_adj = labinc / (1 + 0.5 + 0.3 * nchild) if mar == 1 /*married*/
 	replace labinc_adj = labinc / (1 + 0.3 * nchild) if (nmar == 1 | div == 1) /*single or divorced*/
	
	*Generate categorical variable for education:
	gen educ_cat = 0
	replace educ_cat = 1 if lhs==1
	replace educ_cat = 2 if hs==1
	replace educ_cat = 3 if smcol==1
	replace educ_cat = 4 if col==1
	replace educ_cat = 5 if mcol==1
	
	gen peduc_cat = 0
	replace peduc_cat = 1 if plhs==1
	replace peduc_cat = 2 if phs==1
	replace peduc_cat = 3 if psmcol==1
	replace peduc_cat = 4 if pcol==1
	replace peduc_cat = 5 if pmcol==1
	
	saveold sample`t'_cleaned, replace

}

*Merging samples for all years into one data file
use sample1960_cleaned, clear
foreach t of local years_ex60 {

	append using sample`t'_cleaned

}
saveold sampleall_cleaned, replace

*Computing correlation (Kendall's tau)
**************************************

*Change to output directory
cd "`stata_dir'"
cd "Output"

*Record the Kendall's tau in the file kendall_tau.txt
file open kendall_tau using "kendall_tau.txt", write replace

foreach t of local year_codes {

        cd "`stata_dir'"
        cd "Input"   

        use sample`t'_cleaned, clear
	display "year = " `t'

	*The computatation of Kendall tau takes a while
	*You can comment out the next 4 lines to save time

	ktau educ_cat peduc_cat, stats(taua taub p)

	scalar kendall_tau_a = r(tau_a)

        cd "`stata_dir'"
        cd "Output"

	file write kendall_tau %16.15f (kendall_tau_a) _n
}

file close kendall_tau


*Running regression on educational attainment of spouses
********************************************************

*Change to input directory
cd "`stata_dir'"
cd "Input"

*Load data
use sampleall_cleaned, clear

*Generate variable with years of education
gen educ_years = 0
replace educ_years = 4  if educ== 1
replace educ_years = 8  if educ== 2
replace educ_years = 9  if educ== 3
replace educ_years = 10 if educ== 4
replace educ_years = 11 if educ== 5
replace educ_years = 12 if educ== 6
replace educ_years = 13 if educ== 7
replace educ_years = 14 if educ== 8
replace educ_years = 15 if educ== 9
replace educ_years = 16 if educ== 10
replace educ_years = 17 if educ== 11

gen peduc_years = 0
replace peduc_years = 4  if peduc== 1
replace peduc_years = 8  if peduc== 2
replace peduc_years = 9  if peduc== 3
replace peduc_years = 10 if peduc== 4
replace peduc_years = 11 if peduc== 5
replace peduc_years = 12 if peduc== 6
replace peduc_years = 13 if peduc== 7
replace peduc_years = 14 if peduc== 8
replace peduc_years = 15 if peduc== 9
replace peduc_years = 16 if peduc== 10
replace peduc_years = 17 if peduc== 11

*Generate year dummies and interaction with educ_years variables
foreach t of local years_ex60 {
	gen d_`t' = 0
	replace d_`t' = 1 if year==`t'
	gen educ_years_`t' = 0
	replace educ_years_`t' = educ_years if year==`t'
}


*Run regressions - all years merged
reg peduc_years educ_years d_1970 d_1980 d_1990 d_2000 d_2005 /*
    */ educ_years_1970 educ_years_1980 educ_years_1990 educ_years_2000 educ_years_2005

*Save the results
mat all_coeffs = e(b)

*Change to output directory
cd "`stata_dir'"
cd "Output"

file open gamma_coeff using "gamma_coeff.txt", write replace

file write gamma_coeff %16.15f (all_coeffs[1,7]) _n
file write gamma_coeff %16.15f (all_coeffs[1,8]) _n
file write gamma_coeff %16.15f (all_coeffs[1,9]) _n
file write gamma_coeff %16.15f (all_coeffs[1,10]) _n
file write gamma_coeff %16.15f (all_coeffs[1,11]) _n

file close gamma_coeff

*Change to input directory
cd "`stata_dir'"
cd "Input"

*Erase the temporary files
erase sample1960_cleaned.dta
erase sample1970_cleaned.dta
erase sample1980_cleaned.dta
erase sample1990_cleaned.dta
erase sample2000_cleaned.dta
erase sample2005_cleaned.dta
erase sampleall_cleaned.dta


*Slicing the data into different sub-groups
*******************************************
local year_codes 1960 2005

*Main loop (years)

local counter_tloop = 0 /*loop counter*/

foreach t of local year_codes {

	*Change to input directory
	cd "`stata_dir'"
	cd "Input"

        use sample`t', clear

	*Update counter
	local counter_tloop = `counter_tloop' + 1
        
	*Define the universe of households
	keep if age >= 25 & age <= 54 /*age restriction*/
	keep if eldch <= 19 | eldch == 99 /*keep household with children younger than 19 or N/A children*/
	keep if yngch <= 19 | yngch == 99 /*keep household with children younger than 19 or N/A children*/
	drop if marst == 2 /*2 is spouse missing, so we cannot match these families with sploc and these will automatically have to go*/
	drop if marst == 5 /*5 is widowed, shall we drop them or treat them as singles?*/
	drop if empstat == 0 /*employment status N/A*/

        *Income levels restrictions 
	drop if ftotinc < 0 /*remove negative incomes*/
        drop if ftotinc >= 9999999 /*remove N/A values*/
	
	*Education, marital status, fertility, sex, participation
	gen lhs = educ<6 /*less than highschool*/
	gen hs = educ==6 /*highschool*/
	gen smcol = (educ>6&educ<10) /*some college*/ 	
	gen col = educ==10 /*college*/
	gen mcol = educ==11 /*college  plus*/ 
	
	gen mar = marst == 1 /*2 is spouse missing, so we cannot match these families with sploc and these will automatically have to go*/
	 
	gen nmar = marst == 6 /*never married*/
	
	gen div = (marst == 4 | marst == 3)  /*separated included*/
	
	gen fem = sex == 2

	gen mal = sex == 1
	
	gen kd0 = nchild==0 /* 0 kids*/
	gen kd1 = nchild==1 /* 1 kids*/
	gen kd2 = nchild==2 /* 2 kids*/
	gen kdm2 = nchild>2 /* more than 2 kids*/

	gen wrk = empstat == 1 /*working*/
	gen nwrk = 1 - wrk /*not working*/
	
	*Characteristics of the marital partner
	bysort serial:gen peduc =  educ[sploc]
	bysort serial:gen pempstat =  empstat[sploc]
	bysort serial:gen pincwage = incwage[sploc]

	/*drop if characteristics are missing for married people*/
	drop if mar == 1 & (peduc == . | pempstat ==. | pincwage == .)

	/*check once again whether all married have their spouses with characteristics*/
	gen plhs = peduc<6 /*less than highschool*/
	gen phs = peduc==6 /*highschool*/
	gen psmcol = (peduc>6 & peduc<10) /*some college*/ 	
	gen pcol = peduc==10 /*college*/
	gen pmcol = peduc==11 /*college  plus*/ 
	
	gen pwrk = pempstat == 1
	gen pnwrk = 1 - pwrk

	/*drop married females cause we are going to work only with the married males*/
	drop if mar==1 & fem==1

	/*more restrictions: remove households with subfamilies*/
	gen xfamsize = famsize - nchild
	drop if mar == 1 & xfamsize > 2 /*only the two parents should live with their children*/
	drop if nmar == 1 & xfamsize > 1 /*only the single parent should live with her children*/
	drop if div == 1 & xfamsize > 1 /*only the divorced parent should live with her children*/

	*Equivalence scale
	/*OECD equivalence scale: first adult 1, second adult 0.7, each child 0.5
	  OECD-modified scale:    first adult 1, second adult 0.5, each child 0.3
	  Square root scale:      income / sqrt(family size)                   */

        *Family income: OECD equivalence scales
	gen ftotinc_adj = .
	replace ftotinc_adj = ftotinc / (1 + 0.5 + 0.3 * nchild) if mar == 1 /*married*/
 	replace ftotinc_adj = ftotinc / (1 + 0.3 * nchild) if (nmar == 1 | div == 1) /*single or divorced*/

	*Labor income
	gen labinc =.
	replace labinc = incwage if nmar == 1 | div == 1
	replace labinc = incwage + pincwage if mar == 1

	*Labor income: OECD equivalence scales
	gen labinc_adj =.
	replace labinc_adj = labinc / (1 + 0.5 + 0.3 * nchild) if mar == 1 /*married*/
 	replace labinc_adj = labinc / (1 + 0.3 * nchild) if (nmar == 1 | div == 1) /*single or divorced*/
	
	*Basic types of variables
	gen nmar_m = (nmar==1 & mal==1)
	gen nmar_f = (nmar==1 & fem==1)
	gen div_m = (div==1 & mal==1)
	gen div_f = (div==1 & fem==1)
		
	*Define local for different demographic types of households 
	local marstat mar nmar_m nmar_f div_m div_f
	local education lhs hs smcol col mcol
	local market wrk nwrk
	local kids kd0 kd1 kd2 kdm2
	
	local peducation plhs phs psmcol pcol pmcol
	local pmarket pwrk pnwrk

        *Loop over all variables for income
        foreach var of varlist ftotinc ftotinc_adj {

        *Total income of all households
        quietly sum `var' [aweight=perwt]
	scalar total_`var'_`t' = r(sum)
	scalar mean_`var'_`t' = r(mean)

        *Define centiles of the income distributions and create dummies for the centiles
        _pctile `var' [aweight=perwt], nquantiles(10)
        
	scalar p10_`var' = r(r1) /*thresholds for deciles*/ 
        scalar p20_`var' = r(r2)
	scalar p30_`var' = r(r3)
	scalar p40_`var' = r(r4)
	scalar p50_`var' = r(r5)
	scalar p60_`var' = r(r6)
	scalar p70_`var' = r(r7)
	scalar p80_`var' = r(r8)
	scalar p90_`var' = r(r9)

	gen c1_`var'  =  (`var'<= p10_`var')  /*dummies for deciles*/
	gen c2_`var'  =  (`var' > p10_`var' & `var' <= p20_`var')
	gen c3_`var'  =  (`var' > p20_`var' & `var' <= p30_`var')
	gen c4_`var'  =  (`var' > p30_`var' & `var' <= p40_`var')
	gen c5_`var'  =  (`var' > p40_`var' & `var' <= p50_`var')
	gen c6_`var'  =  (`var' > p50_`var' & `var' <= p60_`var')
	gen c7_`var'  =  (`var' > p60_`var' & `var' <= p70_`var')
	gen c8_`var'  =  (`var' > p70_`var' & `var' <= p80_`var')
	gen c9_`var'  =  (`var' > p80_`var' & `var' <= p90_`var')
	gen c10_`var' =  (`var' > p90_`var')

        }

	local centiles 1 2 3 4 5 6 7 8 9 10

	*Create the local with types of variables

	local i = 0 /*loop counter*/
	
	*Create a numerical variable for the type
	gen household_types_num_`t' = 0

        *Change to output directory
        cd "`stata_dir'"
        cd "Output"

	*Open a file in which we write some results: household types
	if `t' == 1960 {
	
	file open types using "types.txt", write replace
	
	}

	*Open file in which to write the fraction of an (i,j)-household out of the total population of households
        *Open file in which to write the share of total income owned by an average (i,j)-household
        foreach var of varlist ftotinc_adj {

	file open f_`var'_`t' using "f_`var'_`t'.txt", write replace
	file open r_`var'_`t' using "r_`var'_`t'.txt", write replace

        }
	
        *Loop over types of households (i)	
	foreach i1 of local marstat {

		foreach i2 of local education {

			foreach i3 of local market {

				foreach i4 of local kids {

                                        if `i' > 399 { /*this refers to all single and divorced household types*/

					*Update counter 
					local i = `i' + 1
					di "`i'"

					*Calculate the fraction of type-i household in income decile j: type-(i,j), sum_{i,j} f_ij = 1
					foreach j of local centiles {

					foreach var of varlist ftotinc_adj {
					*Create string variable for the household type
					gen `i1'_`i2'_`i3'_`i4' = (`i1' == 1 &`i2' == 1 &`i3' == 1 & `i4' == 1 & c`j'_`var' == 1)
					
					quietly sum `i1'_`i2'_`i3'_`i4' [aweight=perwt]
					scalar f_`var'_`i'_`j'_`t' = r(mean)

					quietly sum `var' if `i1'_`i2'_`i3'_`i4' == 1  [aweight=perwt]
					scalar r_`var'_`i'_`j'_`t' = r(mean) / mean_`var'_`t'

					*Create a temporary macro used only the next if statement
					local r_temp = r_`var'_`i'_`j'_`t'

					*If r is a NaN, replace with a zero
					if `r_temp' == . {

						scalar r_`var'_`i'_`j'_`t' = 0

					}
                                        
					*Need to drop variable because cannot add additional j-index (32 characters limit)
					drop `i1'_`i2'_`i3'_`i4'
					}
					}
					
					*Write the centiles to the text file
					if `t' == 1960 {
					file write types "`i1'_`i2'_`i3'_`i4'" _n
					}
					
					foreach var of varlist ftotinc_adj {
					
					file write f_`var'_`t' %16.15f (f_`var'_`i'_1_`t')  _column(20)  %16.15f (f_`var'_`i'_2_`t') ///
					_column(40) %16.15f (f_`var'_`i'_3_`t')  _column(60) %16.15f (f_`var'_`i'_4_`t') ///
					_column(80) %16.15f (f_`var'_`i'_5_`t')  _column(100) %16.15f (f_`var'_`i'_6_`t') ///
					_column(120) %16.15f (f_`var'_`i'_7_`t')  _column(140) %16.15f (f_`var'_`i'_8_`t') ///
					_column(160) %16.15f (f_`var'_`i'_9_`t')  _column(180) %16.15f (f_`var'_`i'_10_`t') _n

					file write r_`var'_`t' %16.15f (r_`var'_`i'_1_`t')  _column(20)  %16.15f (r_`var'_`i'_2_`t') ///
					_column(40) %16.15f (r_`var'_`i'_3_`t')  _column(60) %16.15f (r_`var'_`i'_4_`t') ///
					_column(80) %16.15f (r_`var'_`i'_5_`t')  _column(100) %16.15f (r_`var'_`i'_6_`t') ///
					_column(120) %16.15f (r_`var'_`i'_7_`t')  _column(140) %16.15f (r_`var'_`i'_8_`t') ///
					_column(160) %16.15f (r_`var'_`i'_9_`t')  _column(180) %16.15f (r_`var'_`i'_10_`t') _n

                                        }

					}


					if `i' < 400 { /*this refers to all married household types*/

					foreach i5 of local peducation {

						foreach i6 of local pmarket {

						*Update counter 
						local i = `i' + 1
						di "`i'"

						*Calculate the fraction of type-i household in income decile j: type-(i,j), sum_{i,j} f_ij = 1
						foreach j of local centiles {

						foreach var of varlist ftotinc_adj {
						*Create string variable for the household type
					        gen `i1'_`i2'_`i5'_`i3'_`i6'_`i4' = (`i1' == 1 &`i2' == 1 &`i3' == 1 & `i4' == 1 & `i5' == 1 & `i6' == 1 & c`j'_`var' == 1)
					
						quietly sum `i1'_`i2'_`i5'_`i3'_`i6'_`i4' [aweight=perwt]
						scalar f_`var'_`i'_`j'_`t' = r(mean)

						quietly sum `var' if `i1'_`i2'_`i5'_`i3'_`i6'_`i4' == 1  [aweight=perwt]
					        scalar r_`var'_`i'_`j'_`t' = r(mean) / mean_`var'_`t'

						*Create a temporary macro used only the next if statement
						local r_temp = r_`var'_`i'_`j'_`t'

						*If r is a NaN, replace with a zero
						if `r_temp' == . {

							scalar r_`var'_`i'_`j'_`t' = 0

						}
                                                
						*Need to drop variable because cannot add additional j-index (32 characters limit)
						drop `i1'_`i2'_`i5'_`i3'_`i6'_`i4'
						}
						}

						*Write the centiles to the text file
						if `t' == 1960 {
						file write types "`i1'_`i2'_`i5'_`i3'_`i6'_`i4'" _n
						}

						foreach var of varlist ftotinc_adj {

						file write f_`var'_`t' %16.15f (f_`var'_`i'_1_`t')  _column(20)  %16.15f (f_`var'_`i'_2_`t') ///
						_column(40) %16.15f (f_`var'_`i'_3_`t')  _column(60) %16.15f (f_`var'_`i'_4_`t') ///
						_column(80) %16.15f (f_`var'_`i'_5_`t')  _column(100) %16.15f (f_`var'_`i'_6_`t') ///
						_column(120) %16.15f (f_`var'_`i'_7_`t')  _column(140) %16.15f (f_`var'_`i'_8_`t') ///
						_column(160) %16.15f (f_`var'_`i'_9_`t')  _column(180) %16.15f (f_`var'_`i'_10_`t') _n

						file write r_`var'_`t' %16.15f (r_`var'_`i'_1_`t')  _column(20)  %16.15f (r_`var'_`i'_2_`t') ///
						_column(40) %16.15f (r_`var'_`i'_3_`t')  _column(60) %16.15f (r_`var'_`i'_4_`t') ///
						_column(80) %16.15f (r_`var'_`i'_5_`t')  _column(100) %16.15f (r_`var'_`i'_6_`t') ///
						_column(120) %16.15f (r_`var'_`i'_7_`t')  _column(140) %16.15f (r_`var'_`i'_8_`t') ///
						_column(160) %16.15f (r_`var'_`i'_9_`t')  _column(180) %16.15f (r_`var'_`i'_10_`t') _n

                                                }
						
						}

					}
					
					}
		
				}

			}

		}

	}

	*Close the text files
	if `t' == 1960 {
	file close types
	}

	foreach var of varlist ftotinc_adj {

	file close f_`var'_`t'

	file close r_`var'_`t'

	}

	*More statistics
	****************

	*More statistics: (i)  Share of wife's labor income out of total household labor income in a contingency table (Table 1)
	
	*Create the variables for the share of wive's income out of total income
	gen wife_income = .
	replace wife_income = pincwage / labinc if mar == 1 & sex == 1 & labinc > 0
	replace wife_income = 0 if mar == 1 & sex == 1 & labinc <= 0

        *Contingency table (fractions of husband-wife pairs by the education of the husband and the wife)
	foreach i1 of local education {
	
		foreach i2 of local peducation {

			*Record the cells	
			quietly sum wife_income if sex == 1 & mar==1 & `i1' == 1 & `i2' == 1  [aweight=perwt]
			scalar wife_income_`i1'_`i2'_`t' = r(mean)

		}

		*Record a vector with the distributions of wive''s income for a give male
		mat wife_income_`i1'_`t' = [wife_income_`i1'_plhs_`t',wife_income_`i1'_phs_`t',wife_income_`i1'_psmcol_`t',wife_income_`i1'_pcol_`t',wife_income_`i1'_pmcol_`t']

	}
	
	*Tables with share of wive's income out of total income by educ of husband and wife
	mat wife_income_`t' = [wife_income_lhs_`t' \ wife_income_hs_`t' \ wife_income_smcol_`t'\ wife_income_col_`t' \ wife_income_mcol_`t']
	matrix colnames wife_income_`t' = plhs phs psmcol pcol pmcol
	matrix rownames wife_income_`t' = lhs hs smcol col mcol

	*Write down results
	local educ_levels 1 2 3 4 5

	file open wife_share_`t' using "wife_share_`t'.txt", write replace

	foreach j of local educ_levels {

	file write wife_share_`t' %16.15f (wife_income_`t'[`j',1])  _column(20)  %16.15f (wife_income_`t'[`j',2]) ///
	_column(40) %16.15f (wife_income_`t'[`j',3])  _column(60) %16.15f (wife_income_`t'[`j',4]) ///
	_column(80) %16.15f (wife_income_`t'[`j',5]) _n

	}

	file close wife_share_`t'

	*More statistics: (ii) Wife's share income and LFP by centiles	
        foreach j of local centiles {

                quietly sum wife_income if c`j'_ftotinc == 1 & mar==1 [aweight=perwt]
		scalar m_wife_income_`j'_`t' = r(mean)

                quietly sum  pwrk if c`j'_ftotinc == 1 & mar == 1 [aweight=perwt]
                scalar wife_pwrk_`j'_`t' = r(mean)

        }

        *Share of wife's income out of total household labor income by deciles
        mat mat_wife_income_`t'= [m_wife_income_1_`t' \ m_wife_income_2_`t' \ m_wife_income_3_`t' \ m_wife_income_4_`t' \ ///
                                   m_wife_income_5_`t' \ m_wife_income_6_`t' \ m_wife_income_7_`t' \ m_wife_income_8_`t' \ ///
                                   m_wife_income_9_`t' \ m_wife_income_10_`t']

        *Share of wives working in each decile of laborincome
        mat wife_pwrk_`t' = [wife_pwrk_1_`t' \ wife_pwrk_2_`t' \ wife_pwrk_3_`t' \ wife_pwrk_4_`t' \ ///
                                   wife_pwrk_5_`t' \ wife_pwrk_6_`t' \ wife_pwrk_7_`t' \ wife_pwrk_8_`t' \ ///
                                   wife_pwrk_9_`t' \ wife_pwrk_10_`t']
	
}

*Write down results
file open income_share_wife using "income_share_wife.txt", write replace
file open lfp_wife using "lfp_wife.txt", write replace

foreach j of local centiles {

file write income_share_wife %16.15f (mat_wife_income_1960[`j',1])  _column(20)  %16.15f (mat_wife_income_2005[`j',1]) _n
file write lfp_wife %16.15f (wife_pwrk_1960[`j',1])  _column(20)  %16.15f (wife_pwrk_2005[`j',1]) _n

}

file close income_share_wife
file close lfp_wife

*More statistics: (iii) Contingency tables
local year_codes 1960 1970 1980 1990 2000 2005

*Main loop (years)

local counter_tloop = 0 /*loop counter*/

foreach t of local year_codes {

	
	*Change to input directory
	cd "`stata_dir'"
	cd "Input"

        use sample`t', clear

	*Update counter
	local counter_tloop = `counter_tloop' + 1
        
	*Define the universe of households
	keep if age >= 25 & age <= 54 /*age restriction*/
	keep if eldch <= 19 | eldch == 99 /*keep household with children younger than 19 or N/A children*/
	keep if yngch <= 19 | yngch == 99 /*keep household with children younger than 19 or N/A children*/
	drop if marst == 2 /*2 is spouse missing, so we cannot match these families with sploc and these will automatically have to go*/
	drop if marst == 5 /*5 is widowed, shall we drop them or treat them as singles?*/
	drop if empstat == 0 /*employment status N/A*/

        *Income levels restrictions 
	drop if ftotinc < 0 /*remove negative incomes*/
        drop if ftotinc >= 9999999 /*remove N/A values*/
	
	*Education, marital status, fertility, sex, participation
	gen lhs = educ<6 /*less than highschool*/
	gen hs = educ==6 /*highschool*/
	gen smcol = (educ>6&educ<10) /*some college*/ 	
	gen col = educ==10 /*college*/
	gen mcol = educ==11 /*college  plus*/ 
	
	gen mar = marst == 1 /*2 is spouse missing, so we cannot match these families with sploc and these will automatically have to go*/

	gen nmar = marst == 6 /*never married*/
	
	gen div = (marst == 4 | marst == 3)  /*separated included*/
	 	
	gen fem = sex == 2

	gen mal = sex == 1
		
	*Characteristics of the marital partner
	bysort serial:gen peduc =  educ[sploc]
	bysort serial:gen pempstat =  empstat[sploc]
	bysort serial:gen pincwage = incwage[sploc]

	/*drop if characteristics are missing for married people*/
	drop if mar == 1 & (peduc == . | pempstat ==. | pincwage == .)

	/*check once again whether all married have their spouses with characteristics*/
	gen plhs = peduc<6 /*less than highschool*/
	gen phs = peduc==6 /*highschool*/
	gen psmcol = (peduc>6 & peduc<10) /*some college*/ 	
	gen pcol = peduc==10 /*college*/
	gen pmcol = peduc==11 /*college  plus*/ 
	
	/*drop married females cause we are going to work only with the married males*/
	drop if mar==1 & fem==1

	/*more restrictions: remove households with subfamilies*/
	gen xfamsize = famsize - nchild
	drop if mar == 1 & xfamsize > 2 /*only the two parents should live with their children*/
	drop if nmar == 1 & xfamsize > 1 /*only the single parent should live with her children*/
	drop if div == 1 & xfamsize > 1 /*only the divorced parent should live with her children*/
			
	*Define local for different demographic types of households 
	local education lhs hs smcol col mcol	
	local peducation plhs phs psmcol pcol pmcol

	*Count the number of married households (from the view point of the husband)
	quietly count if sex == 1 & mar==1 & (plhs==1 | phs==1 | psmcol==1 | pcol==1 |pmcol==1) & (lhs==1 | hs==1 | smcol==1 | col==1 |mcol==1)
	scalar num_mar_`t' = r(N)

        *Contingency table (fractions of husband-wife pairs by the education of the husband and the wife)
	foreach i1 of local education {
	
		foreach i2 of local peducation {

			*Record the cells
			quietly sum `i1' if sex == 1 & mar==1 & `i2' == 1  [aweight=perwt]
			scalar pr_`i1'_`i2'_`t' = r(mean)

			quietly sum `i2' if sex == 1 & mar==1  [aweight=perwt]
			scalar `i2'_`t' = r(mean)

			scalar `i1'_`i2'_`t' = pr_`i1'_`i2'_`t' * `i2'_`t'

		}

		*Record a vector with the distributions of wives for a give male
		mat dist_`i1'_`t' = [`i1'_plhs_`t',`i1'_phs_`t',`i1'_psmcol_`t',`i1'_pcol_`t',`i1'_pmcol_`t']

	}

	*Contingency tables
	mat con_tab_`t' = [dist_lhs_`t' \ dist_hs_`t' \ dist_smcol_`t' \ dist_col_`t' \ dist_mcol_`t']
	matrix colnames con_tab_`t' = plhs phs psmcol pcol pmcol
	matrix rownames con_tab_`t' = lhs hs smcol col mcol

	*Write down results
	local educ_levels 1 2 3 4 5

	*Change to input directory
	cd "`stata_dir'"
	cd "Output"

	file open con_`t' using "con_`t'.txt", write replace

	foreach j of local educ_levels {

	file write con_`t' %16.15f (con_tab_`t'[`j',1])  _column(20)  %16.15f (con_tab_`t'[`j',2]) ///
	_column(40) %16.15f (con_tab_`t'[`j',3])  _column(60) %16.15f (con_tab_`t'[`j',4]) ///
	_column(80) %16.15f (con_tab_`t'[`j',5]) _n

	}

	file close con_`t'

}

*Change to stata directory
cd "`stata_dir'"
