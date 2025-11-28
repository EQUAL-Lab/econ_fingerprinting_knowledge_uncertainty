**** Replication script for the econ fingerprinting manuscript 
** Author: Luis Sanchez (ls2252@cornell.edu) & Martina Occelli (mo386@cornell.edu)
* This code replicates Table 1-2-3 of the manuscript


* Before running the code, replace the library below with your local library
log using "C:\Users\ls2252\OneDrive - Cornell University\Sergio Run\fingerprintingrun.log", replace 

//////////////////////////////  Program start //////////////////////////////

*** Setup: 
clear all
drop _all
program drop _all
set more off
set trace off
set mem 500m
ssc install outreg2, replace

*** Data and variables creation: 
import excel "C:\Users\ls2252\OneDrive - Cornell University\Sergio Run\fingerprint_updated_data FINAL.xlsx", sheet("Sheet1") firstrow


encode region, g(rgn)
encode district, g(dstrct)

rename farm_code farm_id

sum index_final_value, d
scalar indexmean = round(r(mean), 0.01) //Knowledge mean value
dis indexmean

scalar indexmedian = round(r(p50), 0.01) //Knowledge median value
dis indexmedian

gen median_index = 0
replace median_index = 1 if index_final_value >= indexmedian
tab median_index

gen goodsource = 0
replace goodsource = 1 if  source=="Reliable" // Formal source of the seed

replace area = "" if area == "NA"
destring area, gen(plot_size)

recode gender (2=0) // 1=Women, 0=Men

gen firsttime=1 if first_time_seed=="Yes" // No replanting
replace firsttime=0 if first_time_seed=="No" // Yes replanting 

tab farmer_clean_name

gen landraces=0
replace landraces=1 if farmer_clean_name=="landraces"
replace landraces=1 if dna_identified_first=="sacapobres"
replace landraces=1 if dna_identified_first=="generalito"
replace landraces=1 if dna_identified_first=="mantequilla"

gen  trialseed=0
replace trialseed=1 if farmer_clean_name=="ensayo"
replace trialseed=1 if dna_identified_first=="RCB 593"
replace trialseed=1 if dna_identified_first=="SEF 60"
replace trialseed=1 if dna_identified_first=="SEF 64"
replace trialseed=1 if dna_identified_first=="SEF 62"
replace trialseed=1 if dna_identified_first=="SEF 42"
replace trialseed=1 if dna_identified_first=="uran"

gen  other=0
replace other=1 if farmer_clean_name=="other"

gen  dontknow=0
replace dontknow=1 if farmer_clean_name=="dontknow"

gen first_below90=0
replace first_below90=1 if first_gen_distance<0.1

gen first_below95=0
replace first_below95=1 if first_gen_distance<0.05

*** Variables:
global indexval median_index // High knowledge dummy
global covars tricot inder goodsource education gender // List of socioeconomic controls
global checks other dontknow landraces trialseed firsttime //variety-specific controls

****************************** Regs ******************************
*******************************************************************

global mtch match

*** Table 1: Treatment Effects across specifications

preserve
replace match=0 if first_gen_distance>=0.1 // 90% similarity restriction

areg $mtch 1.treated#i.${indexval}, a(rgn) cluster(farm_id)
sum match if e(sample) & treated==0 // mean outcome for control
outreg2 using regs.xls, replace dec(2)

areg $mtch 1.treated#i.${indexval} $checks, a(rgn) cluster(farm_id)
outreg2 using regs.xls, append dec(2) drop($checks)
sum match if e(sample) & treated==0

areg $mtch 1.treated#i.${indexval} $checks $covars, a(rgn) cluster(farm_id)
outreg2 using regs.xls, append dec(2) drop($checks $covars)
sum match if e(sample) & treated==0 

restore


*** Table 2: Match Rates and Precision > 90%

gen match1=match
gen refSD= sd_ref_distance_1*100
gen standsefSD = (sd_ref_distance_1-.0020415)/.0039034 
gen mean_ref_distance = max_ref_distance_1-sd_ref_distance_1 
gen meanSD = mean_ref_distance/sd_ref_distance_1

global unp invFCI_SP standsefSD 

reg $mtch invFCI_SP, cluster(farm_id)
sum $mtch if invFCI_SP==0 & e(sample)
sum $mtch if invFCI_SP==1 & e(sample)

preserve
replace match=0 if first_gen_distance>=0.10
areg match 1.treated#i.${indexval} $unp $checks, a(rgn) cluster(farm_id)
outreg2 using regs.xls, replace dec(2) drop($checks)
sum match first_similarity if e(sample) & treated==0 

areg match1 1.treated#i.${indexval} $unp $checks, a(rgn) cluster(farm_id)
outreg2 using regs.xls, append dec(2) drop($checks)
sum match1 first_similarity if e(sample) & treated==0 

replace match2=0 if second_gen_distance>=0.10
areg match2 1.treated#i.${indexval} $unp $checks, a(rgn) cluster(farm_id)
outreg2 using regs.xls, append dec(2) drop($checks)
sum match2 second_similarity  if e(sample) & treated==0 

replace match3=0 if third_gen_distance>=0.10
areg match3 1.treated#i.${indexval} $unp $checks, a(rgn) cluster(farm_id)
outreg2 using regs.xls, append dec(2) drop($checks)
sum match3 third_similarity if e(sample)  & treated==0 

replace match4=0 if fourth_gen_distance>=0.10
areg match4 1.treated#i.${indexval} $unp $checks, a(rgn) cluster(farm_id)
outreg2 using regs.xls, append dec(2) drop($checks)
sum match4 fourth_similarity  if e(sample)  & treated==0 
restore

*** Table 2: Match Rates and Precision < 90%

preserve
replace match=0 if first_gen_distance>=0.10
areg match 1.treated#i.${indexval} $unp $checks, a(rgn) cluster(farm_id)
outreg2 using regs.xls, replace dec(2) drop($checks)
sum match first_similarity if e(sample) & treated==0 

areg match1 1.treated#i.${indexval} $unp $checks, a(rgn) cluster(farm_id)
outreg2 using regs.xls, append dec(2) drop($checks)
sum match1 first_similarity if e(sample) & treated==0 


areg match2 1.treated#i.${indexval} $unp $checks, a(rgn) cluster(farm_id)
outreg2 using regs.xls, append dec(2) drop($checks)
sum match2 second_similarity  if e(sample) & treated==0 


areg match3 1.treated#i.${indexval} $unp $checks, a(rgn) cluster(farm_id)
outreg2 using regs.xls, append dec(2) drop($checks)
sum match3 third_similarity if e(sample)  & treated==0 

areg match4 1.treated#i.${indexval} $unp $checks, a(rgn) cluster(farm_id)
outreg2 using regs.xls, append dec(2) drop($checks)
sum match4 fourth_similarity  if e(sample)  & treated==0 
restore


*** Table 3: Determinants of knowledge

* Weights: in case we worry that farms with more plots may impact results.
* 		   Adding these weights helps to mantain representatives.

bysort farm_id: gen n_plots = _N
gen ws = 1 / n_plots
destring count_neigh_1km, generate(count_neigh_1km_num) ignore("NA . ,")
destring count_match_5km, generate(count_match_5km_num) ignore("NA . ,")
destring Ext, generate(Ext_num) ignore("NA . ,")


* Knowledge Index as is
areg index_final_value treated $covars [pweight=ws], a(rgn) cluster(farm_id) //knowledge index as is
outreg2 using regs.xls, replace dec(2)
areg index_final_value treated count_neigh_1km_num Ext_num $covars [pweight=ws], a(rgn) cluster(farm_id) //knowledge index as is
outreg2 using regs.xls, append dec(2)
areg index_final_value treated count_match_5km_num Ext_num $covars [pweight=ws], a(rgn) cluster(farm_id) //knowledge index as is
outreg2 using regs.xls, append dec(2)

log close
