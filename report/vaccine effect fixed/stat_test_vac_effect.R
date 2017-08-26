setwd("/Users/sai/")

#50% MDA coverage, correct vaccination
vac_off <- read.csv("~/OneDrive/MORU/Projects/PSA of Savannakhet model/truePSA/output_Inc_Prev_seas/result/vaccination off on correct 1/Intervention_1_2017-07-29 12_47_33.csv")
vac_on <- read.csv("~/OneDrive/MORU/Projects/PSA of Savannakhet model/truePSA/output_Inc_Prev_seas/result/vaccination off on correct 1/Intervention_2_2017-07-29 17_25_17.csv")


#incidence
t.test(vac_off$detected_incidence,vac_on$detected_incidence,paired=TRUE) #0.1187317 0.1242261

inc_difference <- vac_off$detected_incidence-vac_on$detected_incidence

mean(inc_difference)

#prevalence
t.test(vac_off$prevalence,vac_on$prevalence,paired=TRUE) #0.07847784 0.08178142

prev_difference <- vac_off$prevalence-vac_on$prevalence

mean(prev_difference)

#0% MDA coverage, correct vaccination
vac_off <- read.csv("~/OneDrive/MORU/Projects/PSA of Savannakhet model/truePSA/output_Inc_Prev_seas/result/0 MDA coverage/Intervention_1_2017-07-29 22_50_39.csv")
vac_on <- read.csv("~/OneDrive/MORU/Projects/PSA of Savannakhet model/truePSA/output_Inc_Prev_seas/result/0 MDA coverage/Intervention_2_2017-07-30 03_27_19.csv")


#incidence
t.test(vac_off$detected_incidence,vac_on$detected_incidence,paired=TRUE) 
#0%#0.1222528 0.1278804
#50%#0.1187317 0.1242261

inc_difference <- vac_off$detected_incidence-vac_on$detected_incidence

mean(inc_difference)

#prevalence
t.test(vac_off$prevalence,vac_on$prevalence,paired=TRUE) 
#0%#0.08135902 0.08474193
#50%#0.07847784 0.08178142

prev_difference <- vac_off$prevalence-vac_on$prevalence

mean(prev_difference)

