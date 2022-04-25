#Load packages
library(TwoSampleMR)
library(MRInstruments)
library(RadialMR)

#MR Case study: Effect of BMI on CHD with SNPs as IVs
#Extract BMI and CHD data from MR Base
#Clump to remove SNPs in LD and harmonise the SNP direction
bmi_exp_dat <- extract_instruments(outcomes='ieu-a-2') # Body mass index
chd_out_dat <- extract_outcome_data(snps = bmi_exp_dat$SNP, 
                                    outcomes = 'ieu-a-7') # Coronary heart disease
dat <- harmonise_data(bmi_exp_dat, chd_out_dat)

#Conduct the MR analysis 
mr_ivw <- mr(dat, method_list=c("mr_ivw"))
mr_sensitivity <- mr(dat, method_list=c("mr_egger_regression", 
                                        "mr_ivw", 
                                        "mr_weighted_median", 
                                        "mr_weighted_mode"))
res_single <- mr_singlesnp(dat)
res_loo <- mr_leaveoneout(dat)
mr_heterogeneity(dat)

#Plot figures to illustrate the results obtained above
fig_ivw <- mr_scatter_plot(mr_ivw, dat)
fig_sensitivity <- mr_scatter_plot(mr_sensitivity, dat)
fig_forest <- mr_forest_plot(res_single)
fig_loo <- mr_leaveoneout_plot(res_loo)

#Conduct radial MR to identify outlying SNPs
radial <- ivw_radial(dat)
fig_radial <- plot_radial(radial)
