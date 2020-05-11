## Plan to run analysis

# packages
library(lme4)
require(lattice)
library(stringr)
library(merTools)
library(ggplot2)
library(ade4)
library(factoextra)
require(dplyr)
require(drake)
library(lme4)

dir.create("figures", showWarnings = FALSE)

# source all files
source(file.path("R", "Functions.R"))

### plan

plan <- drake_plan(
   upi = read.table('data/UPI.csv',h=TRUE,sep=';', dec = ","),
   bdd = read.table('data/JdD_MB.csv',h=T,sep=';'),
   vec_sp = c("ABAL", "ACPS", "FASY", "PIAB"),
   df_br = format_browsing_rate(bdd, upi),
   df_br_NG = format_browsing_rate_NoGalliv(df_br), 
   list_res_br = run_analysis_browsing_rate(df_br, vec_sp),
   coef_br = get_coef_browsing_rate(list_res_br, vec_sp),
   table_coef_br = format_table_coef_br(list_res_br, vec_sp),
   AIC_br = get_AIC_browsing_rate(list_res, vec_sp),
   DeltaAIC_br = get_deltaAIC_browsing_rate(AIC_br),
   BR_p_site = get_browsing_rate_per_site(df_br, vec_sp),
   pred_BR = get_pred_ci_browsing_rate(list_res_br, df_br, vec_sp),
   plot_BR = plot_browsing_rate(BR_p_site, pred_BR),
   list_res_br_NG = run_analysis_browsing_rate(df_br_NG, vec_sp),
   BR_p_site_NG = get_browsing_rate_per_site(df_br_NG, vec_sp),
   pred_BR_NG = get_pred_ci_browsing_rate(list_res_br_NG, df_br_NG, vec_sp),
   plot_BR_PA_NG = plot_browsing_rate_PA_NoGalliv(BR_p_site_NG, pred_BR_NG),
   df_gr = format_growth(df_br),
   df_gr_NG = format_growth(df_br_NG),
   list_res_gr = run_analysis_growth(df_gr, vec_sp),
   list_res_gr_NG = run_analysis_growth(df_gr_NG, vec_sp),
   data_coef_gr = format_coef_growth(list_res_gr, vec_sp),
   table_coef_gr = format_table_coef_growth(list_res_gr, vec_sp),
   data_coef_gr_NG = format_coef_growth(list_res_gr_NG, vec_sp),
   table_coef_gr_NG = format_table_coef_growth(list_res_gr_NG, vec_sp),
   plot_growth_all = plot_growth(data_coef_gr),
   plot_growth_PA_NG = plot_growth_PA_SG(data_coef_gr_NG),
   ms = rmarkdown::render(
     knitr_in("ms.Rmd"),
     output_file = file_out("ms.pdf"),
     quiet = TRUE
   )
)


# Make plan
make(plan)


# plot plan
config <- drake_config(plan)
vis_drake_graph(config)

### TODO


