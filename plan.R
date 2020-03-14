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

# source all files
source(file.path("R", "Functions.R"))

### plan

plan <- drake_plan(
   bdd = read.table('data/JdD_MB.csv',h=T,sep=';')
   )

# Make plan
make(plan)

# plot plan
config <- drake_config(plan)
vis_drake_graph(config)

### TODO


