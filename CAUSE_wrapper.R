########################################
# D. Kiltschewskij, 05/08/24
#
# CAUSE causality analysis
#
########################################


# libraries
library(data.table)
library(dplyr)
library(cause)
library(genetics.binaRies)
library(ieugwasr)
library(optparse)


## command line arguments
option_list=list(
  make_option("--exp", action="store",default=NA,type="character",help="GWAS file containing SNPs for the exposure. GWAS file should NOT be munged. Must contain SNP, BETA, SE, A1, A2, P and N columns [required]"),
  make_option("--exp_name", action="store",default=NA,type="character",help="Name of the exposure [required]."),
  make_option("--exp_cols", action="store",default=NA,type="character",help="Column names (comma separated) corresponding to the following fields: SNP, EFFECT/BETA, SE, EFFECT_ALLELE, OTHER_ALLELE, P [required]"),
  make_option("--oc", action="store",default=NA,type="character",help="Outcome GWAS file. Should NOT be munged. Must contain SNP, BETA, SE, A1, A2 and P columns [required]"),
  make_option("--oc_name", action="store",default=NA,type="character",help="Name of outcome [required]."),
  make_option("--oc_cols", action="store",default=NA,type="character",help="Column names (comma separated) corresponding to the following fields: SNP, EFFECT/BETA, SE, EFFECT_ALLELE, OTHER_ALLELE, P [required]"),
  make_option("--OR", action="store",default="",type="character",help="Flag indicating whether effect sizes are reported as an odds ratio. Can be either \"exp\", \"oc\" or \"both\". Leave as default if neither GWAS have ORs"),
  make_option("--R2", action="store", default=0.01,type="numeric",help="R2 LD clumping threshold."),
  make_option("--P", action="store", default=1e-3,type="numeric",help="P value LD clumping threshold."),
  make_option("--ref", action="store", default=NA,type="character",help="Path to g1000 LD reference [required]."),
  make_option("--print_all", action="store", default="No",type="character",help="Print all auxillary files? Can be either \"Yes\" or \"No\" (default) [required]."),
  make_option("--path", action="store", default=".",type="character",help="Path to working directory [required].")
)
opt=parse_args(OptionParser(option_list=option_list))
setwd(paste(opt$path))


# convert column flags into vectors
exposure_cols<-strsplit(paste(opt$exp_cols,sep=""),",") %>% unlist
outcome_cols<-strsplit(paste(opt$oc_cols,sep=""),",") %>% unlist


# read in exposure
cat("####################")
cat("\n")
cat(paste("Reading in exposure: ",opt$exp_name, " using the following file: ", opt$exp, sep=""))
cat("\n")
exposure<-fread(opt$exp, data.table=F)


# read in outcome
cat("####################")
cat("\n")
cat(paste("Reading in outcome: ",opt$oc_name, " using the following file: ", opt$oc, sep=""))
cat("\n")
outcome<-fread(opt$oc, data.table=F)


# call columns of interest
cat("####################")
cat("\n")
cat(paste("Calling SNP, EFFECT/BETA, SE, EFFECT_ALLELE, OTHER_ALLELE, P for ",opt$exp,sep=""))
cat("\n")
cat(paste("Using the following column names: ",opt$exp_cols,sep=""))
cat("\n")
exposure<-exposure[,exposure_cols]

cat("####################")
cat("\n")
cat(paste("Calling SNP, EFFECT/BETA, SE, EFFECT_ALLELE, OTHER_ALLELE, P for ",opt$out,sep=""))
cat("\n")
cat(paste("Using the following column names: ",opt$oc_cols,sep=""))
cat("\n")
outcome<-outcome[,outcome_cols]


# generate BETA for binary outcome if required
if (opt$OR=="exp"){
  cat("####################")
  cat("\n")
  cat("Converting exposure ORs to BETAs")
  exposure$BETA<-log(exposure[,exposure_cols[2]]) # position 2 in the column vector should always correspond to the effect size
  exposure<-exposure[,c(1,7,3:6)]
} else if (opt$OR=="oc"){
  cat("####################")
  cat("\n")
  cat("Converting outcome ORs to BETAs")
  cat("\n")
  outcome$BETA<-log(outcome[,outcome_cols[2]])
  outcome<-outcome[,c(1,7,3:6)]
} else if (opt$OR == "both"){
  cat("####################")
  cat("\n")
  cat("Converting exposure and outcome ORs to BETAs")
  cat("\n")
  exposure$BETA<-log(exposure[,exposure_cols[2]])
  exposure<-exposure[,c(1,7,3:6)]
  outcome$BETA<-log(outcome[,outcome_cols[2]])
  outcome<-outcome[,c(1,7,3:6)]
}


# standardise column names
cat("####################")
cat("\n")
cat("Standardising column names")
cat("\n")
colnames(exposure)<-c("SNP","BETA","SE","EFFECT","OTHER","P")
colnames(outcome)<-c("SNP","BETA","SE","EFFECT","OTHER","P")


# merge GWAS
cat("####################")
cat("\n")
cat("Standardising column names")
cat("\n")
harmonised<-gwas_merge(exposure,outcome,
                       snp_name_cols=c("SNP","SNP"),
                       beta_hat_cols = c("BETA","BETA"),
                       se_cols = c("SE","SE"),
                       A1_cols = c("EFFECT","EFFECT"),
                       A2_cols = c("OTHER","OTHER"),
                       pval_cols=c("P","P"))


# calculate nuisance parameters
cat("####################")
cat("\n")
cat("Calculating nuisance parameters")
cat("\n")
set.seed(0)
varlist <- with(harmonised, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(harmonised, varlist)


# LD clump
cat("####################")
cat("\n")
cat("LD clumping/calling top clumped variants")
cat("\n")
r2_thresh<-opt$R2
pval_thresh<-opt$P

clumped<-harmonised %>% 
  rename(rsid = snp,pval = p1) %>% 
  ieugwasr::ld_clump(dat = .,clump_r2 = r2_thresh,clump_p = pval_thresh,plink_bin = genetics.binaRies::get_plink_binary(),pop = "EUR",bfile = opt$ref)
top_vars <- clumped$rsid


# fit CAUSE
cat("####################")
cat("\n")
cat("Fitting CAUSE")
cat("\n")

fit.cause<-function(qa,qb,modelname){
  # fit CAUSE
  res <- cause(X=harmonised,variants=top_vars,param_ests=params,qalpha = qa,qbeta = qb)
  # clean results
  res2<-data.frame(res$elpd)
  res2$exposure<-opt$exp_name
  res2$outcome<-opt$oc_name
  res2$betadist<-paste(modelname)
  res2$p<-pnorm(abs(res2$z),lower.tail = F)
  
  res3<-summary(res, ci_size=0.95)$tab %>% data.frame()
  res3$exposure<-opt$exp_name
  res3$outcome<-opt$oc_name
  res3$betadist<-paste(modelname)

  ## save results
  if (opt$print_all == "Yes"){
    
  # params mix grid 
    write.table(x = params$mix_grid,file=paste(opt$exp_name,"_vs_",opt$oc_name,"_",modelname,"_cause_params.txt",sep=""),sep="\t",row.names = F,quote=F)
    
  # causal estimate plot
    png(filename=paste(opt$exp_name,"_vs_",opt$oc_name,"_",modelname,"_cause.png",sep=""),width = 15,height=3.75,units = "in",res = 800)
    plot(res, type="data")
    png()
    dev.off()
    dev.off()
  }

  # loos
  sink(paste(opt$exp_name,"_vs_",opt$oc_name,"_",modelname,"_cause_loos.txt",sep=""))
  print(res$loos[[2]])
  print(res$loos[[3]])
  sink()

  # model estimates and model summary
  write.table(x = res2,file=paste(opt$exp_name,"_vs_",opt$oc_name,"_",modelname,"_cause_elpd.txt",sep=""),sep="\t",row.names = F,quote=F)
  write.table(x = res3,file=paste(opt$exp_name,"_vs_",opt$oc_name,"_",modelname,"_cause_est.txt",sep=""),sep="\t",row.names = F,quote=F)
}

fit.cause(1,2,"a1_b2")
fit.cause(1,10,"a1_b10")
fit.cause(1,50,"a1_b50")
