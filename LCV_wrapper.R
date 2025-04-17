## LCV example script written by Katie Siewert sourced from https://github.com/lukejoconnor/LCV - adapted by William Reay - further adapted by Dylan Kiltschewskij

########
## Disentangling correlation from causation in the putatitve relationship between significantly genetically correlated blood based biomarkers and pneumonia
## 2020 - William Reay
## 2021 - Dylan Kiltschewskij
########

##Load dependencies
suppressMessages(library(optparse))

##Specify command line inputs
option_list = list(
  make_option("--phenotype_1", action="store", default=NA, type='character',
              help="The name of the first trait to be analysed [required]."),
  make_option("--sumstats_1", action="store", default=NA, type='character',
              help="File name for summary stats for trait 1, gzip compression required."),
  make_option("--phenotype_2", action="store", default=NA, type='character',
              help="The name of the second trait to be analysed [required]."),
  make_option("--sumstats_2", action="store", default=NA, type='character',
              help="File name for summary stats for trait 2, gzip compression required."),
  make_option("--ldscores", action="store", default=NA, type="character",
              help = "Path to LD scores [required]."),
  make_option("--path",action="store",default=NA,type="character",
              help="Path of working directory. Must contain \"RunLCV.R\" and \"MomentFunctions.R\"")
)

opt = parse_args(OptionParser(option_list=option_list))
setwd(paste(opt$path))

#Load in munged summary statistics for FVC and trait to be analysed

cat("#########################")
cat("\n")
cat("Loading summary statistics for ", paste(opt$phenotype_1), "and", paste(opt$phenotype_2))
cat("\n")
cat("#########################")
cat("\n")

#Load trait 1 data 
opt$phenotype_1_df <- read.table(file = paste(opt$sumstats_1,sep=""),
                                    header=TRUE,sep="\t",stringsAsFactors = FALSE, na.strings=c("", "NA"))

opt$phenotype_1_df <- na.omit(opt$phenotype_1_df)

#Load trait 2 data 
opt$phenotype_2_df <- read.table(file = paste(opt$sumstats_2,sep=""),
                                 header=TRUE,sep="\t",stringsAsFactors = FALSE, na.strings=c("", "NA"))

opt$phenotype_2_df <- na.omit(opt$phenotype_2_df)

cat("\n")
cat("#########################")
cat("\n")
cat("GWAS data loaded")
cat("\n")
cat("#########################")
cat("\n")

cat("\n")
cat("#########################")
cat("\n")
cat("Loading LD scores")
cat("\n")
cat("#########################")
cat("\n")

##################################### EDIT ME
#Load LD scores
LD_scores <- read.table(paste(opt$ldscores,sep=""),header=TRUE,sep="\t",stringsAsFactors=FALSE)
##################################### /EDIT ME

cat("\n")
cat("#########################")
cat("\n")
cat("Merging data frames")
cat("\n")
cat("#########################")
cat("\n")


#Merge such that SNPs are annotated with LD scores
Merged_df <- merge(LD_scores,opt$phenotype_1_df,by="SNP")
Annotated_df_combined <- merge(Merged_df,opt$phenotype_2_df,by="SNP")

#Sort by position 
Sorted_df <- Annotated_df_combined[order(Annotated_df_combined[,"CHR"],Annotated_df_combined[,"BP"]),]
Sorted_df$L2 <- as.numeric(Sorted_df$L2)

#Check if any mismatches
mismatch = which(Sorted_df$A1.x!=Sorted_df$A1.y,arr.ind=TRUE)
Sorted_df[mismatch,]$Z.y = Sorted_df[mismatch,]$Z.y*-1
Sorted_df[mismatch,]$A1.y = Sorted_df[mismatch,]$A1.x
Sorted_df[mismatch,]$A2.y = Sorted_df[mismatch,]$A2.x


cat("\n")
cat("#########################")
cat("\n")
cat("Constructing LCV model, output will be sent to text file",  paste(opt$phenotype_1, opt$phenotype_2, "LCV.txt", sep = "_"))
cat("\n")
cat("#########################")
cat("\n")


#Run LCV-need to setwd to directory containing LCV package
source("/g/data/hh57/dkiltschewskij/programs/lcv/RunLCV.R")

file = paste(opt$phenotype_1, "__", opt$phenotype_2, "_LCV.txt", sep = "")

sink(file, append = TRUE)

LCV <- RunLCV(Sorted_df$L2,Sorted_df$Z.x, Sorted_df$Z.y)
sprintf("Estimated posterior gcp=%.2f(%.2f), pvalue=%.1f; estimated rho=%.2f(%.2f)",LCV$gcp.pm, LCV$gcp.pse, LCV$pval.gcpzero.2tailed, LCV$rho.est, LCV$rho.err)

#Sink output to text file

LCV

sink()

#Clear environment after run
rm(list = ls())
