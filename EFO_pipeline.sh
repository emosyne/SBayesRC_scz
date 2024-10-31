#!/bin/bash 

#SBATCH -J SBayesRC_scz
#SBATCH -A MURRAY-SL2-CPU
#SBATCH -p cclake
#SBATCH --nodes=1
#! The Cascade Lake (cclake) nodes have 56 CPUs (cores) each and
#! 3420 MiB of memory per CPU.
#SBATCH --ntasks=56
#! # SBATCH --mem=32G
#SBATCH --time=2:00:00

#! sbatch directives end here (put any additional directives above this line)




###SCRIPT STARTS FROM HERE


module load miniconda/3
conda activate SBayesRC_scz_env

# This is a pipeline how I format any GWAS summary-level data into cojo format, and use SBayesRC method to generate predictors. The method is implemented in [GCTB](https://cnsgenomics.com/software/gctb/#SBayesRCTutorial). 
# We also included example code to run SBayesR, SBayesS, COJO and clumping on this page. 
# They are put in a pipeline with the tool [qsubshcom](https://github.com/zhilizheng/qsubshcom).


# # Preparation

# Resource data is available for local users with the following path setting, and public for [downloading](https://cnsgenomics.com/software/gctb/#Download) too. 


# cd /Users/eosimo/GoogleDrive/WORK/CF_PhD/SBayesRC_scz
cd /home/efo22/SBayesRC_scz

## this is where you placed file cojo_format_v7.R
exedir="./"
# input_files="/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/private_input_files"
input_files="/home/efo22/private_input_files"

## this is where you placed LD files and annotation file. They are available for downloading from GCTB website:
## https://cnsgenomics.com/software/gctb/#Download
annot="${input_files}/SBayes_annots/annot_baseline2.2.txt" # annotation file
ldm1=${input_files}/LD_ref/ukbEUR_Imputed/
ldm2=${input_files}/LD_ref/ukbEUR_HM3/





# For a neat organization, you can put your GWAS file into a folder and name it with the trait. All the outputs of following analysis will be put into the same folder. 
# Several sub-folders will be set for different analysis. 


trait=SCZ
gwas_file=${input_files}/GWAS/PGC3_SCZ_wave3.primary.autosome.public.v3_IMPover8.tsv

## choose the LD reference based on number of SNPs in GWAS data. 
if [ $(wc -l  ${gwas_file}  | awk '{print $1}'  ) -gt  5149563 ]; then   ldm=$ldm1 ; else  ldm=$ldm2; fi

mkdir ${trait} ${trait}/COJO ${trait}/COJO_imp ${trait}/CplusT ${trait}/SBayesR ${trait}/SBayesS  ${trait}/SBayesRC
mkdir ${trait}/SBayesRC/GCTB_v2 

echo " there are " $(wc -l  ${gwas_file} | awk '{print $1}' ) "SNPs in the original GWAS data"



# # COJO format
# GWAS summary level data has various formats. 
# We will need to format it to cojo format to use it as a proper input for GCTB. 
# for example, if the raw gwas file has a header like this:

# > MarkerName A1 A2 Freq LogOR StdErrLogOR P

# you can put each element of the header name into the script per flag:
head -n 1 ${gwas_file}

snpinfo="${input_files}/LD_ref/ukbEUR_HM3/snp.info"
ma_file=${trait}/GWAS

Rscript  ${exedir}/cojo_format_v7.R  \
  --file  ${gwas_file}  \
  --out  ${ma_file}.cojo \
  --SNP  SNP  \
  --A1 A1  \
  --A2  A2 \
  --freq  freq   \
  --pvalue p  \
  --beta  b  \
  --se  se   \
  --samplesize  N \
  --ref ${snpinfo}



# # QC and Impute


gctb --ldm-eigen $ldm --gwas-summary ${ma_file}.cojo --impute-summary --out ${ma_file}_imp.ma --thread 50

# # impute by block
# # for i in {1..591}
# # do
# #   echo "block: $i"
# #   gctb --ldm-eigen $ldm --gwas-summary ${ma_file}.cojo --impute-summary --out ${ma_file}_block --block $i --thread 4
# # done
# # # This command will generate a text file ${ma_file}_block for each block. After all blocks are finished, use the following command to merge all .ma files across block into a single file with genome-wide SNPs:
# # gctb --gwas-summary ${ma_file}_block --merge-block-gwas-summary --out ${ma_file}_allblocks_imp.ma


# # SBayesRC


# gctb   \
# --sbayes RC  \
# --ldm-eigen   ${ldm}   \
# --gwas-summary   ${ma_file}_imp.cojo   \
# --annot  $annot  \
# --out  ${trait}/SBayesRC/GCTB_v2/${ma_file}_imp.cojo_imputed_sbrc_gctb   \
# --thread 4



# # plot SBayesRC effect size vs. marginal effect size

# At last we compare the marginal effect size with the effect size from SBayesRC with a simple plot. 


# plotcmd=" Rscript  ${exedir}/effect_size_plot_for_GCTB_output.R    $trait   ${gwas_file}_imp.cojo.imputed.cojo    ${trait}/SBayesRC/GCTB_v2/${gwas_file}_imp.cojo.imputed_sbrc_gctb   "
# jobname="effect_plot_"${trait}
# plotsub=`qsubshcom "$plotcmd"  1  50G  $jobname  1:00:00  " -wait=$sbrc_gctb_sub  " `


# As an example:

# <img src="Anorexia_01_pgcAN2.2019-07.modified.vcf.tsv_sbrc.txt_compare_marginal_effect_vs_SBayesRC_20231103_10_18.png" width="50%" height="50%" />


# # Genetic Variance


# Rscript  quick_Vg_report_from_snpRes_file.R  ${PGS_file}.nopred  ${predictor}.snpRes


# # profile PGS


# plink \
#    --bfile  $input \
#    --score  $predictor  2 5 8  header sum    \
#    --out ${outdir}/${cohort}_${trait}_SBayesRC



# # Clumping

# Clumping plus threshold method could be used as a baseline model to compare the prediction accuracy.


# snps2exclude=2_duplicated_SNPs_to_exclude.txt          

# clmp_sub=`qsubshcom "plink --bfile UKB_20k/PLINK/ukbEURu_imp_chr{TASK_ID}_v3_impQC_20k \
#   --exclude $snps2exclude \
#   --clump ${trait}/${gwas_file}.cojo \
#   --clump-p1 0.05  \
#   --clump-p2 0.05  \
#   --clump-r2 0.1 \
#   --clump-kb 500 \
#   --clump-field "p" \
#   --out ${trait}/CplusT/${gwas_file}_chr{TASK_ID}_clumped   "  10 150G "clumping" 24:00:00 "  -wait=$formatqsub  -array=1-22 "   `
                


# # COJO
# We can run COJO to see how many significant and independent SNP signal were in the data. 


# cojo_sub=`qsubshcom "gcta-1.94.1  \
# --bfile   UKB_20k/PLINK/ukbEURu_imp_chr{TASK_ID}_v3_impQC_20k \
# --chr   {TASK_ID} \
# --cojo-file   ${trait}/${gwas_file}.cojo  \
# --cojo-slct  \
# --out  ${trait}/COJO/${trait}_chr{TASK_ID}_cojo  "  10 150G "COJO"  24:00:00 "  -wait=$formatqsub  -array=1-22 "   `


# # SBayesR


# {bash, eval =F}
# sbr_gctb_sub=`qsubshcom "gctb  --sbayes R  \
# --ldm-eigen  ${ldm} \
# --gwas-summary  ${trait}/${gwas_file}_imp.cojo.imputed.cojo  \
# --chain-length 5000 --burn-in 2000   --no-mcmc-bin   \
# --out  ${trait}/SBayesR/${trait}_SBayesR_eigen  --thread 10 "  10 150G  "SBayesR" 24:00:00 " -wait=$gctbimputesub"   `


# # SBayesS

# {bash, eval =F}
# sbs_gctb_sub=`qsubshcom "gctb  --sbayes S  \
# --ldm-eigen  ${ldm} \
# --gwas-summary ${trait}/${gwas_file}_imp.cojo.imputed.cojo  \
# --chain-length 5000 --burn-in 2000 --num-chains 3  --no-mcmc-bin   \
# --out  ${trait}/SBayesS/${trait}_SBayesS_eigen --thread 10 "  10 150G "SBayesS" 24:00:00 " -wait=$gctbimputesub"   `




# # Useful links:

# GCTB:  
# > https://cnsgenomics.com/software/gctb/#SBayesRCTutorial

# R version SBayesRC:  
# The same SBayesRC methods is also available in an R version, which is well explained at Zhili's Github.
# Although the output files appear in different format as the GCTB version, and the content of output are different, the core computation is the same. 

# > https://github.com/zhilizheng/SBayesRC  

# qsubshcom:  
# > https://github.com/zhilizheng/qsubshcom

# COJO:
# > https://yanglab.westlake.edu.cn/software/gcta/#COJO

# Clump:
# > https://zzz.bwh.harvard.edu/plink/clump.shtml


