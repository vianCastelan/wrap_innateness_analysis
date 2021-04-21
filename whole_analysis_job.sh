#!/bin/bash

#######################
# innateness analysis #
#######################
#PBS -N whole_analysis_job
#PBS -o /mnt/BioAdHoc/Groups/vd-vijay/vcastelan/innateness_analysis/lmm_test_1000mvg_PC1/wrap/whole_analysis_job_out.txt
#PBS -e /mnt/BioAdHoc/Groups/vd-vijay/vcastelan/innateness_analysis/lmm_test_1000mvg_PC1/wrap/whole_analysis_job_err.txt
#PBS -m abe
#PBS -M vcastelan@lji.org
#PBS -q default
#PBS -l nodes=1:ppn=5
#PBS -l mem=10gb
#PBS -l walltime=48:00:00

OUT_DIRECTORY=${OUT_DIRECTORY:-/mnt/BioAdHoc/Groups/vd-vijay/vcastelan/innateness_analysis/lmm_test_1000mvg_PC1/wrap} #working directory
MASTER=${MASTER:-master_rpmk_table_flavell_lab.tsv} #genes counts table with all cells (YOUR DATA)
BETA=${BETA:-beta_table_1000mvg_PC1_b2.tsv} #table with innateness scores per gene (leave as default to use genes with abs(beta>2) or "results_mouse_beta_table_1000mvg_PC1.tsv" to use all genes

cd ${OUT_DIRECTORY}
echo "Generating beta tables"
python h_innate.py ${MASTER} ${BETA} 
/share/apps/R/3.6.1/bin/Rscript ${OUT_DIRECTORY}/whole_analysis.R

echo "Job completed! \n Check for errors if any."

