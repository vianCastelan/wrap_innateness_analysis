1. Copy the count matrix of the data from which you want to calculate the innateness score to "read_table" folder.
2. Fill all the blanks of whole_analysis_job.sh script (OUT_DIRECTORY,MASTER and BETA fields) -Leave BETA as default unless you want to try all genes-.
3. Fill your desired output directory in the "output_dir" option, in the "whole_analysis.R" script
4. Run the analysis as follows:
qsub whole_analysis_job.sh 

