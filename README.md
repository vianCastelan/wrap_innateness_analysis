1. Clone this repository to your working directory: git clone https://github.com/vianCastelan/wrap_innateness_analysis.git 
2. Enter to the recently cloned folder in your terminal: cd wrap_innateness_analysis
4. Copy the count matrix of the data from which you want to calculate the innateness score to "read_table" folder.
5. Fill all the blanks of whole_analysis_job.sh script (OUT_DIRECTORY,MASTER and BETA fields) -Leave BETA as default unless you want to try all genes-.
6. Fill your desired output directory in the "output_dir" option, in the "whole_analysis.R" script
7. Run the analysis as follows:
qsub whole_analysis_job.sh 

