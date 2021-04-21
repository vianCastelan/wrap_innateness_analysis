#!/home/gascui/miniconda3/bin/python
# libraries
import pandas as pd
import sys
import os

# parameters
innate_path = os.getcwd()
cells_path = '/read_table/'
mouse_beta_path = '/mouse_beta_innateness/'
# read tsv files
if __name__ == '__main__':
    table_path = sys.argv[1]
    cell_types = pd.read_table(innate_path + cells_path + table_path, index_col = 0, delimiter="\t")
    beta_path = sys.argv[2]
    innate = pd.read_table(innate_path + mouse_beta_path + beta_path, index_col = 0, delimiter="\t")
#
print("reading files from:",innate_path+cells_path)
print("The dimension of your expression matrix are:", cell_types.shape)
print("The dimension of your beta matrix are:", innate.shape)
#
cell_types = cell_types.set_index('gene')
cell_types = cell_types.drop(['Excel_gene'], axis=1)
# filtering
vc = innate['gene'].value_counts()
vc[vc > 1].shape
innate_filt = innate[innate['gene'].isin(vc[vc ==1].index.tolist())]
tdf = innate_filt.set_index("gene")[['beta']]
print(tdf.head())
cell_types_filt = cell_types[cell_types.index.isin(tdf.index.tolist())]
print("filtered genes x samples:", cell_types_filt.shape)
cell_types_filt = cell_types_filt.join(tdf, how = 'left')
# multiplication for all genes: woosh!
#‌　 ∧＿∧　
#（。·ω·。)つ━☆·*。
#⊂　　 ノ 　　　·゜+.
#　しーＪ　　　°。+ *´¨)
#　　　　　　　　　.· ´¸.·*´¨) ¸.·*¨)
#　　　　　　　　　　(¸.·´ (¸.·'*  ☆
dfs = []
for c in cell_types_filt.columns.tolist():
    if c != "beta":
        df = cell_types_filt[c] * cell_types_filt["beta"]
        df = df.to_frame(c)
        dfs.append(df)
cell_types_filt_trans = pd.concat(dfs, axis = 1)
# aggregate
aggregate_score_cell_types = cell_types_filt_trans.sum(axis = 0)
aggregate_score_cell_types = aggregate_score_cell_types.to_frame("innateness_score")
aggregate_score_cell_types['category'] = "cell_types"
aggregate_score_cell_types["cell_type"] = aggregate_score_cell_types.index
aggregate_score_cell_types = aggregate_score_cell_types.reset_index(drop = True)
#
print("This is what the output looks like: ")
print(aggregate_score_cell_types.head())
# export data
os.system('ls -l')
fn_out = 'export/innateness_python.tsv'
aggregate_score_cell_types.to_csv(fn_out, sep = '\t', index = False)
# python h_innate.py master_rpmk_table_filtered_2020-10-12_23-32-PM.tsv innateness_beta_mod_unique_2020-10-15_15-23-PM.tsv
