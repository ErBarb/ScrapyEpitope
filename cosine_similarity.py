from scipy import spatial 
from scipy import stats
import random
import numpy as np
import seaborn as sn
import csv
import matplotlib.pyplot as plt
import pandas as pd
import researchpy as rp
from functools import reduce



#Protparam
df_protparam_pred = pd.read_csv('pred_analysis_results/protparam_analysis.csv', index_col=False)
df_expasy_pred = pd.read_csv('pred_analysis_results/expasy_analysis.csv', index_col=False)
# df_toxinpred_pred = pd.read_csv('pred_analysis_results/toxinpred_processed.csv', index_col=False)
# df_aa_pred = pd.read_csv('pred_analysis_results/aa_composition.csv', index_col=False)
# df_atomic_pred = pd.read_csv('pred_analysis_results/atomic_composition.csv')
df_vaxijen_pred = pd.read_csv('pred_analysis_results/vaxijen_analysis.csv', index_col=False)
# df_pepstats_pred = pd.read_csv('pred_analysis_results/pepstats_analysis.csv', index_col=False)
# print(df_protparam_pred.info())
# print(df_expasy_pred.info())
#print(df_vaxijen_pred.describe())

# print(df_atomic_pred.describe())

# pred_aromaticity = df_protparam_pred['aromaticity'].tolist()
# print(df_protparam_pred.describe())
# print(df_protparam_pred.info())

# df_protparam_iedb = pd.read_csv('analysis_iedb_results/protparam_analysis.csv', index_col=False)
# iedb_aromaticity = df_protparam_iedb['aromaticity'].tolist()
# ds1_sample = random.sample(iedb_aromaticity, len(pred_aromaticity))
# print(df_protparam_iedb.describe())
# print(df_protparam_iedb.info())

# corr = df_protparam_pred['instability_index'].corr(df_protparam_iedb['instability_index'])
# corr = df_protparam_pred.corrwith(df_protparam_iedb, axis=0)
#print(corr)

#Expasy
# df_expasy_pred = pd.read_csv('pred_analysis_results/expasy_analysis.csv', index_col=False)
# pred_aliphatic_index = df_expasy_pred['aliphatic_index'].tolist()
# print(df_expasy_pred.describe())
# print(df_expasy_pred.info())

df_protparam_iedb = pd.read_csv('analysis_iedb_results/protparam_analysis.csv', index_col=False)
df_expasy_iedb = pd.read_csv('analysis_iedb_results/expasy_analysis.csv', index_col=False)
#df_pepstats_iedb = pd.read_csv('analysis_iedb_results/pepstats_analysis.csv', index_col=False)
#df_atomic_iedb = pd.read_csv('analysis_iedb_results/atomic_composition.csv')
df_vaxijen_iedb = pd.read_csv('analysis_iedb_results/vaxijen_analysis.csv', index_col=False)
#print(df_aa_iedb.describe())
# print(df_expasy_iedb.describe())
# iedb_aliphatic_index= df_expasy_iedb['aliphatic_index'].tolist()

# print(df_aa_pred.mean())
# print(df_aa_iedb.mean())

# ds1_sample = random.sample(pred_aliphatic_index, len(iedb_aliphatic_index))
# print(df_expasy_iedb.describe())
# print(df_expasy_iedb.info())

# # corr = df_aromaticity_pred['instability_index'].corr(df_aromaticity_iedb['instability_index'])
# corr = df_expasy_pred.corrwith(df_expasy_iedb, axis=0)
# print(corr)

# summary, results = rp.ttest(group1= df_atomic_iedb[df_aa_iedb.columns[3]],
#                             group2= df_atomic_pred[df_atomic_pred.columns[3]])

# print(summary)
# print(results)
# print(stats.ttest_ind(df_atomic_iedb[df_atomic_iedb.columns[0]],
#                 df_atomic_pred[df_atomic_pred.columns[0]]))
# print("Negatively charged residues ", stats.ttest_ind(df_expasy_iedb[df_expasy_iedb.columns[1]],
#                 df_expasy_pred[df_expasy_pred.columns[1]]))
# print("Positively charged residues ", stats.ttest_ind(df_expasy_iedb[df_expasy_iedb.columns[2]],
#                 df_expasy_pred[df_expasy_pred.columns[2]]))
# print("Half-life ", stats.ttest_ind(df_expasy_iedb[df_expasy_iedb.columns[3]],
#                 df_expasy_pred[df_expasy_pred.columns[3]]))
# print("Aliphatic index ", stats.ttest_ind(df_expasy_iedb[df_expasy_iedb.columns[4]],
#                 df_expasy_pred[df_expasy_pred.columns[4]]))



# result = 1 - spatial.distance.cosine(ds1_sample, iedb_aliphatic_index)
# print(result)
pred_df = reduce(lambda x,y: pd.merge(x,y, on='peptide', how='outer'), [df_expasy_pred, df_protparam_pred, df_vaxijen_pred])
#print(pred_df)
iedb_df = reduce(lambda x,y: pd.merge(x,y, on='peptide', how='outer'), [df_expasy_iedb, df_protparam_iedb, df_vaxijen_iedb])
#iedb_df = pd.merge(df_protparam_iedb, df_expasy_iedb, df_aa_iedb, df_atomic_iedb, on='peptide')
#print(iedb_df)

full_iedb = iedb_df[~iedb_df.isnull().any(axis=1)]
full_iedb.drop_duplicates(keep='first', inplace=True)
# #print(full_iedb)

full_pred = pred_df[~pred_df.isnull().any(axis=1)]
full_pred.drop_duplicates(keep='first', inplace=True)
# #print(full_pred)

# ds1 = [5.5, 118.57, 1457.693, 5.24, 0.071, 36.8, 0.357, 0.214, 0.429, 1.021, -1.16]
# ds2 = [20.0, 21.67, 1823.828, 4.37, 0.056, 23.111, 0.111, 0.556, 0.056, -1.206, -1.236]

# result = 1 - spatial.distance.cosine(ds1, ds2)
# print(result)
# indexes = []
# for index, row in full_pred.iterrows():
#     #print(index)
#     pred_epitope_parameters = row.tolist()
#     #print(pred_epitope_parameters)
#     #print(pred_epitope_parameters)
#     count = 0
#     for iedb_index, iedb_row in full_iedb.iterrows():
#         iedb_epitope_parameters = iedb_row.tolist()
#         #print(iedb_epitope_parameters)
#         #print(len(pred_epitope_parameters), len(pred_epitope_parameters))
#         #print(pred_epitope_parameters, iedb_epitope_parameters)
#         result = 1 - spatial.distance.cosine(pred_epitope_parameters[1:-1], iedb_epitope_parameters[1:-1])
#         #print(result)
        
#         if result >= 0.99:
#             count += 1
#             #print(index, pred_epitope_parameters[0],iedb_epitope_parameters[0])
#             # indexes.append(index)
#     indexes.append([index, pred_epitope_parameters[0], count])
    
# for i in indexes:
#     print(i)
# with open('cosine_similarity.csv', 'w') as f: 
#     write = csv.writer(f) 
#     write.writerows(indexes) 

df_cos = pd.read_csv('cosine_similarity.csv', index_col=False)
print(df_cos.describe())
# # sampling_difference = df_expasy_iedb['aliphatic_index'].values - df_expasy_pred['aliphatic_index'].values

# # stats.shapiro(sampling_difference)

# # fig = plt.figure(figsize= (20, 10))
# # ax = fig.add_subplot(111)

# # normality_plot, stat = stats.probplot(sampling_difference, plot= plt, rvalue= True)
# # ax.set_title("Probability plot of sampling difference", fontsize= 20)
# # ax.set

# # plt.show()



# print(len(iedb_aromaticity))
# print(len(pred_aromaticity))

# data = np.array([ds1_sample, pred_aromaticity])
# covMatrix = np.cov(data, bias=False)
# sn.heatmap(covMatrix, annot=True, fmt='g')
# plt.show()

ds1 = [75, 0, 1, 23]
# ds1 = []
# for i in range(0,100):
#     n = random.randint(-500,500)
#     ds1.append(n)
# ds1_sample = random.sample(ds1, 50)

ds2 = [12, 43, 45, 6]
# ds2 = []
# for i in range(0,100):
#     n = random.randint(-500,500)
#     ds2.append(n)
# ds2_sample = random.sample(ds2, 50)

# data = np.array([ds1, ds2])
# covMatrix = np.cov(data, bias=False)
# sn.heatmap(covMatrix, annot=True, fmt='g')
# plt.show()
#print(covMatrix)



# x = [0, 1, 2, 3, 4, 5]
# y = [20, 15, 25, 0, 5, 10]
# correlation, p_value = stats.pearsonr(x, y)
# print(correlation)
