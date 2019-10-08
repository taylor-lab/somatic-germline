import pandas as pd
import numpy as np
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn import preprocessing, cross_validation
from sklearn.metrics import confusion_matrix
from sklearn.metrics import average_precision_score
from sklearn.cross_validation import KFold, cross_val_score
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--classifier_input', type=str, default='data/test_annotated.maf', help='The annotated unique input to classifier')
parser.add_argument('--classifier_output', type=str, default='data/test.classifier_output.maf', help='classifier output')
parser.add_argument('--scripts_dir', type=str, default=os.path.join(os.getcwd()), help='scripts dir')

def get_df_prob( predictions):
  df_predictions = pd.DataFrame(data=predictions, columns=['Prob_0', 'Prob_1'])
  df_predictions['Predicted_label'] = np.where(df_predictions['Prob_1']> df_predictions['Prob_0'], 1, 0)
  df_predictions['Predicted_label_soft'] = np.where(df_predictions['Prob_1']> 0.45, 1, 0)
  df_predictions['prediction_prob'] = df_predictions['Prob_1']
  return(df_predictions)

def get_gene_level_annotations():
  #get annotations for germline genes as oncogenes / tumor suppressors or genes whose function is via gain vs. loss
  gene_level_annotation = pd.read_table(gene_level_input_file, sep = '\t')
  gene_annotation_columns = gene_level_annotation.columns.tolist()
  gene_annotation_columns = gene_annotation_columns.remove('Hugo_Symbol')

  oncogenes = gene_level_annotation[gene_level_annotation['OncoKB Oncogene']==1]['Hugo_Symbol'].unique().tolist()
  oncogenes.append('TP53')
  #POLE and POLD1 are added as oncogenes because we know that only missense mutations in POLE lead to signature 10. 
  oncogenes.append('POLE')
  oncogenes.append('POLD1')
  non_cancergenes = gene_level_annotation[gene_level_annotation['OMIM']==0]['Hugo_Symbol'].unique().tolist()


  tumor_suppressors =  gene_level_annotation[gene_level_annotation['OncoKB TSG']==1]['Hugo_Symbol'].unique().tolist()
  tumor_suppressors.remove('POLE')
  tumor_suppressors.remove('POLD1')
  tumor_suppressors.append('EPCAM')

  #this file contains some additional gene annotations from CB.
  function_map_other_genes = pd.read_table(gene_function_map_input, sep='\t')
  other_gof_genes = function_map_other_genes[function_map_other_genes['oncogenic mechanism']=='gain-of-function']['Hugo_Symbol'].tolist()
  other_lof_genes = function_map_other_genes[function_map_other_genes['oncogenic mechanism']=='loss-of-function']['Hugo_Symbol'].tolist()
  #print gene_level_annotation.head()

  tumor_suppressors =list(set(tumor_suppressors + other_lof_genes))
  oncogenes =list(set(oncogenes + other_gof_genes+['VTCN1', 'YES1', 'XPO1', 'TRAF7']))
  print(sorted(oncogenes[:5]))
  print(sorted(tumor_suppressors[:5]))
  return oncogenes, tumor_suppressors, gene_level_annotation



def main():
  args = parser.parse_args()
  classifier_input = args.classifier_input
  classifier_output = args.classifier_output
  scripts_dir = args.scripts_dir
  global gene_level_input_file, gene_function_map_input
  #this is the path where all the scripts and files exist. 
  annotation_path = os.path.join(scripts_dir, "annotation_files")

  #paths to annotation sources
  gene_level_input_file = os.path.join(annotation_path, "gene_level_annotation.txt")
  gene_function_map_input = os.path.join(annotation_path, "gene_function_other_genes.txt")
  germline_gene_list = os.path.join(annotation_path, "germline_gene_list.txt")
  training_data_file = os.path.join(scripts_dir, "data/classifier_training_data_V1.txt.gz") 

  rare_variants_76genes = pd.read_table(training_data_file, sep="\t", compression = 'gzip')
  labels = rare_variants_76genes['signed_out'].tolist()

  keep_columns = ['ExAC2_AF', 'ExAC2_AF_ASJ', 'clinvar_pathogenic','clinvar_benign', 'clinvar_uncertain','GOLD_STARS', 
  						'Consequence_frameshift_variant','Consequence_splice_region_variant', 'Consequence_missense_variant', 'Consequence_stop_gained', 
                            'Consequence_splice_acceptor_variant', 'Variant_Classification', 
                            'dbnsfp.fathmm.mkl.coding_rankscore', 'dbnsfp.mutationassessor.rankscore',  'dbnsfp.eigen.pc.raw',
                            'cadd.phast_cons.primate', 'dbnsfp.genocanyon.score',  'MinorAlleleFreq' ,'oncogenic', 'is-a-hotspot',
                           'ada_score', 'last_exon_terminal', 'OMIM', 'cell cycle checkpoint', 'HR', 'DSB', 'DNA replication', 'DSBR', 'MMR', 'DNA repair', 
                                'cell cycle', 'response to DNA damage stimulus', 
                                'BER', 'NER', 'OncoKB Oncogene', 'OncoKB TSG','mutation_mechanism_consistency',
                                 'ratio_ASJ', 'splice_dist', 'Consequence_splice_donor_variant'
                           ]

  final_columns = ['ExAC2_AF', 'clinvar_pathogenic', 'Consequence_frameshift_variant',
                            'ExAC2_AF_ASJ',  'Consequence_splice_region_variant', 
                              'Consequence_missense_variant', 
                              'Consequence_stop_gained', 
                            'Consequence_splice_acceptor_variant', 'Consequence_splice_donor_variant' ,  
                           'GOLD_STARS', 'Variant_Classification_Nonsense_Mutation',
                           'clinvar_benign', 'clinvar_uncertain' , 'dbnsfp.fathmm.mkl.coding_rankscore', 
                           'dbnsfp.mutationassessor.rankscore',  'dbnsfp.eigen.pc.raw',
                            'cadd.phast_cons.primate', 
                           'dbnsfp.genocanyon.score','Variant_Classification_Intron', 
                           'Variant_Classification_Splice_Region','Variant_Classification_Splice_Site',
                           #'Variant_Classification_Nonstop_Mutation', 'Variant_Classification_RNA',
                             'MinorAlleleFreq' ,'oncogenic', 'is-a-hotspot',
                           'ada_score', 'last_exon_terminal',
                               'OMIM', 'cell cycle checkpoint', 'HR', 'DSB', 'DNA replication', 'DSBR', 'MMR', 'DNA repair', 
                                'cell cycle', 'response to DNA damage stimulus', 
                                'BER', 'NER', 'OncoKB Oncogene', 'OncoKB TSG','mutation_mechanism_consistency',
                                 'ratio_ASJ', 'splice_dist', 
                               ]

  df_76genes = rare_variants_76genes[keep_columns]
  df_76genes = pd.get_dummies(df_76genes)
  print(df_76genes.shape)
  df_76genes = df_76genes.fillna(value=0)

  # final_columns = keep_columns.remove('Variant_Classification') + ['Variant_Classification_Nonsense_Mutation', 'Variant_Classification_Intron', 
  #                          'Variant_Classification_Splice_Region','Variant_Classification_Splice_Site',] 
  #                          #"Variant_Classification_3'UTR", "Variant_Classification_5'UTR",'Variant_Classification_Nonstop_Mutation', 'Variant_Classification_RNA'

  df_76genes = df_76genes[final_columns]

  #define model and train
  rf_weights = {0: 1, 1: 4.5}
  rf = RandomForestClassifier(n_estimators=300, criterion= 'gini', min_samples_leaf= 1,
                                   min_samples_split= 2, max_features=12, random_state=4, class_weight=rf_weights
                               )
  rf.fit(df_76genes, labels)
  k_fold = KFold(len(labels), n_folds=10)

  cv_score = cross_val_score(rf, df_76genes, labels, cv=k_fold, scoring='precision')
  print("\nPrecision: %0.2f (+/- %0.2f)" % (cv_score.mean(), cv_score.std() * 2))

  cv_score = cross_val_score(rf, df_76genes, labels, cv=k_fold, scoring='recall')
  print("\nRecall: %0.2f (+/- %0.2f)" % (cv_score.mean(), cv_score.std() * 2))
  coeffs = np.array(rf.feature_importances_)

  #get feature importances
  coeffs = np.array(rf.feature_importances_)
  feature_labels = df_76genes.columns.tolist()
  print("\nRandom Forest Feature importances")
  dict_LS_coeffs = {}
  for i in range(len(feature_labels)):
          dict_LS_coeffs[feature_labels[i]] = float(coeffs[i])
      
      
  # for w in sorted(dict_LS_coeffs, key=dict_LS_coeffs.get, reverse=True):
  #     print(w+"\t"+str(dict_LS_coeffs[w]))
      
  list_keepcols = []
  for w in sorted(dict_LS_coeffs, key=dict_LS_coeffs.get, reverse=True):
      if(dict_LS_coeffs[w]>0.01):
          list_keepcols.append(w)
  print("\n\n top features")
  print(list_keepcols)

  ###read data
  df_all = pd.read_table(classifier_input, sep = '\t')
  print(df_all.shape)
  df_all = df_all.rename(columns={'dbnsfp.fathmm-mkl.coding_rankscore': 'dbnsfp.fathmm.mkl.coding_rankscore', 
                                'dbnsfp.eigen-pc.raw': 'dbnsfp.eigen.pc.raw'})

  variant_types = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins',
                   'Intron', 'Missense_Mutation', 'Nonsense_Mutation', 'Splice_Region', 
                   'Splice_Site', 'Translation_Start_Site']
  df_all = df_all[df_all['Variant_Classification'].isin(variant_types)] #remove in V2
  df_all_subset = df_all[keep_columns]

  df_all_subset = pd.get_dummies(df_all_subset)
  df_all_subset = df_all_subset.fillna(value=0)
  try:
    df_all_subset = df_all_subset[final_columns]
  except KeyError:
    for colname in final_columns:
      if(colname not in df_all_subset.columns.tolist()):
        df_all_subset[colname] = 0
    df_all_subset = df_all_subset[final_columns]


  df_predictions_prob = rf.predict_proba(df_all_subset)
  df_predictions_prob = get_df_prob(df_predictions_prob)
  df_predictions = rf.predict(df_all_subset)
  e = pd.Series(df_predictions)
  f = pd.Series(df_predictions_prob['prediction_prob'])
  df_all['prediction'] = e.values
  df_all['prediction_probability'] = f.values
  
  oncogenes, tumor_suppressors, gene_level_annotation = get_gene_level_annotations()
  #identify truncating mutations in oncogene
  df_all['truncating_oncogene'] = np.where((df_all['Hugo_Symbol'].isin(oncogenes))&
                                                (~(df_all['Hugo_Symbol'].isin(tumor_suppressors))) &
                                      (df_all['Variant_Classification'].isin([
                  'Nonsense_Mutation', 'Frame_Shift_Ins', 'Frame_Shift_Del', 'Splice_Site'
              ])), 1, 0)

  df_all['prediction'] = np.where(df_all['truncating_oncogene']==1,
                                           0, df_all['prediction'])

  final_columns_to_report = ['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Alternate_Allele', 'Variant_Type', 'Variant_Classification', 'Hugo_Symbol', 'HGVSc', 'Protein_position', 'HGVSp_Short', 'Normal_Sample', 'n_alt_count', 'n_depth', 'ExAC2_AF', 'ExAC2_AF_ASJ', 'Consequence', 'CLINICAL_SIGNIFICANCE', 'GOLD_STARS', 'Exon_Number', 'oncogenic', 'last_exon_terminal', 'prediction', 'prediction_probability', 'truncating_oncogene',]

  df_all[final_columns_to_report].to_csv(classifier_output, sep = '\t', index = None)

if __name__ == '__main__':
  main()  

