import pandas as pd
import numpy as np
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input_maf', type=str, default='data/test_input.maf', help='The input maf file')
parser.add_argument('--annotated_maf', type=str, default='data/test_annotated.maf', help='Annotated maf file')
parser.add_argument('--scripts_dir', type=str, default=os.path.join(os.getcwd()), help='scripts dir')


def get_gene_list(input_file):
	gene_new_names = {'MLL': 'KMT2A', "MLL3": "KMT2C", "MLL2": "KMT2D", "MLL4": "KMT2B"}
	gene_list_cv = pd.read_table(input_file, sep='\t', header=None)
	gene_list_cv.columns = ['region', 'cytoband', 'gene_details']
	gene_list_cv['gene_name'] = gene_list_cv['gene_details'].str.split(":").str.get(0)
	gene_list_cv = gene_list_cv.replace({'gene_name': gene_new_names})
	gene_list_cv = gene_list_cv.dropna(axis=0)
	cv_gene_list = gene_list_cv['gene_name'].unique().tolist()
	#print(len(cv_gene_list))
	return(cv_gene_list)

def get_myvariant_annotations(inputmaf, input_df):
	print "\n\nrun myvariant info to get annotations"
	#get relevant columns from maf file to annotate for myvariant info
	#Chromosome\tStart_Position\tReference_Allele\tAlternate_Allele\tEnd_Position\tVariant_Type
	variant_info_input_file = inputmaf.replace(".maf",".input_to_variantinfo.maf")
	columns_to_keep = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Alternate_Allele', 'End_Position', 'Variant_Type']
	input_df_subset = input_df[columns_to_keep]
	input_df_subset = input_df_subset.drop_duplicates()
	input_df_subset.to_csv(variant_info_input_file, sep = '\t', index = None)

	#we want to sort without the header, so save header to a different file and sort without header, and finally merge files and delete tmp file
	#faster awk solution but discontinued due to possible reordering of input files
	#variant_info_file_tmp = inputmaf.replace(".maf",".input_to_variantinfo.tmp.maf")
	# str_awk = '''awk  -F "\\t" 'BEGIN {OFS = "\\t"} {print $5, $6, $11, $12, $7, $10 }' '''+ inputmaf +" | tail -n+1 | sort | uniq  > "+variant_info_file_tmp
	# os.system(str_awk)
	# print str_awk
	# str_header = "cat " +scripts_dir+ "myvariantinfo_header.txt " +variant_info_file_tmp+" > "+variant_info_input_file
	# str_rm_tmp = "rm -rf "+variant_info_file_tmp
	# os.system(str_header)
	# os.system(str_rm_tmp)
	#script below calls wrapper around myvariant.info and performs batch-wise annotation. 
	myvariant_path = os.path.join(scripts_dir, "run_variant_info_annotation_biothing.py")
	str_myvariant = "python "+ myvariant_path +" "+variant_info_input_file+" "+myvariant_batch_location
	os.system(str_myvariant)
	print str_myvariant

def merge_myvariant_annotations():
	#read and merge variant info batches
	print "\n\nmerge myvariant info annotations to "+variant_info_output_merge
	files = os.listdir(myvariant_batch_location)
	keep_columns = ['query', 'cadd.encode.h3k4me1', 'dbnsfp.fathmm-mkl.coding_rankscore', 
                         'dbnsfp.mutationassessor.rankscore', 'cadd.encode.h3k27ac',  'dbnsfp.eigen-pc.raw_coding',
                         'cadd.phast_cons.primate', 'dbnsfp.genocanyon.score', 'cadd.encode.exp']

	if os.path.isfile(variant_info_output_merge):
    		print("delete existing variant info merge file")
       		str_rm = "rm -f "+variant_info_output_merge
    		os.system(str_rm)

	for filename in files:
	    batch_input = os.path.join(myvariant_batch_location, filename)
	    batch_data = pd.read_table(batch_input, sep='\t')
	    for column in keep_columns:
	        if(column not in batch_data.columns.tolist()):
	            batch_data[column] = np.NaN
	            print ("\n\nWARNING : myvariant.info batch file "+batch_input+" does not contain the column "+column)
	    batch_data = batch_data[keep_columns]
	    if not os.path.isfile(variant_info_output_merge):
	        batch_data.to_csv(variant_info_output_merge, header =True, index=None, sep='\t')
	    else: # else it exists so append without writing the header
	        batch_data.to_csv(variant_info_output_merge, mode = 'a',header=False, index=None, sep='\t')


def get_oncokb_annotations(inputmaf, input_df):
	print "\n\nrun OncoKB get annotations"
	oncokb_input_file = inputmaf.replace(".maf",".oncokb_input.maf")
	oncokb_output_file = inputmaf.replace(".maf",".oncokb_output.maf")
	
	columns_to_keep = ['Hugo_Symbol','Variant_Classification', 'Protein_position', 'HGVSp_Short']
	input_df_subset = input_df[columns_to_keep]
	input_df_subset = input_df_subset.drop_duplicates()
	input_df_subset['Normal_Sample'] = "SAMPLE"
	input_df_subset.to_csv(oncokb_input_file, sep = '\t', index = None)


	#oncokb_input_file_tmp = inputmaf.replace(".maf",".oncokb_input.tmp.maf")
	# str_awk = '''awk  -F "\\t" 'BEGIN {OFS = "\\t"} {print $1, $9, $32, $17, "Normal_Sample"}' ''' +inputmaf+" | sort | uniq  > "+oncokb_input_file_tmp
	# str_header = "cat " + scripts_dir + "oncokb_header.txt " +oncokb_input_file_tmp+" > "+oncokb_input_file
	cwd = os.getcwd()
	oncokb_cmd = "python MafAnnotator.py -a -i "+cwd+"/"+oncokb_input_file+" -o "+cwd+"/"+oncokb_output_file
	#oncoKB needs to be run from it's own directory so change to that dir and then back
	os.chdir(oncokb_path)
	os.system(oncokb_cmd)
	os.chdir(cwd)
	return oncokb_output_file

def join_annotations(oncokb_file, variantinfo_file, cohort_maf):
	print "\n\njoin different annotation sources"
	oncoKB_maf = pd.read_table(oncokb_file, sep = '\t')
	print oncoKB_maf.shape
	variantinfo_maf = pd.read_table(variantinfo_file, sep = '\t')
	print variantinfo_maf.shape

	#merge with OncoKB
	oncoKB_maf['mutation']= oncoKB_maf['Hugo_Symbol']+":"+oncoKB_maf['HGVSp_Short']+":"+oncoKB_maf['Protein_position'].astype(str)
	cohort_maf['mutation']= cohort_maf['Hugo_Symbol']+":"+cohort_maf['HGVSp_Short']+":"+cohort_maf['Protein_position'].astype(str)
	oncoKB_maf_subset = oncoKB_maf[['mutation', 'oncogenic', 'is-a-hotspot', 'is-a-3d-hotspot']]
	oncoKB_maf_subset = oncoKB_maf_subset[~(pd.isnull(oncoKB_maf_subset['mutation']))]
	cohort_maf = pd.merge(cohort_maf, oncoKB_maf_subset, on = 'mutation', how = 'left')
	cohort_maf['oncogenic'] = np.where(pd.isnull(cohort_maf['mutation']), 0, cohort_maf['oncogenic'])
	cohort_maf['is-a-hotspot'] = np.where(pd.isnull(cohort_maf['mutation']), 0, cohort_maf['is-a-hotspot'])
	cohort_maf['is-a-3d-hotspot'] = np.where(pd.isnull(cohort_maf['mutation']), 0, cohort_maf['is-a-3d-hotspot'])

	cohort_maf = cohort_maf.drop_duplicates()
	cohort_maf['Chromosome'] = cohort_maf['Chromosome'].astype(str)
	#merge with myvariant.info
	cohort_maf['query'] = np.where(cohort_maf['Variant_Type']=='SNP',
                              "chr"+cohort_maf['Chromosome']+":g."+cohort_maf['Start_Position'].astype(str)+cohort_maf['Reference_Allele']+">"+cohort_maf['Alternate_Allele'],
                              "chr"+cohort_maf['Chromosome']+":g."+cohort_maf['Start_Position'].astype(str)+"_"+cohort_maf['End_Position'].astype(str)+"ins"+cohort_maf['Alternate_Allele'] #INS
                              )

	cohort_maf['query'] = np.where(cohort_maf['Variant_Type']=='DEL',
                              "chr"+cohort_maf['Chromosome']+":g."+cohort_maf['Start_Position'].astype(str)+"_"+cohort_maf['End_Position'].astype(str)+"del", #DEL,
                              cohort_maf['query'] 
                              )
	cohort_maf = pd.merge(cohort_maf, variantinfo_maf, on=['query'], how='left')
	print cohort_maf.shape

	#merge with dbScSNV
	#dbscsnv annotations
	print ("perform dbScSNV annotations. This will take a while")
	#dbscSNV = pd.read_csv(dbscSNV_file, compression='gzip',sep = '\t' )

	dbscsnv_files = os.listdir(dbscSNV_folder)
	column_names = ['chr', 'pos', 'ref', 'alt', 'ada_score', 'rf_score']
	dbscSNV = pd.DataFrame(columns = column_names )
	for filename in dbscsnv_files:
	    dbscsnv_file = os.path.join(dbscSNV_folder, filename)
	    if(filename.startswith("dbscSNV1.1.chr")):
	        #print(dbscsnv_file)
	        dbscSNV_1 = pd.read_table(dbscsnv_file, sep='\t', compression = 'gzip')
	        dbscSNV_1 = dbscSNV_1[column_names]
	        dbscSNV = dbscSNV.append(dbscSNV_1, ignore_index=True)
	    
	print(dbscSNV.shape)
	print dbscSNV.head()
	dbscSNV.columns = ['Chromosome','Start_Position', 'Reference_Allele',  'Alternate_Allele', 'ada_score', 'rf_score']
	dbscSNV['Chromosome'] = dbscSNV['Chromosome'].astype(str)
	dbscSNV['Start_Position'] = dbscSNV['Start_Position'].astype(int)
	    
	# print(dbscSNV.shape)
	print dbscSNV.head()
	dbscSNV['Chromosome'] = dbscSNV['Chromosome'].astype(str)
	dbscSNV['Start_Position'] = dbscSNV['Start_Position'].astype(int)
	cohort_maf = pd.merge(cohort_maf, dbscSNV, 
	                 on=['Chromosome', 'Start_Position', 'Reference_Allele',  'Alternate_Allele'], how='left')
		
	return cohort_maf



def Calculate_MAF(cohort_maf):
	print "\n\ncalculate Minor allele frequency"
	sample_list = cohort_maf[['Normal_Sample']].drop_duplicates()
	sample_list['Panel_type'] = sample_list['Normal_Sample'].str.split("-").str.get(-1)

	sample_count_panel = pd.DataFrame({
	        'sample_count_panel' : 
	        sample_list.groupby(['Panel_type'])['Normal_Sample'].nunique()
	    }).reset_index()

	dict_panel_count = dict(zip(sample_count_panel['Panel_type'], sample_count_panel['sample_count_panel']))
	expected_keys = ['IM3', 'IM5', 'IM6']
	for key in expected_keys:
		if key not in dict_panel_count.keys():
			dict_panel_count[key] = 0

	print(dict_panel_count)

	cohort_maf['n_alt_freq'] = cohort_maf['n_alt_count']/cohort_maf['n_depth']
	cohort_maf["allele_count"] =  np.where(cohort_maf['n_alt_freq']>0.75, 2, 1)
	cohort_maf['Hugo_Symbol'] = cohort_maf['Hugo_Symbol'].fillna(".")

	MAF_alterations = pd.DataFrame({
	        'minor_allele_count' : 
	        cohort_maf.groupby(['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Alternate_Allele',
	                              'Hugo_Symbol'
	                            ])['allele_count'].sum(),
	        'median_VAF' : 
	        cohort_maf.groupby(['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Alternate_Allele',
	                              'Hugo_Symbol'
	                            ])['n_alt_freq'].median()
	        
	    }).reset_index()
	print(MAF_alterations.shape)
	MAF_alterations['sample_count_bypanel'] = dict_panel_count['IM3']+dict_panel_count['IM5']+dict_panel_count['IM6']

	#adjust MAF calculation for genes that are only in CV5 or CV6
	MAF_alterations['sample_count_bypanel'] = np.where(MAF_alterations['Hugo_Symbol'].isin(cv6_only_list), 
	                                                   dict_panel_count['IM6'], 
	                                                   MAF_alterations['sample_count_bypanel']
	                                                  )
	MAF_alterations['sample_count_bypanel'] = np.where(MAF_alterations['Hugo_Symbol'].isin(cv5_only_list), 
	                                                   dict_panel_count['IM6'] + dict_panel_count['IM5'], 
	                                                   MAF_alterations['sample_count_bypanel']
	                                                  )
	MAF_alterations['sample_count_bypanel'] = MAF_alterations['sample_count_bypanel'].astype(float)
	MAF_alterations['MinorAlleleFreq'] = MAF_alterations['minor_allele_count']/(2* MAF_alterations['sample_count_bypanel'])
	MAF_alterations = MAF_alterations.sort_values('MinorAlleleFreq', ascending=False)
	MAF_alterations = MAF_alterations[['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Alternate_Allele',
	                                 'median_VAF', 'MinorAlleleFreq' ]]
	return MAF_alterations


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

def downstream_annotations(cohort_maf_uniq, oncogenes, tumor_suppressors):
	try:
	    cohort_maf_uniq['ExAC2_AF_ASJ'] = np.where(cohort_maf_uniq['ExAC2_AF_ASJ']=='.',0, 
	                                               cohort_maf_uniq['ExAC2_AF_ASJ'])
	except TypeError:
	    pass

	cohort_maf_uniq['ExAC2_AF_ASJ'] = cohort_maf_uniq['ExAC2_AF_ASJ'].fillna(0)
	cohort_maf_uniq['ExAC2_AF_ASJ'] = cohort_maf_uniq['ExAC2_AF_ASJ'].astype(float)
	cohort_maf_uniq['ExAC2_AF'] = np.where(pd.isnull(cohort_maf_uniq['ExAC2_AF']), 0, cohort_maf_uniq['ExAC2_AF'])
	cohort_maf_uniq['mutation_mechanism_consistency'] = 0

	cohort_maf_uniq['mutation_mechanism_consistency'] = np.where((cohort_maf_uniq['Variant_Classification'].isin([
	    'Nonsense_Mutation', 'Frame_Shift_Ins', 'Frame_Shift_del', 'Splice_Site', 'Missense_Mutation', 
	    'Splice_Region'
	        ])) &
	      (cohort_maf_uniq['Hugo_Symbol'].isin(tumor_suppressors)), 1,  
	      cohort_maf_uniq['mutation_mechanism_consistency'] ) 

	cohort_maf_uniq['mutation_mechanism_consistency'] = np.where((cohort_maf_uniq['Variant_Classification'].isin([
	    'Missense_Mutation'])) & (cohort_maf_uniq['Hugo_Symbol'].isin(oncogenes)), 1,  
	      cohort_maf_uniq['mutation_mechanism_consistency'] ) 

	list_consequences =  cohort_maf_uniq['Consequence'].unique().tolist()
	list_consequences_uniq = []
	for i in list_consequences:
	    list_items = i.split(",")
	    list_consequences_uniq = list_consequences_uniq + list_items
	list_consequences_uniq = list(set(list_consequences_uniq))
	list_consequences_uniq = [i for i in list_consequences_uniq]
	
	for colname in list_consequences_uniq:
	    cohort_maf_uniq['Consequence_'+colname] = np.where(
	        cohort_maf_uniq['Consequence'].str.find(colname)>-1, 1, 0)
	list_consequences_uniq = ["Consequence_"+i for i in list_consequences_uniq]  




	splice_variants = cohort_maf_uniq[cohort_maf_uniq['Variant_Classification'].isin(['Splice_Site','Splice_Region'])]
	non_splice_variants = cohort_maf_uniq[~(cohort_maf_uniq['Variant_Classification'].isin(['Splice_Site','Splice_Region']))]
	non_splice_variants['splice_dist'] = -200

	splice_variants['temp'] = splice_variants['HGVSc']
	splice_variants['temp'] = splice_variants['temp'].str.replace("+",";")
	splice_variants['temp'] = splice_variants['temp'].str.replace("-",";")

	splice_variants['temp1'] = splice_variants['temp'].str.split(";").str.get(1)
	splice_variants['splice_dist'] = splice_variants['temp1'].str.extract('(\d+)')
	splice_variants['splice_dist'] = splice_variants['splice_dist'].fillna(0)

	splice_variants = splice_variants.drop(['temp', 'temp1'], axis =1 )

	cohort_maf_uniq = pd.concat([splice_variants, non_splice_variants])
	print cohort_maf_uniq.shape

	cohort_maf_uniq['splice_dist'] = cohort_maf_uniq['splice_dist'].astype(int)

	cohort_maf_uniq['mutation'] = cohort_maf_uniq['Hugo_Symbol'] + ":" + cohort_maf_uniq['HGVSc']
	cohort_maf_uniq['ratio_ASJ'] = cohort_maf_uniq['ExAC2_AF_ASJ']/cohort_maf_uniq['ExAC2_AF'].astype(float)
	cohort_maf_uniq['ratio_ASJ'] = cohort_maf_uniq['ratio_ASJ'].astype(float)
	cohort_maf_uniq['ratio_ASJ'] = cohort_maf_uniq['ratio_ASJ'].fillna(0.0)
	cohort_maf_uniq['ExAC2_AF'] = cohort_maf_uniq['ExAC2_AF'].fillna(0.0)

	pathogenic_terms = ['Pathogenic', 'pathogenic', ]
	benign_terms = ['Benign', 'benign']
	uncertain_terms = ['conflicting'] 
	cohort_maf_uniq['clinvar_pathogenic'] = np.where(cohort_maf_uniq['CLINICAL_SIGNIFICANCE'].str.contains('|'.join(pathogenic_terms),na=False),
	                                          1, 0)
	cohort_maf_uniq['clinvar_uncertain'] = np.where(cohort_maf_uniq['CLINICAL_SIGNIFICANCE'].str.contains('|'.join(uncertain_terms), na=False),
	                                          1, 0)
	cohort_maf_uniq['clinvar_benign'] = np.where(cohort_maf_uniq['CLINICAL_SIGNIFICANCE'].str.contains('|'.join(benign_terms), na=False),
	                                          1, 0)

	cohort_maf_uniq['is-a-hotspot'] = cohort_maf_uniq['is-a-hotspot'].fillna('N')
	try:
		cohort_maf_uniq['is-a-hotspot'] = np.where(cohort_maf_uniq['is-a-hotspot']=='Y', 1 , 0)
	except TypeError:
		cohort_maf_uniq['is-a-hotspot'] = 0

	cohort_maf_uniq['oncogenic'] = np.where(cohort_maf_uniq['oncogenic'].isin(['Oncogenic', 'Likely Oncogenic']), 1, 0)

	#annotate last exon - 5' end of exon 
	cohort_maf_uniq['exon_number'] = cohort_maf_uniq['Exon_Number'].str.split('/').str.get(0)
	cohort_maf_uniq['total_exons'] = cohort_maf_uniq['Exon_Number'].str.split('/').str.get(1)
	cohort_maf_uniq['Protein_position'] = cohort_maf_uniq['Protein_position'].str.replace("\?-", "")

	cohort_maf_uniq['protein_position'] = cohort_maf_uniq['Protein_position'].str.split('/').str.get(0)
	cohort_maf_uniq['transcript_len'] = cohort_maf_uniq['Protein_position'].str.split('/').str.get(1)

	cohort_maf_uniq['protein_position'] = np.where(cohort_maf_uniq['protein_position']=='-', 0, cohort_maf_uniq['protein_position'])
	cohort_maf_uniq['protein_position'] = np.where(cohort_maf_uniq['protein_position'].str.find('?')>-1,
	                                             cohort_maf_uniq['protein_position'].str.replace('?', ''), 
	                                                                 cohort_maf_uniq['protein_position'])
	cohort_maf_uniq['protein_position'] = np.where(cohort_maf_uniq['protein_position'].str.find('-')>-1,
	                                             cohort_maf_uniq['protein_position'].str.split("-").str.get(0), 
	                                                                 cohort_maf_uniq['protein_position'])


	cohort_maf_uniq['protein_position'] = cohort_maf_uniq['protein_position'].fillna(0)
	cohort_maf_uniq['transcript_len'] = cohort_maf_uniq['transcript_len'].fillna(0)
	cohort_maf_uniq['protein_position'] = cohort_maf_uniq['protein_position'].astype(int)
	cohort_maf_uniq['transcript_len'] = cohort_maf_uniq['transcript_len'].astype(int)
	cohort_maf_uniq['tail_end'] = np.where((cohort_maf_uniq['protein_position']/cohort_maf_uniq['transcript_len'])>0.95, 1, 0)
	cohort_maf_uniq['last_exon_terminal'] = np.where((cohort_maf_uniq['tail_end']==1) & 
	                                                                   (cohort_maf_uniq['exon_number']==cohort_maf_uniq['total_exons']),
	                                               1, 0)
	return cohort_maf_uniq


def main():
	global args, cv5_only_list, cv6_only_list, myvariant_batch_location, variant_info_output_merge, scripts_dir
	global annotation_path, oncokb_path, gene_level_input_file, gene_function_map_input, dbscSNV_folder
	args = parser.parse_args()
	input_maf = args.input_maf
	annotated_maf = args.annotated_maf
	scripts_dir = args.scripts_dir

	#this is the path where all the scripts and files exist. 
	annotation_path = os.path.join(scripts_dir, "annotation_files")

	#paths to annotation sources
	oncokb_path = os.path.join(scripts_dir, "oncokb-annotator")
	gene_level_input_file = os.path.join(annotation_path, "gene_level_annotation.txt")
	gene_function_map_input = os.path.join(annotation_path, "gene_function_other_genes.txt")
	dbscSNV_folder = os.path.join(annotation_path, "dbscSNV")


	#this directory gets created at location where input file resides
	full_input_path = os.path.abspath(input_maf)
	working_dir = os.path.dirname(full_input_path)
	print ("all intermediate files will get created in this directory : "+working_dir)
	
	myvariant_batch_location = os.path.join(working_dir, "variant_info_batches")
	variant_info_output_merge = os.path.join(working_dir, "merged_variant_info_batches.txt")

	#check if variant_info_batch folder exists, if not, create it
	if(os.path.isdir(myvariant_batch_location)):
		pass
	else:
		print ("creating directory for writing myvariant info batch files")
		create_dir = "mkdir "+myvariant_batch_location
		os.system(create_dir)

	gene_list_cv5_file = os.path.join(annotation_path, "IMPACT_cv5.gene_intervals.list.annotated")
	cv5_gene_list = get_gene_list(gene_list_cv5_file)

	gene_list_cv3_file = os.path.join(annotation_path, "IMPACT_cv3.gene_intervals.list.annotated")
	cv3_gene_list = get_gene_list(gene_list_cv3_file)

	gene_list_cv6_file = os.path.join(annotation_path, "IMPACT_cv6.gene_intervals.list.annotated")
	cv6_gene_list = get_gene_list(gene_list_cv6_file)

	#get list of genes that are only in CV5 and not in CV3
	cv5_only_list = set(cv5_gene_list).difference(set(cv3_gene_list))
	cv6_only_list = set(cv6_gene_list).difference(set(cv5_gene_list))

	input_maf_data = pd.read_table(input_maf, sep = '\t')
	print input_maf_data.shape


	get_myvariant_annotations(input_maf, input_maf_data)
	oncokb_maf  = get_oncokb_annotations(input_maf, input_maf_data)
	merge_myvariant_annotations()

	merged_maf = join_annotations(oncokb_maf, variant_info_output_merge, input_maf_data)
	merged_maf = merged_maf.drop_duplicates()
	merged_maf.to_csv(annotated_maf.replace(".maf", ".extended.maf"), sep = '\t', index = None)

	print merged_maf.shape
	#remove black list variants
	try:
		merged_maf = merged_maf[~(merged_maf['DMP_Blacklist']=='Blacklist')]
	except TypeError:
		pass
	except KeyError:
		pass
	print merged_maf.shape

	#filter variant types to include only those in training data -- this step has been removed 
	variant_types = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins',
	                 'Intron', 'Missense_Mutation', 'Nonsense_Mutation', 'Splice_Region', 
	                 'Splice_Site', 'Translation_Start_Site']

	#merged_maf = merged_maf[merged_maf['Variant_Classification'].isin(variant_types)]

	#get number of patients by different panel types

	MAF_alterations_cohort = Calculate_MAF(merged_maf)

	merged_maf = pd.merge(merged_maf, MAF_alterations_cohort, 
	                      on=['Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Alternate_Allele',],
	                       how='left')

	#merged_maf_rare = merged_maf[(merged_maf['MinorAlleleFreq']<0.02) ]
	merged_maf_rare = merged_maf[merged_maf['ExAC2_AF']<0.02 ]
	merged_maf_uniq = merged_maf_rare.drop_duplicates(subset=['Chromosome', 'Start_Position',
	                                                                 'End_Position', 'Reference_Allele',
	                                                                 'Alternate_Allele'])
	print merged_maf_uniq.shape

	merged_maf_uniq.to_csv(annotated_maf.replace(".maf", ".uniq.maf"), sep = '\t', index = None)

	oncogenes, tumor_suppressors, gene_level_annotation = get_gene_level_annotations()
	merged_maf_uniq = pd.merge(merged_maf_uniq, gene_level_annotation, on='Hugo_Symbol', how='left')

	final_maf = downstream_annotations(merged_maf_uniq, oncogenes, tumor_suppressors)

	final_maf.to_csv(annotated_maf, sep = '\t', index = None)

if __name__ == '__main__':
	main()	
