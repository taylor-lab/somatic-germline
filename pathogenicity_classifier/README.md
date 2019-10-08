# Germline-Classifier
This is the source code for the random forest based classifier described in Srinivasan P, Bandlamudi C, et al. In preparation. This was used to classify rare germline variants called across 17,152 advanced cancer patients sequenced using the MSK-IMPACT assay.

Requirements:
Refer to requirements.txt to find other python dependencies including myvariant python package from https://docs.myvariant.info/en/latest/doc/packages.html#myvariant-python-module

Inputs to Scripts :
To run the classifier, you will need an input file in maf format with the following columns :

#Variant level info

Chromosome, Start_Position,End_Position, Reference_Allele, Alternate_Allele, Variant_Type, Variant_Classification,Consequence, 
Hugo_Symbol, Protein_position, HGVSp_Short, Exon_Number

#Patient level info

Normal_Sample, n_alt_count, n_depth, 

#annotations from Clinvar and gnomAD

ExAC2_AF, ExAC2_AF_ASJ, CLINICAL_SIGNIFICANCE, GOLD_STARS,  



If you don't have variant level and ClinVar, gnomAD annotations, please perform preliminary annotation with VariantEffectPredictor
The coordinates are expected to conform to VEP standards in order to allow for correct merging with curated training data

# Usage

Step 1 : python data_proprocessing.py --input_maf input_file --annotated_maf annotated_output_maf 
  This step adds additional annotations from myvariant.info, OncoKB and dbScSNV and also adds adds more computed annotations.
  
Step 2 : python germline_classifier.py --classifier_input output_from_previous_step --classifier_output classifier_output
  This will train the Random Forest model on provided training data and perform predictions on input file as save as output file.

Test data is provided in data folder. 

# Example
#create all relevant columns and annotations

python data_proprocessing.py --input_maf data/test_input.maf --annotated_maf data/test_annotated.maf

#run classifier

python germline_classifier.py --classifier_input data/test_annotated.maf --classifier_output data/test.classifier_output

# Output
Predictions from classifier are stored in column "prediction" as a binary value. The column "prediction_probability" shows the probability of being called pathogenic for every variant. 

Note : the scripts extract rare variants based on gnomAD frequencies and minor allele frequency within provided aggregate maf. It also extracts unique variants and subsets to variant classes that exclude UTR and silent mutations. 


