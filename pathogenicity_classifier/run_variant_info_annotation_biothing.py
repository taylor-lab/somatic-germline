import sys
import myvariant
import os
import pandas as pd
batch_size = 1000
import json
import gzip
from biothings_client import get_client

line_count = 0
batch_count = 0
query_batch = []
mv = get_client('variant')
input_file = sys.argv[1]
output_dir = sys.argv[2]
print ("running myvariant info .."+input_file)
with open(input_file) as fin:
    Line = fin.readline()
    Line = fin.readline()
    while Line:
        elements = Line.strip().split("\t")
        end_position = elements[4]
        variant_type = elements[5]
        if variant_type=='SNP':
            strquery = "chr"+elements[0]+":g."+elements[1]+elements[2]+">"+elements[3]
        elif variant_type == 'INS':
            #chr2:g.17142_17143insA
            strquery = "chr"+elements[0]+":g."+elements[1]+"_"+end_position+"ins"+elements[3]
        elif variant_type == 'DEL':
            #chrMT:g.8271_8279del
            strquery = "chr"+elements[0]+":g."+elements[1]+"_"+end_position+"del"

        if line_count % batch_size == 0 :
            if(len(query_batch)>0):

                batch_count = batch_count + 1
                output = mv.getvariants(query_batch,
                    fields="dbnsfp,clinvar,evs,cadd,gwassnps,cosmic,docm,snpedia,emv,grasp,civic,cgi", as_dataframe=True)
                output_file = output_dir+ "/batch"+"_"+str(batch_count)+".txt"
		output.to_csv(output_file, sep="\t", encoding='utf-8')
            query_batch = []
            query_batch.append(strquery)

        else:
            query_batch.append(strquery)
        line_count = line_count + 1
        Line = fin.readline()
	
batch_count = batch_count + 1
output = mv.getvariants(query_batch, 
		fields="dbnsfp,clinvar,evs,cadd,gwassnps,cosmic,docm,snpedia,emv,grasp,civic,cgi", as_dataframe=True)
output_file = output_dir+ "/batch"+"_"+str(batch_count)+".txt"
output.to_csv(output_file, sep="\t", encoding='utf-8')
