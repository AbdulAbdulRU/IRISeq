import HTSeq
import os
import sys
import pandas as pd
import pickle

def gene_count_reference(gtf_file, output_folder):
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # read in the gtf file, and then construct the genome interval for exons, genes, and gene end dictionary
    gtf_file = HTSeq.GFF_Reader(gtf_file, end_included=True)
    gene_annotat_file = output_folder + "/gene_name_annotate.txt"
    exon_annotat_file = output_folder + "/exon_name_annotate.txt"
    
    gene_annotat = open(gene_annotat_file, "w")
    exon_annotat = open(exon_annotat_file, "w")
    
    exons = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    genes = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    exon_only = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
    
    gene_end = {}
    gene_start = {}
    exon_n = 0
    gene_n = 0
    transcript_n = 0
    gene_count = 0

    
    print("Start generating exon genomic arrays....")
    print("Start generating gene genomic arrays....")
    print("Start calculating transcript start and end of genes....")

    for feature in gtf_file:
        if feature.type == "exon":
            exon_n += 1
            exons[ feature.iv ] += feature.attr["gene_id"]
            exon_only[ feature.iv ] += feature.attr["exon_id"]
            
            message = (feature.attr["gene_id"] + "," + feature.attr["gene_type"] + "," 
                       + feature.attr["exon_id"] + "," + feature.attr["gene_name"]  + "\n")
            exon_annotat.write(message)
            
        elif feature.type == "gene":
            gene_n +=1
            genes[ feature.iv ] += feature.attr["gene_id"]
            gene_count += 1
            
            # for human and mouse gtf file
            message = (feature.attr["gene_id"] + "," + feature.attr["gene_type"] + "," 
                       + "gene" + "," + feature.attr["gene_name"] + "," + str(gene_count) + "\n")

            gene_annotat.write(message)
            
        elif feature.type == "transcript":
            transcript_n += 1
            #print "feature gene name: ", feature.attr["gene_id"]
            if feature.attr["gene_id"] in gene_end.keys():
                gene_end[ feature.attr["gene_id"] ].add(feature.iv.end_d)
                gene_start[ feature.attr["gene_id"] ].add(feature.iv.start_d)
            else:
                gene_end[ feature.attr["gene_id"] ] = set()
                gene_start[ feature.attr["gene_id"] ] = set()
                gene_end[ feature.attr["gene_id"] ].add(feature.iv.end_d)
                gene_start[ feature.attr["gene_id"] ].add(feature.iv.start_d)

    print("Detected gene number: ", gene_n)
    print("Detected transcript number: ", transcript_n)
    print("Detected exon number: ", exon_n)
    
    gene_annotat.close()
    exon_annotat.close()
    
    gene_annotat = pd.read_csv(gene_annotat_file, header=None)
    exon_annotate = pd.read_csv(exon_annotat_file, header=None)

    #print("print WBGENE id:", gene_annotat.loc["WBGene00004947", 4])
    #print("Print transcript end, ", len(gene_end))
    print("Save all files...")
    gene_reference = {}
    gene_reference["genes"] = genes
    gene_reference["exons"] = exons
    gene_reference["exon_only"] = exon_only
    gene_reference["gene_end"] = gene_end
    gene_reference["gene_start"] = gene_start
    gene_reference["gene_annotat"] = gene_annotat
    gene_reference["exon_annotat"] = exon_annotate
    
    output_file = output_folder + "/" + "Gene_reference.pickle"
    with open(output_file, "wb") as f:
        pickle.dump(gene_reference, f)

    print("All analysis done~")
'''    
if __name__ == "__main__":
    gtf_file = sys.argv[1]
    output_folder = sys.argv[2]
    gene_count_reference(gtf_file, output_folder)
'''