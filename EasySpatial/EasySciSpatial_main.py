#! /ru-auth/local/home/jcao/anaconda3/bin/python
#SBATCH --partition=hpc,cao,cao_bigmem
#SBATCH --nodes=1 
#SBATCH --ntasks=2 #Annotate the number of cores
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --job-name="JC_20230716_EasySci_test" # annotate job name
#SBATCH -o /ru-auth/local/home/jcao/Report_slurm/output.JC_20230716_EasySci_test
#SBATCH --mail-user=jcao@rockefeller.edu
#SBATCH --mail-type=ALL

import sys
import subprocess
import os
import gzip
from multiprocessing import Pool

# Define Easysci script folder
EasySpatial_script_folder = '/ru-auth/local/home/jcao/Script/EasySpatial/'
sys.path.append(EasySpatial_script_folder)
from Spatial_UMI_barcode_extraction import extract_spatial_barcode_files
from Fastq_trim_multi_files import Fastq_trim_files
from STAR import Fastq_star_alignment_multi_files
from Sam_filter_multi_files import Sam_filter_files
from Sam_rm_dup_barcode_UMI_multi_files import rm_dup_files
from Sam_gene_counting_multi_files import scRNA_count_parallel
from Summary_gene_count_multi_files import Gene_count_summary
from Generate_adata import generate_adata_from_gene_count
from File_functions import *
from Count_reads import *

# Define the input folder, output folder and core number
fastq_folder = "/ru-auth/local/home/jcao/projects/JC_20200628_sciRNAseqDev/Raw_data/JC_20230714_test_spatial_pipeline/Fastq/sub_sampled/"
sample_ID = "/ru-auth/local/home/jcao/projects/JC_20200628_sciRNAseqDev/Raw_data/JC_20230714_test_spatial_pipeline/Fastq/sample_ID.txt"
output_folder = "/ru-auth/local/home/jcao/projects/JC_20200628_sciRNAseqDev/Raw_data/JC_20230714_test_spatial_pipeline/output"
core = 2

print("\n******* Fastq_folder:", fastq_folder)
print("\n******* Sample ID:", sample_ID)
print("\n******* Output folder:", output_folder)
print("\n******* Core number:", core)

# Define the location of common tools for sequencing processing
star_path="/rugpfs/fs0/cao_lab/scratch/asziraki/anaconda3/envs/original_pipeline/bin/STAR"
cutadapt_path = "/ru-auth/local/home/jcao/anaconda3_new/envs/cutadaptenv/bin/cutadapt"
samtools_path = "/rugpfs/fs0/cao_lab/scratch/asziraki/anaconda3/envs/original_pipeline/bin/samtools"

# define the index and gtf files for alignment and gene counting
index="/ru-auth/local/home/jcao/store/Reference/Index/STAR_hg19_mm10_RNAseq"
gtf_file="/ru-auth/local/home/jcao/store/Reference/GTF/rmchr.gencode.v19.chr_patch_hapl_scaff.annotation.gencode.vM12.chr_patch_hapl_scaff.annotation.gtf.gz"
gtf_annotation_file = "/ru-auth/local/home/jcao/store/Reference/Gene_annotation/hg19_mm10_Gene_reference.pickle"

# define the location of the beads barcode files
list_barcode_1_file="/ru-auth/local/home/aabdul/restored_projects_CaoStore/restored/projects/Making_Dictionary_1/AA_20230618_spatial_data_barcode_reference/Spatial_R2_barcode_1.pickle"
list_barcode_2_file="/ru-auth/local/home/aabdul/restored_projects_CaoStore/restored/projects/Making_Dictionary_1/AA_20230618_spatial_data_barcode_reference/Spatial_R2_barcode_2.pickle"
list_barcode_3_file="/ru-auth/local/home/aabdul/restored_projects_CaoStore/restored/projects/Making_Dictionary_1/AA_20230618_spatial_data_barcode_reference/Spatial_R2_barcode_3.pickle"
list_barcode_4_file="/ru-auth/local/home/aabdul/restored_projects_CaoStore/restored/projects/Making_Dictionary_1/AA_20230618_spatial_data_barcode_reference/Spatial_R2_barcode_4_bead1.pickle"
mismatch_distance=1

# Define intermediate output folders
Fastq_barcode_attached = output_folder + "Fastq_barcode_attached/"
dir_make(Fastq_barcode_attached)

Fastq_trimmed = output_folder + "Fastq_trimmed/"
dir_make(Fastq_trimmed)

Sam_STAR = output_folder + "Sam_STAR/"
dir_make(Sam_STAR)

Sam_filtered = output_folder + "Sam_filtered"
dir_make(Sam_filtered)

Sam_rmdup = output_folder + "Sam_rmdup"
dir_make(Sam_rmdup)

Bed_gene_count = output_folder + "Bed_gene_count"
dir_make(Bed_gene_count)

Summary_gene_count = output_folder + "Summary_gene_count"
dir_make(Summary_gene_count)

Adata_folder = output_folder + "Adata"
dir_make(Adata_folder)

################# Update the name of the fastq files
input_command = f"for sample in $(cat {sample_ID}); do echo changing name $sample; mv {fastq_folder}/*$sample*R1*.fastq.gz {fastq_folder}/$sample.R1.fastq.gz; mv {fastq_folder}/*$sample*R2*.fastq.gz {fastq_folder}/$sample.R2.fastq.gz; mv {fastq_folder}/*$sample*R3*.fastq.gz {fastq_folder}/$sample.R3.fastq.gz; done"
result = subprocess.run(input_command, shell=True, text=True, )
print(result)

################# Extract barocde for attahment to UMI
extract_spatial_barcode_files(fastq_folder, sample_ID, Fastq_barcode_attached, core, list_barcode_1_file, list_barcode_2_file, list_barcode_3_file, list_barcode_4_file)

################# Trimming the read2
Fastq_trim_files(Fastq_barcode_attached, sample_ID, Fastq_trimmed, core)

################# Align the read2
Fastq_star_alignment_multi_files(Fastq_trimmed, sample_ID, Sam_STAR, core, index, star_path)

################# Filter the reads
Sam_filter_files(Sam_STAR, sample_ID, Sam_filtered, core)

################# Remove duplicated reads
rm_dup_files(Sam_filtered, sample_ID, Sam_rmdup, core)

################# Generate gene count files
scRNA_count_parallel(Sam_rmdup, sample_ID, Bed_gene_count, gtf_annotation_file, core)

################# Generate summarized gene count files
Gene_count_summary(Bed_gene_count, sample_ID, Summary_gene_count)

################# Generate adata object from the gene count data
generate_adata_from_gene_count(Summary_gene_count, Bed_gene_count + "/gene_anno.csv", Adata_folder, 50)

################# Now we are going to calculate the number of reads following each step of processing
df_count = pd.DataFrame({"Sample_name" : read_csv_to_list(sample_ID)})
df_count["Count_Fastq_raw"] = Fastq_count_reads_files(fastq_folder, sample_ID)
df_tmp = Count_Align_STAR_files(Sam_STAR, sample_ID)
df_count["Count_Fastq_filtered"] = Fastq_count_reads_files(Fastq_barcode_attached, sample_ID)
df_count["Count_Fastq_trimmed"] = list(df_tmp[0])
df_count["Count_Sam_mapped"] = list(df_tmp[1] + df_tmp[2])
df_count["Count_Sam_unique"] = list(df_tmp[1])
#df_count["Count_Fastq_filtered"] = Fastq_count_reads_files(Fastq_barcode_attached, sample_ID)
#df_count["Count_Sam_mapped"] = SAM_count_mapped_reads_files(Sam_STAR, sample_ID)
#df_count["Count_Sam_unique"] = SAM_count_reads_files(Sam_filtered, sample_ID)
df_count["Count_Sam_rmdup"] = SAM_count_reads_files(Sam_rmdup, sample_ID)
df_count["Count_Sam_annotated"] = count_mapped_reads_files(Bed_gene_count, sample_ID)
df_count["Count_Cell_reads"] = Count_cell_reads(Adata_folder)
df_count.to_csv(Adata_folder + "/df_reads.csv")

################# Remove all intermediate folders
print("Remove all intermediate folders....")
dir_rm(Fastq_barcode_attached)
dir_rm(Fastq_trimmed)
dir_rm(Sam_STAR)
dir_rm(Sam_filtered)
dir_rm(Sam_rmdup)
dir_rm(Bed_gene_count)
dir_rm(Summary_gene_count)

print("All analysis done")
