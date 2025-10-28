import sys
import gzip
import pickle
from multiprocessing import Pool
from functools import partial

def UMI_attach_read2_barcode_list(sample, input_folder, output_folder, ligation_barcode_list, RT_barcode_list, mismatch_rate = 1):
    Read1 = input_folder + "/" + sample + ".R1.fastq.gz"
    Read2 = input_folder + "/" + sample + ".R3.fastq.gz"
    Read3 = input_folder + "/" + sample + ".R2.fastq.gz"
    output_file = output_folder + "/" + sample + ".R2.fastq.gz"
    mismatch_rate = int(mismatch_rate)
    f1 = gzip.open(Read1, 'rt')
    f2 = gzip.open(Read2, 'rt')
    f3 = gzip.open(output_file, 'wt')
    f4 = gzip.open(Read3, 'rt')
    
    line1 = f1.readline()
    line2 = f2.readline()
    line3 = f4.readline()
    total_line = 0
    filtered_line = 0
    
    while (line1):
        total_line += 1
        line1 = f1.readline()
        line3 = f4.readline()

        tmp_lig = line3[0:10]

        if tmp_lig in ligation_barcode_list:
            ligation_bc_match = ligation_barcode_list[tmp_lig]
            target_RT = line1[8:18]

            if target_RT in RT_barcode_list:
                barcode = RT_barcode_list[target_RT]
                filtered_line += 1
                UMI = line1[:8]
                first_line = '@' + ligation_bc_match + barcode + ',' + UMI + ',' + line2[1:]
                f3.write(first_line)

                second_line = f2.readline()
                f3.write(second_line)

                third_line = f2.readline()
                f3.write(third_line)

                four_line = f2.readline()
                f3.write(four_line)

                line2 = f2.readline()
            
            else:
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                
        else:
            line2 = f2.readline()
            line2 = f2.readline()
            line2 = f2.readline()
            line2 = f2.readline()

        line1 = f1.readline()
        line1 = f1.readline()
        line1 = f1.readline()

        line3 = f4.readline() 
        line3 = f4.readline()
        line3 = f4.readline()

    f1.close()
    f2.close()
    f3.close()
    f4.close()
    print("sample name: %s, total line: %f, filtered line: %f, filter rate: %f" 
          %(sample, total_line, filtered_line, float(filtered_line) / float(total_line)))

def attach_UMI_files(input_folder, sampleID, output_folder, ligation_barcode_file, RT_barcode_file, core):
    init_message = '''
    --------------------------start attaching UMI-----------------------------
    input folder: %s
    sample ID: %s
    output_folder: %s
    ligation barcode file: %s
    RT barcode file: %s
    ___________________________________________________________________________
    ''' %(input_folder, sampleID, output_folder, ligation_barcode_file, RT_barcode_file)
    
    print(init_message)
    
    print("Load ligation barcode dictionary...")
    
    barcodes = open(ligation_barcode_file, "rb")
    ligation_barcode_list = pickle.load(barcodes)
    barcodes.close()
    
    print("Load RT barcode dictionary...")
    
    barcodes = open(RT_barcode_file, "rb")
    RT_barcode_list = pickle.load(barcodes)
    barcodes.close()
    
    sample_file = open(sampleID)
    sample_list = []
    for line in sample_file:
        sample = line.strip()
        sample_list.append(sample)
    sample_file.close()
    #print(sample_list)
    
    p = Pool(processes = int(core))
    func = partial(UMI_attach_read2_barcode_list, input_folder = input_folder, output_folder=output_folder, ligation_barcode_list = ligation_barcode_list, RT_barcode_list=RT_barcode_list, mismatch_rate = 1)
    result = p.map(func, sample_list)
    p.close()
    p.join()
    
    com_message = '''~~~~~~~~~~~~~~~UMI attachment done~~~~~~~~~~~~~~~~~~'''
    print(com_message)
'''
if __name__ == "__main__":
    input_folder = sys.argv[1]
    sampleID = sys.argv[2]
    output_folder = sys.argv[3]
    ligation_barcode_file = sys.argv[4]
    RT_barcode_file = sys.argv[5]
    core=sys.argv[6]
    attach_UMI_files(input_folder, sampleID, output_folder, ligation_barcode_file, RT_barcode_file, core)

'''
