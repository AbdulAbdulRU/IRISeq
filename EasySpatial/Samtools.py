# define a function for filtering reads with samtools
def Sam_filter_sort(input_file, output_file,  
                    samtools_path="/rugpfs/fs0/cao_lab/scratch/asziraki/anaconda3/envs/original_pipeline/bin/samtools"):
    input_command = f"{samtools_path} view -bh -q 30 -F 4 {input_file}|{samtools_path} sort -@ 10 -|{samtools_path} view -h ->{output_file}"
    print(input_command)
    result = subprocess.check_output(input_command, shell=True, text=True)
    print(result)
