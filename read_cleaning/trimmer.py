###
### usage: time python trimmer.py &> messages.txt
###

import os, datetime, sys

def printt(label):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(label)))

    return None

def trimmomatic_caller(sample):

    executable='time java -jar {}trimmomatic-0.39.jar PE -threads {} -phred33 '.format(trimmomatic_path,number_threads)
    options=' ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'.format(adapter_file)

    inputs = os.listdir(raw_fastq_dir + sample)
    inputs.sort()

    input1 = raw_fastq_dir + sample + '/' + inputs[0]
    input2 = raw_fastq_dir + sample + '/' + inputs[1]

    output_dir = clean_fastq_dir + sample + '/'
    os.mkdir(output_dir)

    output1 = output_dir + sample + '_R1_clean.fastq.gz'
    output2 = output_dir + sample + '_R2_clean.fastq.gz'

    garbage1 = output_dir + sample + '_R1_garbage.fastq.gz'
    garbage2 = output_dir + sample + '_R2_garbage.fastq.gz'

    input_files = input1 + ' ' + input2
    output_files = output1 + ' ' + garbage1 + ' ' + output2 + ' ' + garbage2

    command = executable + input_files + ' ' + output_files + options

    printt('about to clean {}'.format(sample))
    print('')
    print(command)
    print('')
    os.system(command)
    print('')

    return None

# 0. user defined variablessample
raw_fastq_dir = '/home/adrian/projects/nautholsvik/data/decode_220406/'
clean_fastq_dir = '/home/adrian/projects/nautholsvik/data/clean_fastq/'
trimmomatic_path = '/home/adrian/software/Trimmomatic-0.39/'
adapter_file = trimmomatic_path + 'adapters/TruSeq3-PE-2.fa'
number_threads = 20

# 1. recover samples
working_dirs = os.listdir(raw_fastq_dir)
working_dirs.sort()
print(working_dirs, len(working_dirs))

# 2. iterate Trimmomatic
for working_dir in working_dirs:
    trimmomatic_caller(working_dir)
