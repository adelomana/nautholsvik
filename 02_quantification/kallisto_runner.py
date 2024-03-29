###
### usage: time python kallisto_runner.py &> messages.kallisto.txt
###

import sys, datetime, os

def kallisto_caller(label):

    printt('about to quantify {}'.format(label))

    sample_output_dir = results_dir + label
    executable = 'time kallisto quant'
    options = ' -i {} -o {} --bias -t {} -b {} {} '.format(transcriptome_index, sample_output_dir, threads, boots, strand_flag)

    # define fastq files
    found_files = os.listdir(clean_fastq_dir + label)
    all_fastq_roots = []
    for element in found_files:
        root = element.split('_R')[0]
        all_fastq_roots.append(root)
    fastq_pair_labels = list(set(all_fastq_roots))

    fastq_files = ''
    for element in fastq_pair_labels:
        pair_string = '{} {} '.format(clean_fastq_dir + label + '/' + element + '_R1_clean.fastq.gz', clean_fastq_dir + label + '/' + element + '_R2_clean.fastq.gz')
        fastq_files = fastq_files + pair_string
    # end define fastq files

    command = executable + options + fastq_files

    print('')
    print(command)
    os.system(command)
    print('')

    return None

def printt(label):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(label)))

    return None

#
# 0. user defined variables
#
boots = 100
threads = 20
clean_fastq_dir = '/home/adrian/projects/nautholsvik/data/clean_fastq/'
results_dir = '/home/adrian/projects/nautholsvik/results/kallisto/kallisto.{}/'.format(boots)
transcriptome_index = '/home/adrian/software/kallisto/ensembl_v96/hsa/transcriptome.idx'

strand_flag = '--rf-stranded'       # [quant] processed 130,967,601 reads, 58,495,796 reads pseudoaligned
strand_flag = '--fr-stranded'       # [quant] processed 130,967,601 reads, 57,583,752 reads pseudoaligned
strand_flag = ''                    # [quant] processed 130,967,601 reads, 116,270,207 reads pseudoaligned

#
# 1. recover labels
#
printt('recover labels...')

labels = next(os.walk(clean_fastq_dir))[1]
labels.sort()
print(labels, len(labels))

#
# 2. call kallisto quant
#
if os.path.exists(results_dir) == False:
    os.mkdir(results_dir)

for label in labels:
    kallisto_caller(label)
