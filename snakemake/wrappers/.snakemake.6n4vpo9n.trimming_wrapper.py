
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00outputq\x03csnakemake.io\nOutputFiles\nq\x04)\x81q\x05(X,\x00\x00\x00preprocessed/SRR970578_trimmed_trimmed.fq.gzq\x06X*\x00\x00\x00preprocessed_step1/SRR970578_trimmed.fq.gzq\x07X9\x00\x00\x00preprocessed_step1/SRR970578.fastq.gz_trimming_report.txtq\x08e}q\t(X\x08\x00\x00\x00pass2_fqq\nh\x06X\x08\x00\x00\x00pass1_fqq\x0bh\x07X\x0f\x00\x00\x00pass1_fq_reportq\x0ch\x08X\x06\x00\x00\x00_namesq\r}q\x0e(h\nK\x00N\x86q\x0fh\x0bK\x01N\x86q\x10h\x0cK\x02N\x86q\x11uubX\t\x00\x00\x00wildcardsq\x12csnakemake.io\nWildcards\nq\x13)\x81q\x14X\t\x00\x00\x00SRR970578q\x15a}q\x16(h\r}q\x17X\x06\x00\x00\x00sampleq\x18K\x00N\x86q\x19sX\x06\x00\x00\x00sampleq\x1ah\x15ubX\x07\x00\x00\x00threadsq\x1bK\x01X\x06\x00\x00\x00configq\x1c}q\x1dX\x0b\x00\x00\x00config_pathq\x1eX\x19\x00\x00\x00configs/hg38_SRP029589.pyq\x1fsX\t\x00\x00\x00resourcesq csnakemake.io\nResources\nq!)\x81q"(K\x01K\x01e}q#(X\x06\x00\x00\x00_nodesq$K\x01X\x06\x00\x00\x00_coresq%K\x01h\r}q&(h$K\x00N\x86q\'h%K\x01N\x86q(uubX\x03\x00\x00\x00logq)csnakemake.io\nLog\nq*)\x81q+}q,h\r}q-sbX\x04\x00\x00\x00ruleq.X\x10\x00\x00\x00perform_trimmingq/X\x05\x00\x00\x00inputq0csnakemake.io\nInputFiles\nq1)\x81q2X\x1d\x00\x00\x00sratofastq/SRR970578.fastq.gzq3a}q4(X\x02\x00\x00\x00R1q5h3h\r}q6h5K\x00N\x86q7subX\x06\x00\x00\x00paramsq8csnakemake.io\nParams\nq9)\x81q:(K\x05X\x12\x00\x00\x00preprocessed_step1q;K\x12X\x0c\x00\x00\x00preprocessedq<X\x07\x00\x00\x00defaultq=K&e}q>(X\x0c\x00\x00\x00phred_cutoffq?K\x05X\t\x00\x00\x00pass1_dirq@h;X\n\x00\x00\x00min_lengthqAK\x12X\t\x00\x00\x00pass2_dirqBh<X\x07\x00\x00\x00adapterqCh=X\n\x00\x00\x00max_lengthqDK&h\r}qE(h?K\x00N\x86qFh@K\x01N\x86qGhAK\x02N\x86qHhBK\x03N\x86qIhCK\x04N\x86qJhDK\x05N\x86qKuubub.'); from snakemake.logging import logger; logger.printshellcmds = True
######## Original script #########
from snakemake.shell import shell
from riboraptor.kmer import fastq_kmer_histogram
import operator

# CTG is a common one
# AGATCG.. is a Truseq ribo 3' illumina
# TGGAAT.. is Truseq small rna

PREFERRED_KMERS = [
    'CTGTAGGCACCATCAAT', 'AGATCGGAAGAGCACACGTCT', 'TGGAATTCTCGGGTGCCAAGG',
    'CTGTAGGCAC'
]


def get_top_kmer(kmer_series_dict):
    """Return a kmer if it's percentage exceed 30%"""
    # Start from the longest kmer and stop
    # at where this criterion is met
    for kmer_length, kmer_series in sorted(
            kmer_series_dict.items(), key=operator.itemgetter(0),
            reverse=True):
        kmer_list = kmer_series.index.tolist()
        # Are any of the top 4 kmers from our PREFFERED_KMERS list?
        index_counter = 0
        for kmer in kmer_list:
            index_counter += 1
            if kmer in PREFERRED_KMERS:
                return kmer
            if index_counter >= 5:
                break
        # 30 is an arbitrary chosen threshold, but it works in *most* cases
        over_represented = kmer_series[kmer_series >= 30]
        if len(over_represented) >= 1:
            return over_represented.index.tolist()[0]
    return None


params = snakemake.params
pass1_dir = snakemake.params.pass1_dir
pass2_dir = snakemake.params.pass2_dir
output_1 = snakemake.output.pass1_fq
output_2 = snakemake.output.pass2_fq
adapter = snakemake.params.adapter
output_1_report = output_1 + '_trimming_report.txt'
output_2_report = output_2 + '_trimming_report.txt'

# Do first pass
if adapter is None or adapter == 'default':
    shell(r'''trim_galore -o {pass1_dir} --length {params.min_length} \
          -q {params.phred_cutoff} {snakemake.input.R1}''')
    # Are there any  over-represented sequences?
    # If yes, do a second pass
    # since no adater was provided
    histogram = fastq_kmer_histogram(output_1)
    adapter2 = get_top_kmer(histogram)
    if adapter2 is None:
        # Else just copy
        shell(r'''cp -r {output_1} {output_2}''')
        shell(
            r'''echo "No adapter found in second pass" > {output_2_report}''')
    else:
        shell(r'''trim_galore -o {pass2_dir} --length {params.min_length} \
            -a {adapter2} \
            -q {params.phred_cutoff} {output_1}''')

else:
    shell(r'''trim_galore -o {pass1_dir} --length {params.min_length} \
          -a {adapter} \
          -q {params.phred_cutoff} {snakemake.input.R1}''')
    shell(r'''cp -r {output_1} {output_2}''')
    shell(
        r'''echo "Used user provided adapter for one pass only (no second pass)" > {output_2_report}'''
    )
