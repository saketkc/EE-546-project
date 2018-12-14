
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\t\x00\x00\x00resourcesq\x03csnakemake.io\nResources\nq\x04)\x81q\x05(K\x01K\x01e}q\x06(X\x06\x00\x00\x00_coresq\x07K\x01X\x06\x00\x00\x00_namesq\x08}q\t(h\x07K\x00N\x86q\nX\x06\x00\x00\x00_nodesq\x0bK\x01N\x86q\x0cuh\x0bK\x01ubX\x04\x00\x00\x00ruleq\rX\x10\x00\x00\x00perform_trimmingq\x0eX\x07\x00\x00\x00threadsq\x0fK\x01X\x06\x00\x00\x00outputq\x10csnakemake.io\nOutputFiles\nq\x11)\x81q\x12(X9\x00\x00\x00preprocessed_step1/SRR970588.fastq.gz_trimming_report.txtq\x13X*\x00\x00\x00preprocessed_step1/SRR970588_trimmed.fq.gzq\x14X,\x00\x00\x00preprocessed/SRR970588_trimmed_trimmed.fq.gzq\x15e}q\x16(h\x08}q\x17(X\x0f\x00\x00\x00pass1_fq_reportq\x18K\x00N\x86q\x19X\x08\x00\x00\x00pass1_fqq\x1aK\x01N\x86q\x1bX\x08\x00\x00\x00pass2_fqq\x1cK\x02N\x86q\x1duh\x1ah\x14h\x1ch\x15h\x18h\x13ubX\t\x00\x00\x00wildcardsq\x1ecsnakemake.io\nWildcards\nq\x1f)\x81q X\t\x00\x00\x00SRR970588q!a}q"(h\x08}q#X\x06\x00\x00\x00sampleq$K\x00N\x86q%sX\x06\x00\x00\x00sampleq&h!ubX\x06\x00\x00\x00configq\'}q(X\x0b\x00\x00\x00config_pathq)X\x19\x00\x00\x00configs/hg38_SRP029589.pyq*sX\x05\x00\x00\x00inputq+csnakemake.io\nInputFiles\nq,)\x81q-X\x1d\x00\x00\x00sratofastq/SRR970588.fastq.gzq.a}q/(h\x08}q0X\x02\x00\x00\x00R1q1K\x00N\x86q2sh1h.ubX\x03\x00\x00\x00logq3csnakemake.io\nLog\nq4)\x81q5}q6h\x08}q7sbX\x06\x00\x00\x00paramsq8csnakemake.io\nParams\nq9)\x81q:(X\x12\x00\x00\x00preprocessed_step1q;X\x0c\x00\x00\x00preprocessedq<X\x07\x00\x00\x00defaultq=K\x05K\x12K&e}q>(X\t\x00\x00\x00pass1_dirq?h;X\t\x00\x00\x00pass2_dirq@h<X\x07\x00\x00\x00adapterqAh=X\x0c\x00\x00\x00phred_cutoffqBK\x05X\n\x00\x00\x00min_lengthqCK\x12h\x08}qD(h?K\x00N\x86qEh@K\x01N\x86qFhAK\x02N\x86qGhBK\x03N\x86qHhCK\x04N\x86qIX\n\x00\x00\x00max_lengthqJK\x05N\x86qKuhJK&ubub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
