######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\t\x00\x00\x00wildcardsq\x03csnakemake.io\nWildcards\nq\x04)\x81q\x05X\n\x00\x00\x00SRR5680922q\x06a}q\x07(X\x06\x00\x00\x00_namesq\x08}q\tX\x06\x00\x00\x00sampleq\nK\x00N\x86q\x0bsX\x06\x00\x00\x00sampleq\x0ch\x06ubX\x04\x00\x00\x00ruleq\rX\x10\x00\x00\x00perform_trimmingq\x0eX\t\x00\x00\x00resourcesq\x0fcsnakemake.io\nResources\nq\x10)\x81q\x11(K\x01K\x01e}q\x12(h\x08}q\x13(X\x06\x00\x00\x00_coresq\x14K\x00N\x86q\x15X\x06\x00\x00\x00_nodesq\x16K\x01N\x86q\x17uh\x16K\x01h\x14K\x01ubX\x03\x00\x00\x00logq\x18csnakemake.io\nLog\nq\x19)\x81q\x1a}q\x1bh\x08}q\x1csbX\x05\x00\x00\x00inputq\x1dcsnakemake.io\nInputFiles\nq\x1e)\x81q\x1fX\x1e\x00\x00\x00sratofastq/SRR5680922.fastq.gzq a}q!(h\x08}q"X\x02\x00\x00\x00R1q#K\x00N\x86q$sh#h ubX\x06\x00\x00\x00outputq%csnakemake.io\nOutputFiles\nq&)\x81q\'(X-\x00\x00\x00preprocessed/SRR5680922_trimmed_trimmed.fq.gzq(X+\x00\x00\x00preprocessed_step1/SRR5680922_trimmed.fq.gzq)X:\x00\x00\x00preprocessed_step1/SRR5680922.fastq.gz_trimming_report.txtq*e}q+(h\x08}q,(X\x08\x00\x00\x00pass2_fqq-K\x00N\x86q.X\x08\x00\x00\x00pass1_fqq/K\x01N\x86q0X\x0f\x00\x00\x00pass1_fq_reportq1K\x02N\x86q2uh-h(h/h)h1h*ubX\x07\x00\x00\x00threadsq3K\x01X\x06\x00\x00\x00paramsq4csnakemake.io\nParams\nq5)\x81q6(X\x0c\x00\x00\x00preprocessedq7X\x12\x00\x00\x00preprocessed_step1q8K\x05X\x07\x00\x00\x00defaultq9K\x12K&e}q:(X\t\x00\x00\x00pass1_dirq;h8X\t\x00\x00\x00pass2_dirq<h7h\x08}q=(h<K\x00N\x86q>h;K\x01N\x86q?X\x0c\x00\x00\x00phred_cutoffq@K\x02N\x86qAX\x07\x00\x00\x00adapterqBK\x03N\x86qCX\n\x00\x00\x00min_lengthqDK\x04N\x86qEX\n\x00\x00\x00max_lengthqFK\x05N\x86qGuh@K\x05hBh9hDK\x12hFK&ubX\x06\x00\x00\x00configqH}qIX\x0b\x00\x00\x00config_pathqJX\x19\x00\x00\x00configs/hg38_SRP109126.pyqKsub.'
)
from snakemake.logging import logger
logger.printshellcmds = True
######## Original script #########
from snakemake.shell import shell
from riboraptor.kmer import fastq_kmer_histogram
import operator


def get_top_kmer(kmer_series_dict):
    """Return a kmer if it's percentage exceed 30%"""
    # Start from the longest kmer and stop
    # at where this criterion is met
    for kmer_length, kmer_series in sorted(
            kmer_series_dict.items(), key=operator.itemgetter(0),
            reverse=True):
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
