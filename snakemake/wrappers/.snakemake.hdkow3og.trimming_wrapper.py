
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x04\x00\x00\x00ruleq\x03X\x10\x00\x00\x00perform_trimmingq\x04X\x05\x00\x00\x00inputq\x05csnakemake.io\nInputFiles\nq\x06)\x81q\x07X\x1d\x00\x00\x00sratofastq/SRR970578.fastq.gzq\x08a}q\t(X\x06\x00\x00\x00_namesq\n}q\x0bX\x02\x00\x00\x00R1q\x0cK\x00N\x86q\rsh\x0ch\x08ubX\t\x00\x00\x00wildcardsq\x0ecsnakemake.io\nWildcards\nq\x0f)\x81q\x10X\t\x00\x00\x00SRR970578q\x11a}q\x12(X\x06\x00\x00\x00sampleq\x13h\x11h\n}q\x14X\x06\x00\x00\x00sampleq\x15K\x00N\x86q\x16subX\x06\x00\x00\x00outputq\x17csnakemake.io\nOutputFiles\nq\x18)\x81q\x19(X*\x00\x00\x00preprocessed_step1/SRR970578_trimmed.fq.gzq\x1aX,\x00\x00\x00preprocessed/SRR970578_trimmed_trimmed.fq.gzq\x1bX9\x00\x00\x00preprocessed_step1/SRR970578.fastq.gz_trimming_report.txtq\x1ce}q\x1d(X\x08\x00\x00\x00pass1_fqq\x1eh\x1aX\x08\x00\x00\x00pass2_fqq\x1fh\x1bh\n}q (h\x1eK\x00N\x86q!h\x1fK\x01N\x86q"X\x0f\x00\x00\x00pass1_fq_reportq#K\x02N\x86q$uh#h\x1cubX\t\x00\x00\x00resourcesq%csnakemake.io\nResources\nq&)\x81q\'(K\x01K\x01e}q((X\x06\x00\x00\x00_nodesq)K\x01h\n}q*(X\x06\x00\x00\x00_coresq+K\x00N\x86q,h)K\x01N\x86q-uh+K\x01ubX\x07\x00\x00\x00threadsq.K\x01X\x03\x00\x00\x00logq/csnakemake.io\nLog\nq0)\x81q1}q2h\n}q3sbX\x06\x00\x00\x00paramsq4csnakemake.io\nParams\nq5)\x81q6(K\x05X\x0c\x00\x00\x00preprocessedq7K&K\x12X\x07\x00\x00\x00defaultq8X\x12\x00\x00\x00preprocessed_step1q9e}q:(X\x0c\x00\x00\x00phred_cutoffq;K\x05X\t\x00\x00\x00pass2_dirq<h7X\n\x00\x00\x00max_lengthq=K&X\n\x00\x00\x00min_lengthq>K\x12h\n}q?(h;K\x00N\x86q@h<K\x01N\x86qAh=K\x02N\x86qBh>K\x03N\x86qCX\x07\x00\x00\x00adapterqDK\x04N\x86qEX\t\x00\x00\x00pass1_dirqFK\x05N\x86qGuhDh8hFh9ubX\x06\x00\x00\x00configqH}qIX\x0b\x00\x00\x00config_pathqJX\x19\x00\x00\x00configs/hg38_SRP029589.pyqKsub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
