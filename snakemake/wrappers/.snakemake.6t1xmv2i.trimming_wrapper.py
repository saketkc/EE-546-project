
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00outputq\x03csnakemake.io\nOutputFiles\nq\x04)\x81q\x05(X9\x00\x00\x00preprocessed_step1/SRR970593.fastq.gz_trimming_report.txtq\x06X*\x00\x00\x00preprocessed_step1/SRR970593_trimmed.fq.gzq\x07X,\x00\x00\x00preprocessed/SRR970593_trimmed_trimmed.fq.gzq\x08e}q\t(X\x06\x00\x00\x00_namesq\n}q\x0b(X\x0f\x00\x00\x00pass1_fq_reportq\x0cK\x00N\x86q\rX\x08\x00\x00\x00pass1_fqq\x0eK\x01N\x86q\x0fX\x08\x00\x00\x00pass2_fqq\x10K\x02N\x86q\x11uh\x0ch\x06h\x0eh\x07h\x10h\x08ubX\t\x00\x00\x00resourcesq\x12csnakemake.io\nResources\nq\x13)\x81q\x14(K\x01K\x01e}q\x15(X\x06\x00\x00\x00_nodesq\x16K\x01h\n}q\x17(X\x06\x00\x00\x00_coresq\x18K\x00N\x86q\x19h\x16K\x01N\x86q\x1auh\x18K\x01ubX\x07\x00\x00\x00threadsq\x1bK\x01X\x03\x00\x00\x00logq\x1ccsnakemake.io\nLog\nq\x1d)\x81q\x1e}q\x1fh\n}q sbX\x05\x00\x00\x00inputq!csnakemake.io\nInputFiles\nq")\x81q#X\x1d\x00\x00\x00sratofastq/SRR970593.fastq.gzq$a}q%(h\n}q&X\x02\x00\x00\x00R1q\'K\x00N\x86q(sh\'h$ubX\x06\x00\x00\x00paramsq)csnakemake.io\nParams\nq*)\x81q+(K\x12X\x12\x00\x00\x00preprocessed_step1q,X\x07\x00\x00\x00defaultq-K\x05K&X\x0c\x00\x00\x00preprocessedq.e}q/(X\n\x00\x00\x00min_lengthq0K\x12h\n}q1(h0K\x00N\x86q2X\t\x00\x00\x00pass1_dirq3K\x01N\x86q4X\x07\x00\x00\x00adapterq5K\x02N\x86q6X\x0c\x00\x00\x00phred_cutoffq7K\x03N\x86q8X\n\x00\x00\x00max_lengthq9K\x04N\x86q:X\t\x00\x00\x00pass2_dirq;K\x05N\x86q<uh5h-h3h,h7K\x05h;h.h9K&ubX\t\x00\x00\x00wildcardsq=csnakemake.io\nWildcards\nq>)\x81q?X\t\x00\x00\x00SRR970593q@a}qA(h\n}qBX\x06\x00\x00\x00sampleqCK\x00N\x86qDsX\x06\x00\x00\x00sampleqEh@ubX\x04\x00\x00\x00ruleqFX\x10\x00\x00\x00perform_trimmingqGX\x06\x00\x00\x00configqH}qIX\x0b\x00\x00\x00config_pathqJX\x19\x00\x00\x00configs/hg38_SRP029589.pyqKsub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
