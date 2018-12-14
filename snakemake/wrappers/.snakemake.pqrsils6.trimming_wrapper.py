
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00paramsq\x03csnakemake.io\nParams\nq\x04)\x81q\x05(K\x05X\x0c\x00\x00\x00preprocessedq\x06X\x07\x00\x00\x00defaultq\x07K\x12K&X\x12\x00\x00\x00preprocessed_step1q\x08e}q\t(X\x0c\x00\x00\x00phred_cutoffq\nK\x05X\t\x00\x00\x00pass2_dirq\x0bh\x06X\x06\x00\x00\x00_namesq\x0c}q\r(h\nK\x00N\x86q\x0eh\x0bK\x01N\x86q\x0fX\x07\x00\x00\x00adapterq\x10K\x02N\x86q\x11X\n\x00\x00\x00min_lengthq\x12K\x03N\x86q\x13X\n\x00\x00\x00max_lengthq\x14K\x04N\x86q\x15X\t\x00\x00\x00pass1_dirq\x16K\x05N\x86q\x17uh\x10h\x07h\x12K\x12h\x14K&h\x16h\x08ubX\x05\x00\x00\x00inputq\x18csnakemake.io\nInputFiles\nq\x19)\x81q\x1aX\x1d\x00\x00\x00sratofastq/SRR970565.fastq.gzq\x1ba}q\x1c(h\x0c}q\x1dX\x02\x00\x00\x00R1q\x1eK\x00N\x86q\x1fsh\x1eh\x1bubX\t\x00\x00\x00resourcesq csnakemake.io\nResources\nq!)\x81q"(K\x01K\x01e}q#(X\x06\x00\x00\x00_nodesq$K\x01h\x0c}q%(X\x06\x00\x00\x00_coresq&K\x00N\x86q\'h$K\x01N\x86q(uh&K\x01ubX\x06\x00\x00\x00outputq)csnakemake.io\nOutputFiles\nq*)\x81q+(X,\x00\x00\x00preprocessed/SRR970565_trimmed_trimmed.fq.gzq,X*\x00\x00\x00preprocessed_step1/SRR970565_trimmed.fq.gzq-X9\x00\x00\x00preprocessed_step1/SRR970565.fastq.gz_trimming_report.txtq.e}q/(X\x08\x00\x00\x00pass2_fqq0h,X\x08\x00\x00\x00pass1_fqq1h-h\x0c}q2(h0K\x00N\x86q3h1K\x01N\x86q4X\x0f\x00\x00\x00pass1_fq_reportq5K\x02N\x86q6uh5h.ubX\t\x00\x00\x00wildcardsq7csnakemake.io\nWildcards\nq8)\x81q9X\t\x00\x00\x00SRR970565q:a}q;(h\x0c}q<X\x06\x00\x00\x00sampleq=K\x00N\x86q>sX\x06\x00\x00\x00sampleq?h:ubX\x04\x00\x00\x00ruleq@X\x10\x00\x00\x00perform_trimmingqAX\x03\x00\x00\x00logqBcsnakemake.io\nLog\nqC)\x81qD}qEh\x0c}qFsbX\x07\x00\x00\x00threadsqGK\x01X\x06\x00\x00\x00configqH}qIX\x0b\x00\x00\x00config_pathqJX\x19\x00\x00\x00configs/hg38_SRP029589.pyqKsub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
