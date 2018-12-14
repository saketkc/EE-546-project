
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00paramsq\x03csnakemake.io\nParams\nq\x04)\x81q\x05(K&X\x07\x00\x00\x00defaultq\x06X\x12\x00\x00\x00preprocessed_step1q\x07K\x12X\x0c\x00\x00\x00preprocessedq\x08K\x05e}q\t(X\n\x00\x00\x00max_lengthq\nK&X\n\x00\x00\x00min_lengthq\x0bK\x12X\x06\x00\x00\x00_namesq\x0c}q\r(h\nK\x00N\x86q\x0eh\x0bK\x03N\x86q\x0fX\t\x00\x00\x00pass1_dirq\x10K\x02N\x86q\x11X\x07\x00\x00\x00adapterq\x12K\x01N\x86q\x13X\t\x00\x00\x00pass2_dirq\x14K\x04N\x86q\x15X\x0c\x00\x00\x00phred_cutoffq\x16K\x05N\x86q\x17uh\x10h\x07h\x12h\x06h\x14h\x08h\x16K\x05ubX\x06\x00\x00\x00outputq\x18csnakemake.io\nOutputFiles\nq\x19)\x81q\x1a(X*\x00\x00\x00preprocessed_step1/SRR970588_trimmed.fq.gzq\x1bX9\x00\x00\x00preprocessed_step1/SRR970588.fastq.gz_trimming_report.txtq\x1cX,\x00\x00\x00preprocessed/SRR970588_trimmed_trimmed.fq.gzq\x1de}q\x1e(X\x08\x00\x00\x00pass1_fqq\x1fh\x1bX\x0f\x00\x00\x00pass1_fq_reportq h\x1ch\x0c}q!(h\x1fK\x00N\x86q"h K\x01N\x86q#X\x08\x00\x00\x00pass2_fqq$K\x02N\x86q%uh$h\x1dubX\x03\x00\x00\x00logq&csnakemake.io\nLog\nq\')\x81q(}q)h\x0c}q*sbX\x07\x00\x00\x00threadsq+K\x01X\x05\x00\x00\x00inputq,csnakemake.io\nInputFiles\nq-)\x81q.X\x1d\x00\x00\x00sratofastq/SRR970588.fastq.gzq/a}q0(X\x02\x00\x00\x00R1q1h/h\x0c}q2h1K\x00N\x86q3subX\x04\x00\x00\x00ruleq4X\x10\x00\x00\x00perform_trimmingq5X\t\x00\x00\x00resourcesq6csnakemake.io\nResources\nq7)\x81q8(K\x01K\x01e}q9(X\x06\x00\x00\x00_coresq:K\x01h\x0c}q;(h:K\x00N\x86q<X\x06\x00\x00\x00_nodesq=K\x01N\x86q>uh=K\x01ubX\t\x00\x00\x00wildcardsq?csnakemake.io\nWildcards\nq@)\x81qAX\t\x00\x00\x00SRR970588qBa}qC(X\x06\x00\x00\x00sampleqDhBh\x0c}qEX\x06\x00\x00\x00sampleqFK\x00N\x86qGsubX\x06\x00\x00\x00configqH}qIX\x0b\x00\x00\x00config_pathqJX\x19\x00\x00\x00configs/hg38_SRP029589.pyqKsub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
