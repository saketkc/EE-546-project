
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\t\x00\x00\x00resourcesq\x03csnakemake.io\nResources\nq\x04)\x81q\x05(K\x01K\x01e}q\x06(X\x06\x00\x00\x00_namesq\x07}q\x08(X\x06\x00\x00\x00_nodesq\tK\x00N\x86q\nX\x06\x00\x00\x00_coresq\x0bK\x01N\x86q\x0cuh\tK\x01h\x0bK\x01ubX\x06\x00\x00\x00paramsq\rcsnakemake.io\nParams\nq\x0e)\x81q\x0f(X\x12\x00\x00\x00preprocessed_step1q\x10K\x12K\x05K&X\x0c\x00\x00\x00preprocessedq\x11X\x07\x00\x00\x00defaultq\x12e}q\x13(X\t\x00\x00\x00pass1_dirq\x14h\x10X\n\x00\x00\x00min_lengthq\x15K\x12X\x0c\x00\x00\x00phred_cutoffq\x16K\x05h\x07}q\x17(h\x14K\x00N\x86q\x18h\x15K\x01N\x86q\x19h\x16K\x02N\x86q\x1aX\n\x00\x00\x00max_lengthq\x1bK\x03N\x86q\x1cX\t\x00\x00\x00pass2_dirq\x1dK\x04N\x86q\x1eX\x07\x00\x00\x00adapterq\x1fK\x05N\x86q uh\x1bK&h\x1dh\x11h\x1fh\x12ubX\x05\x00\x00\x00inputq!csnakemake.io\nInputFiles\nq")\x81q#X\x1d\x00\x00\x00sratofastq/SRR970543.fastq.gzq$a}q%(h\x07}q&X\x02\x00\x00\x00R1q\'K\x00N\x86q(sh\'h$ubX\x06\x00\x00\x00configq)}q*X\x0b\x00\x00\x00config_pathq+X\x19\x00\x00\x00configs/hg38_SRP029589.pyq,sX\x03\x00\x00\x00logq-csnakemake.io\nLog\nq.)\x81q/}q0h\x07}q1sbX\t\x00\x00\x00wildcardsq2csnakemake.io\nWildcards\nq3)\x81q4X\t\x00\x00\x00SRR970543q5a}q6(h\x07}q7X\x06\x00\x00\x00sampleq8K\x00N\x86q9sX\x06\x00\x00\x00sampleq:h5ubX\x07\x00\x00\x00threadsq;K\x01X\x04\x00\x00\x00ruleq<X\x10\x00\x00\x00perform_trimmingq=X\x06\x00\x00\x00outputq>csnakemake.io\nOutputFiles\nq?)\x81q@(X,\x00\x00\x00preprocessed/SRR970543_trimmed_trimmed.fq.gzqAX*\x00\x00\x00preprocessed_step1/SRR970543_trimmed.fq.gzqBX9\x00\x00\x00preprocessed_step1/SRR970543.fastq.gz_trimming_report.txtqCe}qD(h\x07}qE(X\x08\x00\x00\x00pass2_fqqFK\x00N\x86qGX\x08\x00\x00\x00pass1_fqqHK\x01N\x86qIX\x0f\x00\x00\x00pass1_fq_reportqJK\x02N\x86qKuhFhAhHhBhJhCubub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
