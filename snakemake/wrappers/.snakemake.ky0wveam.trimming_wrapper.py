######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x04\x00\x00\x00ruleq\x03X\x10\x00\x00\x00perform_trimmingq\x04X\t\x00\x00\x00wildcardsq\x05csnakemake.io\nWildcards\nq\x06)\x81q\x07X\n\x00\x00\x00SRR5680923q\x08a}q\t(X\x06\x00\x00\x00sampleq\nh\x08X\x06\x00\x00\x00_namesq\x0b}q\x0cX\x06\x00\x00\x00sampleq\rK\x00N\x86q\x0esubX\x05\x00\x00\x00inputq\x0fcsnakemake.io\nInputFiles\nq\x10)\x81q\x11X\x1e\x00\x00\x00sratofastq/SRR5680923.fastq.gzq\x12a}q\x13(X\x02\x00\x00\x00R1q\x14h\x12h\x0b}q\x15h\x14K\x00N\x86q\x16subX\t\x00\x00\x00resourcesq\x17csnakemake.io\nResources\nq\x18)\x81q\x19(K\x01K\x01e}q\x1a(X\x06\x00\x00\x00_coresq\x1bK\x01X\x06\x00\x00\x00_nodesq\x1cK\x01h\x0b}q\x1d(h\x1bK\x00N\x86q\x1eh\x1cK\x01N\x86q\x1fuubX\x07\x00\x00\x00threadsq K\x01X\x06\x00\x00\x00paramsq!csnakemake.io\nParams\nq")\x81q#(K\x12X\x07\x00\x00\x00defaultq$K&K\x05X\x0c\x00\x00\x00preprocessedq%X\x12\x00\x00\x00preprocessed_step1q&e}q\'(X\n\x00\x00\x00min_lengthq(K\x12X\x07\x00\x00\x00adapterq)h$X\n\x00\x00\x00max_lengthq*K&X\x0c\x00\x00\x00phred_cutoffq+K\x05X\t\x00\x00\x00pass2_dirq,h%h\x0b}q-(h(K\x00N\x86q.h)K\x01N\x86q/h*K\x02N\x86q0h+K\x03N\x86q1h,K\x04N\x86q2X\t\x00\x00\x00pass1_dirq3K\x05N\x86q4uh3h&ubX\x03\x00\x00\x00logq5csnakemake.io\nLog\nq6)\x81q7}q8h\x0b}q9sbX\x06\x00\x00\x00configq:}q;X\x0b\x00\x00\x00config_pathq<X\x19\x00\x00\x00configs/hg38_SRP109126.pyq=sX\x06\x00\x00\x00outputq>csnakemake.io\nOutputFiles\nq?)\x81q@(X:\x00\x00\x00preprocessed_step1/SRR5680923.fastq.gz_trimming_report.txtqAX+\x00\x00\x00preprocessed_step1/SRR5680923_trimmed.fq.gzqBX-\x00\x00\x00preprocessed/SRR5680923_trimmed_trimmed.fq.gzqCe}qD(X\x0f\x00\x00\x00pass1_fq_reportqEhAX\x08\x00\x00\x00pass2_fqqFhCh\x0b}qG(hEK\x00N\x86qHX\x08\x00\x00\x00pass1_fqqIK\x01N\x86qJhFK\x02N\x86qKuhIhBubub.'
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
