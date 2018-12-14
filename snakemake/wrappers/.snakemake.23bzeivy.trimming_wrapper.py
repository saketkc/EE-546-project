######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00outputq\x03csnakemake.io\nOutputFiles\nq\x04)\x81q\x05(X:\x00\x00\x00preprocessed_step1/SRR5680922.fastq.gz_trimming_report.txtq\x06X-\x00\x00\x00preprocessed/SRR5680922_trimmed_trimmed.fq.gzq\x07X+\x00\x00\x00preprocessed_step1/SRR5680922_trimmed.fq.gzq\x08e}q\t(X\x0f\x00\x00\x00pass1_fq_reportq\nh\x06X\x08\x00\x00\x00pass2_fqq\x0bh\x07X\x06\x00\x00\x00_namesq\x0c}q\r(h\nK\x00N\x86q\x0eh\x0bK\x01N\x86q\x0fX\x08\x00\x00\x00pass1_fqq\x10K\x02N\x86q\x11uh\x10h\x08ubX\x04\x00\x00\x00ruleq\x12X\x10\x00\x00\x00perform_trimmingq\x13X\t\x00\x00\x00wildcardsq\x14csnakemake.io\nWildcards\nq\x15)\x81q\x16X\n\x00\x00\x00SRR5680922q\x17a}q\x18(X\x06\x00\x00\x00sampleq\x19h\x17h\x0c}q\x1aX\x06\x00\x00\x00sampleq\x1bK\x00N\x86q\x1csubX\x06\x00\x00\x00paramsq\x1dcsnakemake.io\nParams\nq\x1e)\x81q\x1f(K\x05X\x12\x00\x00\x00preprocessed_step1q K&X\x0c\x00\x00\x00preprocessedq!X\x07\x00\x00\x00defaultq"K\x12e}q#(X\x0c\x00\x00\x00phred_cutoffq$K\x05X\x07\x00\x00\x00adapterq%h"X\n\x00\x00\x00max_lengthq&K&X\t\x00\x00\x00pass2_dirq\'h!X\t\x00\x00\x00pass1_dirq(h X\n\x00\x00\x00min_lengthq)K\x12h\x0c}q*(h$K\x00N\x86q+h(K\x01N\x86q,h&K\x02N\x86q-h\'K\x03N\x86q.h%K\x04N\x86q/h)K\x05N\x86q0uubX\x05\x00\x00\x00inputq1csnakemake.io\nInputFiles\nq2)\x81q3X\x1e\x00\x00\x00sratofastq/SRR5680922.fastq.gzq4a}q5(X\x02\x00\x00\x00R1q6h4h\x0c}q7h6K\x00N\x86q8subX\x06\x00\x00\x00configq9}q:X\x0b\x00\x00\x00config_pathq;X\x19\x00\x00\x00configs/hg38_SRP109126.pyq<sX\x03\x00\x00\x00logq=csnakemake.io\nLog\nq>)\x81q?}q@h\x0c}qAsbX\x07\x00\x00\x00threadsqBK\x01X\t\x00\x00\x00resourcesqCcsnakemake.io\nResources\nqD)\x81qE(K\x01K\x01e}qF(X\x06\x00\x00\x00_coresqGK\x01h\x0c}qH(hGK\x00N\x86qIX\x06\x00\x00\x00_nodesqJK\x01N\x86qKuhJK\x01ubub.'
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
