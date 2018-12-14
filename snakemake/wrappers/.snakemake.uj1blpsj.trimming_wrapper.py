######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05X\x1e\x00\x00\x00sratofastq/SRR5680922.fastq.gzq\x06a}q\x07(X\x06\x00\x00\x00_namesq\x08}q\tX\x02\x00\x00\x00R1q\nK\x00N\x86q\x0bsh\nh\x06ubX\x07\x00\x00\x00threadsq\x0cK\x01X\t\x00\x00\x00resourcesq\rcsnakemake.io\nResources\nq\x0e)\x81q\x0f(K\x01K\x01e}q\x10(h\x08}q\x11(X\x06\x00\x00\x00_nodesq\x12K\x00N\x86q\x13X\x06\x00\x00\x00_coresq\x14K\x01N\x86q\x15uh\x12K\x01h\x14K\x01ubX\x04\x00\x00\x00ruleq\x16X\x10\x00\x00\x00perform_trimmingq\x17X\x06\x00\x00\x00paramsq\x18csnakemake.io\nParams\nq\x19)\x81q\x1a(K&K\x05X\x07\x00\x00\x00defaultq\x1bX\x12\x00\x00\x00preprocessed_step1q\x1cK\x12X\x0c\x00\x00\x00preprocessedq\x1de}q\x1e(X\t\x00\x00\x00pass2_dirq\x1fh\x1dh\x08}q (h\x1fK\x05N\x86q!X\x0c\x00\x00\x00phred_cutoffq"K\x01N\x86q#X\x07\x00\x00\x00adapterq$K\x02N\x86q%X\n\x00\x00\x00max_lengthq&K\x00N\x86q\'X\n\x00\x00\x00min_lengthq(K\x04N\x86q)X\t\x00\x00\x00pass1_dirq*K\x03N\x86q+uh"K\x05h$h\x1bh&K&h(K\x12h*h\x1cubX\x06\x00\x00\x00configq,}q-X\x0b\x00\x00\x00config_pathq.X\x19\x00\x00\x00configs/hg38_SRP109126.pyq/sX\x03\x00\x00\x00logq0csnakemake.io\nLog\nq1)\x81q2}q3h\x08}q4sbX\x06\x00\x00\x00outputq5csnakemake.io\nOutputFiles\nq6)\x81q7(X-\x00\x00\x00preprocessed/SRR5680922_trimmed_trimmed.fq.gzq8X:\x00\x00\x00preprocessed_step1/SRR5680922.fastq.gz_trimming_report.txtq9X+\x00\x00\x00preprocessed_step1/SRR5680922_trimmed.fq.gzq:e}q;(X\x08\x00\x00\x00pass2_fqq<h8X\x0f\x00\x00\x00pass1_fq_reportq=h9h\x08}q>(h<K\x00N\x86q?h=K\x01N\x86q@X\x08\x00\x00\x00pass1_fqqAK\x02N\x86qBuhAh:ubX\t\x00\x00\x00wildcardsqCcsnakemake.io\nWildcards\nqD)\x81qEX\n\x00\x00\x00SRR5680922qFa}qG(h\x08}qHX\x06\x00\x00\x00sampleqIK\x00N\x86qJsX\x06\x00\x00\x00sampleqKhFubub.'
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
