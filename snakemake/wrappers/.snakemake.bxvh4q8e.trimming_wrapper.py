######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00outputq\x03csnakemake.io\nOutputFiles\nq\x04)\x81q\x05(X+\x00\x00\x00preprocessed_step1/SRR5680922_trimmed.fq.gzq\x06X:\x00\x00\x00preprocessed_step1/SRR5680922.fastq.gz_trimming_report.txtq\x07X-\x00\x00\x00preprocessed/SRR5680922_trimmed_trimmed.fq.gzq\x08e}q\t(X\x08\x00\x00\x00pass1_fqq\nh\x06X\x06\x00\x00\x00_namesq\x0b}q\x0c(h\nK\x00N\x86q\rX\x0f\x00\x00\x00pass1_fq_reportq\x0eK\x01N\x86q\x0fX\x08\x00\x00\x00pass2_fqq\x10K\x02N\x86q\x11uh\x0eh\x07h\x10h\x08ubX\x04\x00\x00\x00ruleq\x12X\x10\x00\x00\x00perform_trimmingq\x13X\x06\x00\x00\x00configq\x14}q\x15X\x0b\x00\x00\x00config_pathq\x16X\x19\x00\x00\x00configs/hg38_SRP109126.pyq\x17sX\x03\x00\x00\x00logq\x18csnakemake.io\nLog\nq\x19)\x81q\x1a}q\x1bh\x0b}q\x1csbX\x05\x00\x00\x00inputq\x1dcsnakemake.io\nInputFiles\nq\x1e)\x81q\x1fX\x1e\x00\x00\x00sratofastq/SRR5680922.fastq.gzq a}q!(h\x0b}q"X\x02\x00\x00\x00R1q#K\x00N\x86q$sh#h ubX\x07\x00\x00\x00threadsq%K\x01X\t\x00\x00\x00resourcesq&csnakemake.io\nResources\nq\')\x81q((K\x01K\x01e}q)(h\x0b}q*(X\x06\x00\x00\x00_nodesq+K\x00N\x86q,X\x06\x00\x00\x00_coresq-K\x01N\x86q.uh+K\x01h-K\x01ubX\x06\x00\x00\x00paramsq/csnakemake.io\nParams\nq0)\x81q1(X\x07\x00\x00\x00defaultq2X\x0c\x00\x00\x00preprocessedq3K\x12K\x05K&X\x12\x00\x00\x00preprocessed_step1q4e}q5(X\x07\x00\x00\x00adapterq6h2X\t\x00\x00\x00pass2_dirq7h3X\n\x00\x00\x00min_lengthq8K\x12h\x0b}q9(h6K\x00N\x86q:h7K\x01N\x86q;X\n\x00\x00\x00max_lengthq<K\x04N\x86q=X\x0c\x00\x00\x00phred_cutoffq>K\x03N\x86q?h8K\x02N\x86q@X\t\x00\x00\x00pass1_dirqAK\x05N\x86qBuh>K\x05h<K&hAh4ubX\t\x00\x00\x00wildcardsqCcsnakemake.io\nWildcards\nqD)\x81qEX\n\x00\x00\x00SRR5680922qFa}qG(h\x0b}qHX\x06\x00\x00\x00sampleqIK\x00N\x86qJsX\x06\x00\x00\x00sampleqKhFubub.'
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
