######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05X\x1e\x00\x00\x00sratofastq/SRR5680920.fastq.gzq\x06a}q\x07(X\x02\x00\x00\x00R1q\x08h\x06X\x06\x00\x00\x00_namesq\t}q\nh\x08K\x00N\x86q\x0bsubX\x06\x00\x00\x00configq\x0c}q\rX\x0b\x00\x00\x00config_pathq\x0eX\x19\x00\x00\x00configs/hg38_SRP109126.pyq\x0fsX\t\x00\x00\x00wildcardsq\x10csnakemake.io\nWildcards\nq\x11)\x81q\x12X\n\x00\x00\x00SRR5680920q\x13a}q\x14(X\x06\x00\x00\x00sampleq\x15h\x13h\t}q\x16X\x06\x00\x00\x00sampleq\x17K\x00N\x86q\x18subX\x07\x00\x00\x00threadsq\x19K\x01X\x03\x00\x00\x00logq\x1acsnakemake.io\nLog\nq\x1b)\x81q\x1c}q\x1dh\t}q\x1esbX\x06\x00\x00\x00outputq\x1fcsnakemake.io\nOutputFiles\nq )\x81q!(X-\x00\x00\x00preprocessed/SRR5680920_trimmed_trimmed.fq.gzq"X+\x00\x00\x00preprocessed_step1/SRR5680920_trimmed.fq.gzq#X:\x00\x00\x00preprocessed_step1/SRR5680920.fastq.gz_trimming_report.txtq$e}q%(X\x08\x00\x00\x00pass2_fqq&h"X\x08\x00\x00\x00pass1_fqq\'h#X\x0f\x00\x00\x00pass1_fq_reportq(h$h\t}q)(h&K\x00N\x86q*h\'K\x01N\x86q+h(K\x02N\x86q,uubX\x06\x00\x00\x00paramsq-csnakemake.io\nParams\nq.)\x81q/(K\x12X\x12\x00\x00\x00preprocessed_step1q0K&K\x05X\x0c\x00\x00\x00preprocessedq1X\x07\x00\x00\x00defaultq2e}q3(X\n\x00\x00\x00min_lengthq4K\x12X\t\x00\x00\x00pass1_dirq5h0X\n\x00\x00\x00max_lengthq6K&X\x0c\x00\x00\x00phred_cutoffq7K\x05X\t\x00\x00\x00pass2_dirq8h1X\x07\x00\x00\x00adapterq9h2h\t}q:(h4K\x00N\x86q;h5K\x01N\x86q<h6K\x02N\x86q=h7K\x03N\x86q>h8K\x04N\x86q?h9K\x05N\x86q@uubX\t\x00\x00\x00resourcesqAcsnakemake.io\nResources\nqB)\x81qC(K\x01K\x01e}qD(X\x06\x00\x00\x00_coresqEK\x01X\x06\x00\x00\x00_nodesqFK\x01h\t}qG(hFK\x00N\x86qHhEK\x01N\x86qIuubX\x04\x00\x00\x00ruleqJX\x10\x00\x00\x00perform_trimmingqKub.'
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
