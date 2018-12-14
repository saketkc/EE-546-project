######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00paramsq\x03csnakemake.io\nParams\nq\x04)\x81q\x05(K\x05X\x07\x00\x00\x00defaultq\x06K\x12K&X\x0c\x00\x00\x00preprocessedq\x07X\x12\x00\x00\x00preprocessed_step1q\x08e}q\t(X\x0c\x00\x00\x00phred_cutoffq\nK\x05X\x07\x00\x00\x00adapterq\x0bh\x06X\n\x00\x00\x00min_lengthq\x0cK\x12X\t\x00\x00\x00pass2_dirq\rh\x07X\t\x00\x00\x00pass1_dirq\x0eh\x08X\x06\x00\x00\x00_namesq\x0f}q\x10(h\nK\x00N\x86q\x11h\x0bK\x01N\x86q\x12h\x0cK\x02N\x86q\x13h\x0eK\x05N\x86q\x14h\rK\x04N\x86q\x15X\n\x00\x00\x00max_lengthq\x16K\x03N\x86q\x17uh\x16K&ubX\x03\x00\x00\x00logq\x18csnakemake.io\nLog\nq\x19)\x81q\x1a}q\x1bh\x0f}q\x1csbX\x04\x00\x00\x00ruleq\x1dX\x10\x00\x00\x00perform_trimmingq\x1eX\x06\x00\x00\x00configq\x1f}q X\x0b\x00\x00\x00config_pathq!X\x19\x00\x00\x00configs/hg38_SRP109126.pyq"sX\x06\x00\x00\x00outputq#csnakemake.io\nOutputFiles\nq$)\x81q%(X-\x00\x00\x00preprocessed/SRR5680920_trimmed_trimmed.fq.gzq&X+\x00\x00\x00preprocessed_step1/SRR5680920_trimmed.fq.gzq\'X:\x00\x00\x00preprocessed_step1/SRR5680920.fastq.gz_trimming_report.txtq(e}q)(X\x08\x00\x00\x00pass2_fqq*h&X\x08\x00\x00\x00pass1_fqq+h\'h\x0f}q,(h*K\x00N\x86q-h+K\x01N\x86q.X\x0f\x00\x00\x00pass1_fq_reportq/K\x02N\x86q0uh/h(ubX\x07\x00\x00\x00threadsq1K\x01X\x05\x00\x00\x00inputq2csnakemake.io\nInputFiles\nq3)\x81q4X\x1e\x00\x00\x00sratofastq/SRR5680920.fastq.gzq5a}q6(X\x02\x00\x00\x00R1q7h5h\x0f}q8h7K\x00N\x86q9subX\t\x00\x00\x00resourcesq:csnakemake.io\nResources\nq;)\x81q<(K\x01K\x01e}q=(h\x0f}q>(X\x06\x00\x00\x00_nodesq?K\x00N\x86q@X\x06\x00\x00\x00_coresqAK\x01N\x86qBuh?K\x01hAK\x01ubX\t\x00\x00\x00wildcardsqCcsnakemake.io\nWildcards\nqD)\x81qEX\n\x00\x00\x00SRR5680920qFa}qG(h\x0f}qHX\x06\x00\x00\x00sampleqIK\x00N\x86qJsX\x06\x00\x00\x00sampleqKhFubub.'
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
