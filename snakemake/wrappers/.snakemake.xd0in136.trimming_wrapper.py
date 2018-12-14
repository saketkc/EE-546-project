######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x04\x00\x00\x00ruleq\x03X\x10\x00\x00\x00perform_trimmingq\x04X\x06\x00\x00\x00paramsq\x05csnakemake.io\nParams\nq\x06)\x81q\x07(K\x05X\x0c\x00\x00\x00preprocessedq\x08K\x12K&X\x12\x00\x00\x00preprocessed_step1q\tX\x07\x00\x00\x00defaultq\ne}q\x0b(X\n\x00\x00\x00max_lengthq\x0cK&X\t\x00\x00\x00pass2_dirq\rh\x08X\n\x00\x00\x00min_lengthq\x0eK\x12X\x0c\x00\x00\x00phred_cutoffq\x0fK\x05X\t\x00\x00\x00pass1_dirq\x10h\tX\x07\x00\x00\x00adapterq\x11h\nX\x06\x00\x00\x00_namesq\x12}q\x13(h\x0cK\x03N\x86q\x14h\rK\x01N\x86q\x15h\x0eK\x02N\x86q\x16h\x0fK\x00N\x86q\x17h\x10K\x04N\x86q\x18h\x11K\x05N\x86q\x19uubX\t\x00\x00\x00wildcardsq\x1acsnakemake.io\nWildcards\nq\x1b)\x81q\x1cX\n\x00\x00\x00SRR5680916q\x1da}q\x1e(X\x06\x00\x00\x00sampleq\x1fh\x1dh\x12}q X\x06\x00\x00\x00sampleq!K\x00N\x86q"subX\x05\x00\x00\x00inputq#csnakemake.io\nInputFiles\nq$)\x81q%X\x1e\x00\x00\x00sratofastq/SRR5680916.fastq.gzq&a}q\'(X\x02\x00\x00\x00R1q(h&h\x12}q)h(K\x00N\x86q*subX\x06\x00\x00\x00configq+}q,X\x0b\x00\x00\x00config_pathq-X\x19\x00\x00\x00configs/hg38_SRP109126.pyq.sX\x06\x00\x00\x00outputq/csnakemake.io\nOutputFiles\nq0)\x81q1(X-\x00\x00\x00preprocessed/SRR5680916_trimmed_trimmed.fq.gzq2X:\x00\x00\x00preprocessed_step1/SRR5680916.fastq.gz_trimming_report.txtq3X+\x00\x00\x00preprocessed_step1/SRR5680916_trimmed.fq.gzq4e}q5(X\x08\x00\x00\x00pass2_fqq6h2X\x0f\x00\x00\x00pass1_fq_reportq7h3X\x08\x00\x00\x00pass1_fqq8h4h\x12}q9(h6K\x00N\x86q:h7K\x01N\x86q;h8K\x02N\x86q<uubX\t\x00\x00\x00resourcesq=csnakemake.io\nResources\nq>)\x81q?(K\x01K\x01e}q@(X\x06\x00\x00\x00_nodesqAK\x01X\x06\x00\x00\x00_coresqBK\x01h\x12}qC(hAK\x00N\x86qDhBK\x01N\x86qEuubX\x07\x00\x00\x00threadsqFK\x01X\x03\x00\x00\x00logqGcsnakemake.io\nLog\nqH)\x81qI}qJh\x12}qKsbub.'
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
