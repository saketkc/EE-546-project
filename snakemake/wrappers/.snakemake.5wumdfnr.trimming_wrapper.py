######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x07\x00\x00\x00threadsq\x03K\x01X\x06\x00\x00\x00configq\x04}q\x05X\x0b\x00\x00\x00config_pathq\x06X\x19\x00\x00\x00configs/hg38_SRP109126.pyq\x07sX\x06\x00\x00\x00paramsq\x08csnakemake.io\nParams\nq\t)\x81q\n(K&X\x12\x00\x00\x00preprocessed_step1q\x0bK\x05K\x12X\x0c\x00\x00\x00preprocessedq\x0cX\x07\x00\x00\x00defaultq\re}q\x0e(X\n\x00\x00\x00max_lengthq\x0fK&X\t\x00\x00\x00pass2_dirq\x10h\x0cX\x0c\x00\x00\x00phred_cutoffq\x11K\x05X\n\x00\x00\x00min_lengthq\x12K\x12X\t\x00\x00\x00pass1_dirq\x13h\x0bX\x07\x00\x00\x00adapterq\x14h\rX\x06\x00\x00\x00_namesq\x15}q\x16(h\x0fK\x00N\x86q\x17h\x13K\x01N\x86q\x18h\x11K\x02N\x86q\x19h\x12K\x03N\x86q\x1ah\x10K\x04N\x86q\x1bh\x14K\x05N\x86q\x1cuubX\x03\x00\x00\x00logq\x1dcsnakemake.io\nLog\nq\x1e)\x81q\x1f}q h\x15}q!sbX\t\x00\x00\x00resourcesq"csnakemake.io\nResources\nq#)\x81q$(K\x01K\x01e}q%(X\x06\x00\x00\x00_nodesq&K\x01X\x06\x00\x00\x00_coresq\'K\x01h\x15}q((h\'K\x00N\x86q)h&K\x01N\x86q*uubX\x04\x00\x00\x00ruleq+X\x10\x00\x00\x00perform_trimmingq,X\t\x00\x00\x00wildcardsq-csnakemake.io\nWildcards\nq.)\x81q/X\n\x00\x00\x00SRR5680923q0a}q1(X\x06\x00\x00\x00sampleq2h0h\x15}q3X\x06\x00\x00\x00sampleq4K\x00N\x86q5subX\x06\x00\x00\x00outputq6csnakemake.io\nOutputFiles\nq7)\x81q8(X+\x00\x00\x00preprocessed_step1/SRR5680923_trimmed.fq.gzq9X-\x00\x00\x00preprocessed/SRR5680923_trimmed_trimmed.fq.gzq:X:\x00\x00\x00preprocessed_step1/SRR5680923.fastq.gz_trimming_report.txtq;e}q<(X\x08\x00\x00\x00pass1_fqq=h9X\x0f\x00\x00\x00pass1_fq_reportq>h;X\x08\x00\x00\x00pass2_fqq?h:h\x15}q@(h=K\x00N\x86qAh?K\x01N\x86qBh>K\x02N\x86qCuubX\x05\x00\x00\x00inputqDcsnakemake.io\nInputFiles\nqE)\x81qFX\x1e\x00\x00\x00sratofastq/SRR5680923.fastq.gzqGa}qH(X\x02\x00\x00\x00R1qIhGh\x15}qJhIK\x00N\x86qKsubub.'
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
