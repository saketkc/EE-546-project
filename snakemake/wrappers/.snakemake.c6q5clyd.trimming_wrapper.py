######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00paramsq\x03csnakemake.io\nParams\nq\x04)\x81q\x05(X\x0c\x00\x00\x00preprocessedq\x06X\x12\x00\x00\x00preprocessed_step1q\x07K\x12X\x07\x00\x00\x00defaultq\x08K&K\x05e}q\t(X\t\x00\x00\x00pass2_dirq\nh\x06X\t\x00\x00\x00pass1_dirq\x0bh\x07X\x0c\x00\x00\x00phred_cutoffq\x0cK\x05X\x07\x00\x00\x00adapterq\rh\x08X\x06\x00\x00\x00_namesq\x0e}q\x0f(h\nK\x00N\x86q\x10h\x0bK\x01N\x86q\x11h\x0cK\x05N\x86q\x12h\rK\x03N\x86q\x13X\n\x00\x00\x00max_lengthq\x14K\x04N\x86q\x15X\n\x00\x00\x00min_lengthq\x16K\x02N\x86q\x17uh\x14K&h\x16K\x12ubX\t\x00\x00\x00resourcesq\x18csnakemake.io\nResources\nq\x19)\x81q\x1a(K\x01K\x01e}q\x1b(X\x06\x00\x00\x00_coresq\x1cK\x01X\x06\x00\x00\x00_nodesq\x1dK\x01h\x0e}q\x1e(h\x1cK\x00N\x86q\x1fh\x1dK\x01N\x86q uubX\x07\x00\x00\x00threadsq!K\x01X\x05\x00\x00\x00inputq"csnakemake.io\nInputFiles\nq#)\x81q$X\x1e\x00\x00\x00sratofastq/SRR5680922.fastq.gzq%a}q&(X\x02\x00\x00\x00R1q\'h%h\x0e}q(h\'K\x00N\x86q)subX\x06\x00\x00\x00outputq*csnakemake.io\nOutputFiles\nq+)\x81q,(X+\x00\x00\x00preprocessed_step1/SRR5680922_trimmed.fq.gzq-X:\x00\x00\x00preprocessed_step1/SRR5680922.fastq.gz_trimming_report.txtq.X-\x00\x00\x00preprocessed/SRR5680922_trimmed_trimmed.fq.gzq/e}q0(X\x08\x00\x00\x00pass1_fqq1h-X\x0f\x00\x00\x00pass1_fq_reportq2h.X\x08\x00\x00\x00pass2_fqq3h/h\x0e}q4(h1K\x00N\x86q5h2K\x01N\x86q6h3K\x02N\x86q7uubX\x03\x00\x00\x00logq8csnakemake.io\nLog\nq9)\x81q:}q;h\x0e}q<sbX\x06\x00\x00\x00configq=}q>X\x0b\x00\x00\x00config_pathq?X\x19\x00\x00\x00configs/hg38_SRP109126.pyq@sX\x04\x00\x00\x00ruleqAX\x10\x00\x00\x00perform_trimmingqBX\t\x00\x00\x00wildcardsqCcsnakemake.io\nWildcards\nqD)\x81qEX\n\x00\x00\x00SRR5680922qFa}qG(X\x06\x00\x00\x00sampleqHhFh\x0e}qIX\x06\x00\x00\x00sampleqJK\x00N\x86qKsubub.'
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
