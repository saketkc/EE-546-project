######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x03\x00\x00\x00logq\x03csnakemake.io\nLog\nq\x04)\x81q\x05}q\x06X\x06\x00\x00\x00_namesq\x07}q\x08sbX\x06\x00\x00\x00outputq\tcsnakemake.io\nOutputFiles\nq\n)\x81q\x0b(X+\x00\x00\x00preprocessed_step1/SRR5680923_trimmed.fq.gzq\x0cX-\x00\x00\x00preprocessed/SRR5680923_trimmed_trimmed.fq.gzq\rX:\x00\x00\x00preprocessed_step1/SRR5680923.fastq.gz_trimming_report.txtq\x0ee}q\x0f(X\x08\x00\x00\x00pass1_fqq\x10h\x0cX\x08\x00\x00\x00pass2_fqq\x11h\rh\x07}q\x12(h\x10K\x00N\x86q\x13h\x11K\x01N\x86q\x14X\x0f\x00\x00\x00pass1_fq_reportq\x15K\x02N\x86q\x16uh\x15h\x0eubX\x06\x00\x00\x00paramsq\x17csnakemake.io\nParams\nq\x18)\x81q\x19(X\x0c\x00\x00\x00preprocessedq\x1aK&K\x05X\x07\x00\x00\x00defaultq\x1bK\x12X\x12\x00\x00\x00preprocessed_step1q\x1ce}q\x1d(X\t\x00\x00\x00pass2_dirq\x1eh\x1aX\t\x00\x00\x00pass1_dirq\x1fh\x1ch\x07}q (h\x1eK\x00N\x86q!h\x1fK\x05N\x86q"X\x0c\x00\x00\x00phred_cutoffq#K\x02N\x86q$X\x07\x00\x00\x00adapterq%K\x03N\x86q&X\n\x00\x00\x00min_lengthq\'K\x04N\x86q(X\n\x00\x00\x00max_lengthq)K\x01N\x86q*uh#K\x05h%h\x1bh\'K\x12h)K&ubX\x04\x00\x00\x00ruleq+X\x10\x00\x00\x00perform_trimmingq,X\t\x00\x00\x00wildcardsq-csnakemake.io\nWildcards\nq.)\x81q/X\n\x00\x00\x00SRR5680923q0a}q1(X\x06\x00\x00\x00sampleq2h0h\x07}q3X\x06\x00\x00\x00sampleq4K\x00N\x86q5subX\x07\x00\x00\x00threadsq6K\x01X\x05\x00\x00\x00inputq7csnakemake.io\nInputFiles\nq8)\x81q9X\x1e\x00\x00\x00sratofastq/SRR5680923.fastq.gzq:a}q;(X\x02\x00\x00\x00R1q<h:h\x07}q=h<K\x00N\x86q>subX\t\x00\x00\x00resourcesq?csnakemake.io\nResources\nq@)\x81qA(K\x01K\x01e}qB(X\x06\x00\x00\x00_coresqCK\x01h\x07}qD(hCK\x00N\x86qEX\x06\x00\x00\x00_nodesqFK\x01N\x86qGuhFK\x01ubX\x06\x00\x00\x00configqH}qIX\x0b\x00\x00\x00config_pathqJX\x19\x00\x00\x00configs/hg38_SRP109126.pyqKsub.'
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
