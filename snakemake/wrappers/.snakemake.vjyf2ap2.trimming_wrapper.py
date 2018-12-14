
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05X\x1d\x00\x00\x00sratofastq/SRR970550.fastq.gzq\x06a}q\x07(X\x06\x00\x00\x00_namesq\x08}q\tX\x02\x00\x00\x00R1q\nK\x00N\x86q\x0bsh\nh\x06ubX\x06\x00\x00\x00configq\x0c}q\rX\x0b\x00\x00\x00config_pathq\x0eX\x19\x00\x00\x00configs/hg38_SRP029589.pyq\x0fsX\x06\x00\x00\x00outputq\x10csnakemake.io\nOutputFiles\nq\x11)\x81q\x12(X*\x00\x00\x00preprocessed_step1/SRR970550_trimmed.fq.gzq\x13X9\x00\x00\x00preprocessed_step1/SRR970550.fastq.gz_trimming_report.txtq\x14X,\x00\x00\x00preprocessed/SRR970550_trimmed_trimmed.fq.gzq\x15e}q\x16(h\x08}q\x17(X\x08\x00\x00\x00pass1_fqq\x18K\x00N\x86q\x19X\x0f\x00\x00\x00pass1_fq_reportq\x1aK\x01N\x86q\x1bX\x08\x00\x00\x00pass2_fqq\x1cK\x02N\x86q\x1duh\x18h\x13h\x1ah\x14h\x1ch\x15ubX\t\x00\x00\x00resourcesq\x1ecsnakemake.io\nResources\nq\x1f)\x81q (K\x01K\x01e}q!(h\x08}q"(X\x06\x00\x00\x00_coresq#K\x00N\x86q$X\x06\x00\x00\x00_nodesq%K\x01N\x86q&uh#K\x01h%K\x01ubX\x03\x00\x00\x00logq\'csnakemake.io\nLog\nq()\x81q)}q*h\x08}q+sbX\x06\x00\x00\x00paramsq,csnakemake.io\nParams\nq-)\x81q.(K\x12X\x07\x00\x00\x00defaultq/K&X\x0c\x00\x00\x00preprocessedq0K\x05X\x12\x00\x00\x00preprocessed_step1q1e}q2(X\n\x00\x00\x00min_lengthq3K\x12X\x07\x00\x00\x00adapterq4h/X\t\x00\x00\x00pass1_dirq5h1h\x08}q6(h3K\x00N\x86q7h4K\x01N\x86q8h5K\x05N\x86q9X\t\x00\x00\x00pass2_dirq:K\x03N\x86q;X\x0c\x00\x00\x00phred_cutoffq<K\x04N\x86q=X\n\x00\x00\x00max_lengthq>K\x02N\x86q?uh:h0h<K\x05h>K&ubX\x07\x00\x00\x00threadsq@K\x01X\x04\x00\x00\x00ruleqAX\x10\x00\x00\x00perform_trimmingqBX\t\x00\x00\x00wildcardsqCcsnakemake.io\nWildcards\nqD)\x81qEX\t\x00\x00\x00SRR970550qFa}qG(h\x08}qHX\x06\x00\x00\x00sampleqIK\x00N\x86qJsX\x06\x00\x00\x00sampleqKhFubub.'); from snakemake.logging import logger; logger.printshellcmds = True
######## Original script #########
from snakemake.shell import shell
from riboraptor.kmer import fastq_kmer_histogram
import operator

# CTG is a common one
# AGATCG.. is a Truseq ribo 3' illumina
# TGGAAT.. is Truseq small rna

PREFERRED_KMERS = [
    'CTGTAGGCACCATCAAT', 'AGATCGGAAGAGCACACGTCT', 'TGGAATTCTCGGGTGCCAAGG',
    'CTGTAGGCAC'
]


def get_top_kmer(kmer_series_dict):
    """Return a kmer if it's percentage exceed 30%"""
    # Start from the longest kmer and stop
    # at where this criterion is met
    for kmer_length, kmer_series in sorted(
            kmer_series_dict.items(), key=operator.itemgetter(0),
            reverse=True):
        kmer_list = kmer_series.index.tolist()
        # Are any of the top 4 kmers from our PREFFERED_KMERS list?
        index_counter = 0
        for kmer in kmer_list:
            index_counter += 1
            if kmer in PREFERRED_KMERS:
                return kmer
            if index_counter >= 5:
                break
        # 30 is an arbitrary chosen threshold, but it works in *most* cases
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
