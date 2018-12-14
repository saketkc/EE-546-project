
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00outputq\x03csnakemake.io\nOutputFiles\nq\x04)\x81q\x05(X*\x00\x00\x00preprocessed_step1/SRR970543_trimmed.fq.gzq\x06X9\x00\x00\x00preprocessed_step1/SRR970543.fastq.gz_trimming_report.txtq\x07X,\x00\x00\x00preprocessed/SRR970543_trimmed_trimmed.fq.gzq\x08e}q\t(X\x08\x00\x00\x00pass1_fqq\nh\x06X\x0f\x00\x00\x00pass1_fq_reportq\x0bh\x07X\x06\x00\x00\x00_namesq\x0c}q\r(h\nK\x00N\x86q\x0eh\x0bK\x01N\x86q\x0fX\x08\x00\x00\x00pass2_fqq\x10K\x02N\x86q\x11uh\x10h\x08ubX\x07\x00\x00\x00threadsq\x12K\x01X\t\x00\x00\x00resourcesq\x13csnakemake.io\nResources\nq\x14)\x81q\x15(K\x01K\x01e}q\x16(X\x06\x00\x00\x00_coresq\x17K\x01h\x0c}q\x18(h\x17K\x00N\x86q\x19X\x06\x00\x00\x00_nodesq\x1aK\x01N\x86q\x1buh\x1aK\x01ubX\x06\x00\x00\x00configq\x1c}q\x1dX\x0b\x00\x00\x00config_pathq\x1eX\x19\x00\x00\x00configs/hg38_SRP029589.pyq\x1fsX\x06\x00\x00\x00paramsq csnakemake.io\nParams\nq!)\x81q"(K&K\x12X\x12\x00\x00\x00preprocessed_step1q#X\x0c\x00\x00\x00preprocessedq$K\x05X\x07\x00\x00\x00defaultq%e}q&(X\n\x00\x00\x00max_lengthq\'K&X\n\x00\x00\x00min_lengthq(K\x12X\t\x00\x00\x00pass1_dirq)h#X\x07\x00\x00\x00adapterq*h%X\x0c\x00\x00\x00phred_cutoffq+K\x05h\x0c}q,(h\'K\x00N\x86q-h(K\x01N\x86q.h)K\x02N\x86q/h*K\x05N\x86q0h+K\x04N\x86q1X\t\x00\x00\x00pass2_dirq2K\x03N\x86q3uh2h$ubX\x05\x00\x00\x00inputq4csnakemake.io\nInputFiles\nq5)\x81q6X\x1d\x00\x00\x00sratofastq/SRR970543.fastq.gzq7a}q8(X\x02\x00\x00\x00R1q9h7h\x0c}q:h9K\x00N\x86q;subX\x04\x00\x00\x00ruleq<X\x10\x00\x00\x00perform_trimmingq=X\x03\x00\x00\x00logq>csnakemake.io\nLog\nq?)\x81q@}qAh\x0c}qBsbX\t\x00\x00\x00wildcardsqCcsnakemake.io\nWildcards\nqD)\x81qEX\t\x00\x00\x00SRR970543qFa}qG(X\x06\x00\x00\x00sampleqHhFh\x0c}qIX\x06\x00\x00\x00sampleqJK\x00N\x86qKsubub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
