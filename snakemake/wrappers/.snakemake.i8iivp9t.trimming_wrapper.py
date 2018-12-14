
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00configq\x03}q\x04X\x0b\x00\x00\x00config_pathq\x05X\x19\x00\x00\x00configs/hg38_SRP029589.pyq\x06sX\t\x00\x00\x00resourcesq\x07csnakemake.io\nResources\nq\x08)\x81q\t(K\x01K\x01e}q\n(X\x06\x00\x00\x00_coresq\x0bK\x01X\x06\x00\x00\x00_nodesq\x0cK\x01X\x06\x00\x00\x00_namesq\r}q\x0e(h\x0bK\x00N\x86q\x0fh\x0cK\x01N\x86q\x10uubX\x04\x00\x00\x00ruleq\x11X\x10\x00\x00\x00perform_trimmingq\x12X\x03\x00\x00\x00logq\x13csnakemake.io\nLog\nq\x14)\x81q\x15}q\x16h\r}q\x17sbX\x06\x00\x00\x00outputq\x18csnakemake.io\nOutputFiles\nq\x19)\x81q\x1a(X,\x00\x00\x00preprocessed/SRR970593_trimmed_trimmed.fq.gzq\x1bX*\x00\x00\x00preprocessed_step1/SRR970593_trimmed.fq.gzq\x1cX9\x00\x00\x00preprocessed_step1/SRR970593.fastq.gz_trimming_report.txtq\x1de}q\x1e(X\x08\x00\x00\x00pass2_fqq\x1fh\x1bX\x08\x00\x00\x00pass1_fqq h\x1cX\x0f\x00\x00\x00pass1_fq_reportq!h\x1dh\r}q"(h\x1fK\x00N\x86q#h K\x01N\x86q$h!K\x02N\x86q%uubX\x06\x00\x00\x00paramsq&csnakemake.io\nParams\nq\')\x81q((K&X\x12\x00\x00\x00preprocessed_step1q)K\x12X\x07\x00\x00\x00defaultq*X\x0c\x00\x00\x00preprocessedq+K\x05e}q,(X\n\x00\x00\x00max_lengthq-K&X\t\x00\x00\x00pass1_dirq.h)X\n\x00\x00\x00min_lengthq/K\x12X\x07\x00\x00\x00adapterq0h*X\t\x00\x00\x00pass2_dirq1h+X\x0c\x00\x00\x00phred_cutoffq2K\x05h\r}q3(h-K\x00N\x86q4h.K\x01N\x86q5h/K\x02N\x86q6h0K\x03N\x86q7h1K\x04N\x86q8h2K\x05N\x86q9uubX\x05\x00\x00\x00inputq:csnakemake.io\nInputFiles\nq;)\x81q<X\x1d\x00\x00\x00sratofastq/SRR970593.fastq.gzq=a}q>(X\x02\x00\x00\x00R1q?h=h\r}q@h?K\x00N\x86qAsubX\x07\x00\x00\x00threadsqBK\x01X\t\x00\x00\x00wildcardsqCcsnakemake.io\nWildcards\nqD)\x81qEX\t\x00\x00\x00SRR970593qFa}qG(X\x06\x00\x00\x00sampleqHhFh\r}qIX\x06\x00\x00\x00sampleqJK\x00N\x86qKsubub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
