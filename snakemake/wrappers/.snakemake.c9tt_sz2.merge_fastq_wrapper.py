######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\t\x00\x00\x00resourcesq\x03csnakemake.io\nResources\nq\x04)\x81q\x05(K\x01K\x01e}q\x06(X\x06\x00\x00\x00_nodesq\x07K\x01X\x06\x00\x00\x00_coresq\x08K\x01X\x06\x00\x00\x00_namesq\t}q\n(h\x07K\x00N\x86q\x0bh\x08K\x01N\x86q\x0cuubX\x04\x00\x00\x00ruleq\rX\x0b\x00\x00\x00merge_fastqq\x0eX\x03\x00\x00\x00logq\x0fcsnakemake.io\nLog\nq\x10)\x81q\x11}q\x12h\t}q\x13sbX\x07\x00\x00\x00threadsq\x14K\x01X\t\x00\x00\x00wildcardsq\x15csnakemake.io\nWildcards\nq\x16)\x81q\x17X\n\x00\x00\x00SRX2536403q\x18a}q\x19(X\x06\x00\x00\x00sampleq\x1ah\x18h\t}q\x1bX\x06\x00\x00\x00sampleq\x1cK\x00N\x86q\x1dsubX\x05\x00\x00\x00inputq\x1ecsnakemake.io\nInputFiles\nq\x1f)\x81q (X\x1e\x00\x00\x00sratofastq/SRR5227288.fastq.gzq!X\x1e\x00\x00\x00sratofastq/SRR5227288.fastq.gzq"e}q#(X\r\x00\x00\x00dynamic_inputq$csnakemake.io\nNamedlist\nq%)\x81q&h!a}q\'h\t}q(sbh\t}q)(h$K\x00K\x01\x86q*X\t\x00\x00\x00all_fastqq+K\x01K\x02\x86q,uh+h%)\x81q-h"a}q.h\t}q/sbubX\x06\x00\x00\x00outputq0csnakemake.io\nOutputFiles\nq1)\x81q2X \x00\x00\x00fastq_merged/SRX2536403.fastq.gzq3a}q4h\t}q5sbX\x06\x00\x00\x00paramsq6csnakemake.io\nParams\nq7)\x81q8}q9h\t}q:sbX\x06\x00\x00\x00configq;}q<X\x0b\x00\x00\x00config_pathq=X\x19\x00\x00\x00configs/hg38_SRP098789.pyq>sub.'
)
from snakemake.logging import logger
logger.printshellcmds = True
######## Original script #########
import os
from snakemake.shell import shell

input = (' ').join(snakemake.input.dynamic_input)

shell(r'''cat {input} > {snakemake.output}''')
