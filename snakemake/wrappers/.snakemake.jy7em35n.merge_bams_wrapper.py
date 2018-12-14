######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x07\x00\x00\x00threadsq\x03K\x08X\t\x00\x00\x00resourcesq\x04csnakemake.io\nResources\nq\x05)\x81q\x06(K\x01K\x08e}q\x07(X\x06\x00\x00\x00_nodesq\x08K\x01X\x06\x00\x00\x00_namesq\t}q\n(h\x08K\x00N\x86q\x0bX\x06\x00\x00\x00_coresq\x0cK\x01N\x86q\ruh\x0cK\x08ubX\x05\x00\x00\x00inputq\x0ecsnakemake.io\nInputFiles\nq\x0f)\x81q\x10X\x17\x00\x00\x00bams_srr/SRR5227288.bamq\x11a}q\x12h\t}q\x13sbX\x03\x00\x00\x00logq\x14csnakemake.io\nLog\nq\x15)\x81q\x16}q\x17h\t}q\x18sbX\x06\x00\x00\x00configq\x19}q\x1aX\x0b\x00\x00\x00config_pathq\x1bX\x19\x00\x00\x00configs/hg38_SRP098789.pyq\x1csX\t\x00\x00\x00wildcardsq\x1dcsnakemake.io\nWildcards\nq\x1e)\x81q\x1fX\n\x00\x00\x00SRX2536403q a}q!(X\x06\x00\x00\x00sampleq"h h\t}q#X\x06\x00\x00\x00sampleq$K\x00N\x86q%subX\x06\x00\x00\x00outputq&csnakemake.io\nOutputFiles\nq\')\x81q(X\x13\x00\x00\x00bams/SRX2536403.bamq)a}q*h\t}q+sbX\x06\x00\x00\x00paramsq,csnakemake.io\nParams\nq-)\x81q.X\x04\x00\x00\x00/tmpq/a}q0(h\t}q1X\x07\x00\x00\x00tmp_dirq2K\x00N\x86q3sh2h/ubX\x04\x00\x00\x00ruleq4X\n\x00\x00\x00merge_bamsq5ub.'
)
from snakemake.logging import logger
logger.printshellcmds = True
######## Original script #########
import os
import tempfile
from snakemake.shell import shell

if len(snakemake.input) > 1:
    with tempfile.TemporaryDirectory(dir=snakemake.params.tmp_dir) as temp_dir:
        cmd = ' -in '.join(snakemake.input)
        shell(r'''bamtools merge -in {cmd} -out {snakemake.output}.unsorted \
              && samtools sort -@ {snakemake.threads} \
              -T {temp_dir}/{snakemake.wildcards.sample}_merge_bam \
              -o {snakemake.output} {snakemake.output}.unsorted \
              && samtools index {snakemake.output} \
              && yes | rm -rf {snakemake.output}.unsorted''')
elif len(snakemake.input) == 1:
    source = os.path.abspath(str(snakemake.input[0]))
    destination = os.path.abspath(str(snakemake.output))
    shell('''cp {source} {destination} && cp {source}.bai {destination}.bai''')
