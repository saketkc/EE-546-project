######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\t\x00\x00\x00wildcardsq\x03csnakemake.io\nWildcards\nq\x04)\x81q\x05X\n\x00\x00\x00SRX2536403q\x06a}q\x07(X\x06\x00\x00\x00_namesq\x08}q\tX\x06\x00\x00\x00sampleq\nK\x00N\x86q\x0bsX\x06\x00\x00\x00sampleq\x0ch\x06ubX\x06\x00\x00\x00paramsq\rcsnakemake.io\nParams\nq\x0e)\x81q\x0fX\x04\x00\x00\x00/tmpq\x10a}q\x11(h\x08}q\x12X\x07\x00\x00\x00tmp_dirq\x13K\x00N\x86q\x14sh\x13h\x10ubX\x06\x00\x00\x00configq\x15}q\x16X\x0b\x00\x00\x00config_pathq\x17X\x19\x00\x00\x00configs/hg38_SRP098789.pyq\x18sX\x06\x00\x00\x00outputq\x19csnakemake.io\nOutputFiles\nq\x1a)\x81q\x1bX\x13\x00\x00\x00bams/SRX2536403.bamq\x1ca}q\x1dh\x08}q\x1esbX\x07\x00\x00\x00threadsq\x1fK\x08X\t\x00\x00\x00resourcesq csnakemake.io\nResources\nq!)\x81q"(K\x08K\x01e}q#(X\x06\x00\x00\x00_coresq$K\x08h\x08}q%(h$K\x00N\x86q&X\x06\x00\x00\x00_nodesq\'K\x01N\x86q(uh\'K\x01ubX\x05\x00\x00\x00inputq)csnakemake.io\nInputFiles\nq*)\x81q+X\x17\x00\x00\x00bams_srr/SRR5227288.bamq,a}q-h\x08}q.sbX\x04\x00\x00\x00ruleq/X\n\x00\x00\x00merge_bamsq0X\x03\x00\x00\x00logq1csnakemake.io\nLog\nq2)\x81q3}q4h\x08}q5sbub.'
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
