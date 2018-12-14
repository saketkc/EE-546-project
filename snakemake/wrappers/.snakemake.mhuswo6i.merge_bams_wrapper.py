######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\t\x00\x00\x00resourcesq\x03csnakemake.io\nResources\nq\x04)\x81q\x05(K\x01K\x10e}q\x06(X\x06\x00\x00\x00_nodesq\x07K\x01X\x06\x00\x00\x00_namesq\x08}q\t(h\x07K\x00N\x86q\nX\x06\x00\x00\x00_coresq\x0bK\x01N\x86q\x0cuh\x0bK\x10ubX\t\x00\x00\x00wildcardsq\rcsnakemake.io\nWildcards\nq\x0e)\x81q\x0fX\t\x00\x00\x00SRX668706q\x10a}q\x11(X\x06\x00\x00\x00sampleq\x12h\x10h\x08}q\x13X\x06\x00\x00\x00sampleq\x14K\x00N\x86q\x15subX\x04\x00\x00\x00ruleq\x16X\n\x00\x00\x00merge_bamsq\x17X\x05\x00\x00\x00inputq\x18csnakemake.io\nInputFiles\nq\x19)\x81q\x1a(X\x17\x00\x00\x00bams_srr/SRR1535528.bamq\x1bX\x17\x00\x00\x00bams_srr/SRR1535529.bamq\x1cX\x17\x00\x00\x00bams_srr/SRR1535530.bamq\x1de}q\x1eh\x08}q\x1fsbX\x06\x00\x00\x00paramsq csnakemake.io\nParams\nq!)\x81q"X\x04\x00\x00\x00/tmpq#a}q$(X\x07\x00\x00\x00tmp_dirq%h#h\x08}q&h%K\x00N\x86q\'subX\x06\x00\x00\x00outputq(csnakemake.io\nOutputFiles\nq))\x81q*X\x12\x00\x00\x00bams/SRX668706.bamq+a}q,h\x08}q-sbX\x03\x00\x00\x00logq.csnakemake.io\nLog\nq/)\x81q0}q1h\x08}q2sbX\x07\x00\x00\x00threadsq3K\x10X\x06\x00\x00\x00configq4}q5X\x0b\x00\x00\x00config_pathq6X\x19\x00\x00\x00configs/hg38_SRP045214.pyq7sub.'
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
