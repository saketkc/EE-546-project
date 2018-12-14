
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\t\x00\x00\x00wildcardsq\x03csnakemake.io\nWildcards\nq\x04)\x81q\x05X\t\x00\x00\x00SRX399811q\x06a}q\x07(X\x06\x00\x00\x00_namesq\x08}q\tX\x06\x00\x00\x00sampleq\nK\x00N\x86q\x0bsX\x06\x00\x00\x00sampleq\x0ch\x06ubX\x05\x00\x00\x00inputq\rcsnakemake.io\nInputFiles\nq\x0e)\x81q\x0f(X\x17\x00\x00\x00bams_srr/SRR1062285.bamq\x10X\x17\x00\x00\x00bams_srr/SRR1062286.bamq\x11X\x17\x00\x00\x00bams_srr/SRR1062287.bamq\x12X\x17\x00\x00\x00bams_srr/SRR1062288.bamq\x13X\x17\x00\x00\x00bams_srr/SRR1062289.bamq\x14X\x17\x00\x00\x00bams_srr/SRR1062290.bamq\x15X\x17\x00\x00\x00bams_srr/SRR1062291.bamq\x16X\x17\x00\x00\x00bams_srr/SRR1062292.bamq\x17X\x17\x00\x00\x00bams_srr/SRR1062293.bamq\x18e}q\x19h\x08}q\x1asbX\t\x00\x00\x00resourcesq\x1bcsnakemake.io\nResources\nq\x1c)\x81q\x1d(K\x01K\x01e}q\x1e(h\x08}q\x1f(X\x06\x00\x00\x00_nodesq K\x00N\x86q!X\x06\x00\x00\x00_coresq"K\x01N\x86q#uh"K\x01h K\x01ubX\x03\x00\x00\x00logq$csnakemake.io\nLog\nq%)\x81q&}q\'h\x08}q(sbX\x04\x00\x00\x00ruleq)X\n\x00\x00\x00merge_bamsq*X\x07\x00\x00\x00threadsq+K\x01X\x06\x00\x00\x00configq,}q-X\x0b\x00\x00\x00config_pathq.X\x1b\x00\x00\x00configs/GRCz10_SRP034750.pyq/sX\x06\x00\x00\x00paramsq0csnakemake.io\nParams\nq1)\x81q2X\x04\x00\x00\x00/tmpq3a}q4(h\x08}q5X\x07\x00\x00\x00tmp_dirq6K\x00N\x86q7sh6h3ubX\x06\x00\x00\x00outputq8csnakemake.io\nOutputFiles\nq9)\x81q:X\x12\x00\x00\x00bams/SRX399811.bamq;a}q<h\x08}q=sbub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
