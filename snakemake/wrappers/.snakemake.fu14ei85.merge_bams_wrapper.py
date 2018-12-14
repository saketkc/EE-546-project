
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00outputq\x03csnakemake.io\nOutputFiles\nq\x04)\x81q\x05X\x12\x00\x00\x00bams/SRX399801.bamq\x06a}q\x07X\x06\x00\x00\x00_namesq\x08}q\tsbX\x03\x00\x00\x00logq\ncsnakemake.io\nLog\nq\x0b)\x81q\x0c}q\rh\x08}q\x0esbX\t\x00\x00\x00resourcesq\x0fcsnakemake.io\nResources\nq\x10)\x81q\x11(K\x01K\x01e}q\x12(X\x06\x00\x00\x00_coresq\x13K\x01X\x06\x00\x00\x00_nodesq\x14K\x01h\x08}q\x15(h\x13K\x01N\x86q\x16h\x14K\x00N\x86q\x17uubX\x06\x00\x00\x00configq\x18}q\x19X\x0b\x00\x00\x00config_pathq\x1aX\x1b\x00\x00\x00configs/GRCz10_SRP034750.pyq\x1bsX\x07\x00\x00\x00threadsq\x1cK\x01X\x04\x00\x00\x00ruleq\x1dX\n\x00\x00\x00merge_bamsq\x1eX\x05\x00\x00\x00inputq\x1fcsnakemake.io\nInputFiles\nq )\x81q!(X\x17\x00\x00\x00bams_srr/SRR1062216.bamq"X\x17\x00\x00\x00bams_srr/SRR1062217.bamq#X\x17\x00\x00\x00bams_srr/SRR1062218.bamq$X\x17\x00\x00\x00bams_srr/SRR1062219.bamq%X\x17\x00\x00\x00bams_srr/SRR1062220.bamq&X\x17\x00\x00\x00bams_srr/SRR1062221.bamq\'X\x17\x00\x00\x00bams_srr/SRR1062222.bamq(X\x17\x00\x00\x00bams_srr/SRR1062223.bamq)X\x17\x00\x00\x00bams_srr/SRR1062224.bamq*X\x17\x00\x00\x00bams_srr/SRR1062225.bamq+e}q,h\x08}q-sbX\t\x00\x00\x00wildcardsq.csnakemake.io\nWildcards\nq/)\x81q0X\t\x00\x00\x00SRX399801q1a}q2(h\x08}q3X\x06\x00\x00\x00sampleq4K\x00N\x86q5sX\x06\x00\x00\x00sampleq6h1ubX\x06\x00\x00\x00paramsq7csnakemake.io\nParams\nq8)\x81q9X\x04\x00\x00\x00/tmpq:a}q;(h\x08}q<X\x07\x00\x00\x00tmp_dirq=K\x00N\x86q>sh=h:ubub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
