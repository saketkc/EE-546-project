
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x07\x00\x00\x00threadsq\x03K\x01X\t\x00\x00\x00resourcesq\x04csnakemake.io\nResources\nq\x05)\x81q\x06(K\x01K\x01e}q\x07(X\x06\x00\x00\x00_namesq\x08}q\t(X\x06\x00\x00\x00_nodesq\nK\x00N\x86q\x0bX\x06\x00\x00\x00_coresq\x0cK\x01N\x86q\ruh\nK\x01h\x0cK\x01ubX\x03\x00\x00\x00logq\x0ecsnakemake.io\nLog\nq\x0f)\x81q\x10}q\x11h\x08}q\x12sbX\x04\x00\x00\x00ruleq\x13X\n\x00\x00\x00merge_bamsq\x14X\x06\x00\x00\x00paramsq\x15csnakemake.io\nParams\nq\x16)\x81q\x17X\x04\x00\x00\x00/tmpq\x18a}q\x19(X\x07\x00\x00\x00tmp_dirq\x1ah\x18h\x08}q\x1bh\x1aK\x00N\x86q\x1csubX\x05\x00\x00\x00inputq\x1dcsnakemake.io\nInputFiles\nq\x1e)\x81q\x1f(X\x17\x00\x00\x00bams_srr/SRR1062767.bamq X\x17\x00\x00\x00bams_srr/SRR1062768.bamq!X\x17\x00\x00\x00bams_srr/SRR1062769.bamq"X\x17\x00\x00\x00bams_srr/SRR1062770.bamq#X\x17\x00\x00\x00bams_srr/SRR1062771.bamq$X\x17\x00\x00\x00bams_srr/SRR1062772.bamq%X\x17\x00\x00\x00bams_srr/SRR1062773.bamq&X\x17\x00\x00\x00bams_srr/SRR1062774.bamq\'X\x17\x00\x00\x00bams_srr/SRR1062775.bamq(X\x17\x00\x00\x00bams_srr/SRR1062776.bamq)X\x17\x00\x00\x00bams_srr/SRR1062777.bamq*X\x17\x00\x00\x00bams_srr/SRR1062778.bamq+X\x17\x00\x00\x00bams_srr/SRR1062779.bamq,X\x17\x00\x00\x00bams_srr/SRR1062780.bamq-X\x17\x00\x00\x00bams_srr/SRR1062781.bamq.X\x17\x00\x00\x00bams_srr/SRR1062782.bamq/X\x17\x00\x00\x00bams_srr/SRR1062783.bamq0X\x17\x00\x00\x00bams_srr/SRR1062784.bamq1X\x17\x00\x00\x00bams_srr/SRR1062785.bamq2e}q3h\x08}q4sbX\x06\x00\x00\x00outputq5csnakemake.io\nOutputFiles\nq6)\x81q7X\x12\x00\x00\x00bams/SRX399826.bamq8a}q9h\x08}q:sbX\x06\x00\x00\x00configq;}q<X\x0b\x00\x00\x00config_pathq=X\x1b\x00\x00\x00configs/GRCz10_SRP034750.pyq>sX\t\x00\x00\x00wildcardsq?csnakemake.io\nWildcards\nq@)\x81qAX\t\x00\x00\x00SRX399826qBa}qC(X\x06\x00\x00\x00sampleqDhBh\x08}qEX\x06\x00\x00\x00sampleqFK\x00N\x86qGsubub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
