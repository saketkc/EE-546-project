
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x04\x00\x00\x00ruleq\x03X\n\x00\x00\x00merge_bamsq\x04X\x07\x00\x00\x00threadsq\x05K\x01X\x03\x00\x00\x00logq\x06csnakemake.io\nLog\nq\x07)\x81q\x08}q\tX\x06\x00\x00\x00_namesq\n}q\x0bsbX\t\x00\x00\x00wildcardsq\x0ccsnakemake.io\nWildcards\nq\r)\x81q\x0eX\t\x00\x00\x00SRX399805q\x0fa}q\x10(X\x06\x00\x00\x00sampleq\x11h\x0fh\n}q\x12X\x06\x00\x00\x00sampleq\x13K\x00N\x86q\x14subX\t\x00\x00\x00resourcesq\x15csnakemake.io\nResources\nq\x16)\x81q\x17(K\x01K\x01e}q\x18(h\n}q\x19(X\x06\x00\x00\x00_coresq\x1aK\x00N\x86q\x1bX\x06\x00\x00\x00_nodesq\x1cK\x01N\x86q\x1duh\x1aK\x01h\x1cK\x01ubX\x06\x00\x00\x00configq\x1e}q\x1fX\x0b\x00\x00\x00config_pathq X\x1b\x00\x00\x00configs/GRCz10_SRP034750.pyq!sX\x05\x00\x00\x00inputq"csnakemake.io\nInputFiles\nq#)\x81q$(X\x17\x00\x00\x00bams_srr/SRR1062257.bamq%X\x17\x00\x00\x00bams_srr/SRR1062258.bamq&X\x17\x00\x00\x00bams_srr/SRR1062259.bamq\'X\x17\x00\x00\x00bams_srr/SRR1062260.bamq(X\x17\x00\x00\x00bams_srr/SRR1062261.bamq)X\x17\x00\x00\x00bams_srr/SRR1062262.bamq*X\x17\x00\x00\x00bams_srr/SRR1062263.bamq+X\x17\x00\x00\x00bams_srr/SRR1062264.bamq,X\x17\x00\x00\x00bams_srr/SRR1062265.bamq-X\x17\x00\x00\x00bams_srr/SRR1062266.bamq.X\x17\x00\x00\x00bams_srr/SRR1062267.bamq/X\x17\x00\x00\x00bams_srr/SRR1062268.bamq0X\x17\x00\x00\x00bams_srr/SRR1062269.bamq1X\x17\x00\x00\x00bams_srr/SRR1062270.bamq2e}q3h\n}q4sbX\x06\x00\x00\x00outputq5csnakemake.io\nOutputFiles\nq6)\x81q7X\x12\x00\x00\x00bams/SRX399805.bamq8a}q9h\n}q:sbX\x06\x00\x00\x00paramsq;csnakemake.io\nParams\nq<)\x81q=X\x04\x00\x00\x00/tmpq>a}q?(X\x07\x00\x00\x00tmp_dirq@h>h\n}qAh@K\x00N\x86qBsubub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
