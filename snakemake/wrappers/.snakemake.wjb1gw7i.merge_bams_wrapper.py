
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(X\x17\x00\x00\x00bams_srr/SRR1062197.bamq\x06X\x17\x00\x00\x00bams_srr/SRR1062198.bamq\x07X\x17\x00\x00\x00bams_srr/SRR1062199.bamq\x08X\x17\x00\x00\x00bams_srr/SRR1062200.bamq\tX\x17\x00\x00\x00bams_srr/SRR1062201.bamq\nX\x17\x00\x00\x00bams_srr/SRR1062202.bamq\x0bX\x17\x00\x00\x00bams_srr/SRR1062203.bamq\x0cX\x17\x00\x00\x00bams_srr/SRR1062204.bamq\rX\x17\x00\x00\x00bams_srr/SRR1062205.bamq\x0eX\x17\x00\x00\x00bams_srr/SRR1062206.bamq\x0fe}q\x10X\x06\x00\x00\x00_namesq\x11}q\x12sbX\x04\x00\x00\x00ruleq\x13X\n\x00\x00\x00merge_bamsq\x14X\x03\x00\x00\x00logq\x15csnakemake.io\nLog\nq\x16)\x81q\x17}q\x18h\x11}q\x19sbX\t\x00\x00\x00resourcesq\x1acsnakemake.io\nResources\nq\x1b)\x81q\x1c(K\x01K\x01e}q\x1d(h\x11}q\x1e(X\x06\x00\x00\x00_coresq\x1fK\x00N\x86q X\x06\x00\x00\x00_nodesq!K\x01N\x86q"uh!K\x01h\x1fK\x01ubX\t\x00\x00\x00wildcardsq#csnakemake.io\nWildcards\nq$)\x81q%X\t\x00\x00\x00SRX399799q&a}q\'(h\x11}q(X\x06\x00\x00\x00sampleq)K\x00N\x86q*sX\x06\x00\x00\x00sampleq+h&ubX\x06\x00\x00\x00outputq,csnakemake.io\nOutputFiles\nq-)\x81q.X\x12\x00\x00\x00bams/SRX399799.bamq/a}q0h\x11}q1sbX\x07\x00\x00\x00threadsq2K\x01X\x06\x00\x00\x00configq3}q4X\x0b\x00\x00\x00config_pathq5X\x1b\x00\x00\x00configs/GRCz10_SRP034750.pyq6sX\x06\x00\x00\x00paramsq7csnakemake.io\nParams\nq8)\x81q9X\x04\x00\x00\x00/tmpq:a}q;(X\x07\x00\x00\x00tmp_dirq<h:h\x11}q=h<K\x00N\x86q>subub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
