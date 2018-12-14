
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(X\x17\x00\x00\x00bams_srr/SRR1062455.bamq\x06X\x17\x00\x00\x00bams_srr/SRR1062456.bamq\x07X\x17\x00\x00\x00bams_srr/SRR1062457.bamq\x08X\x17\x00\x00\x00bams_srr/SRR1062458.bamq\tX\x17\x00\x00\x00bams_srr/SRR1062459.bamq\nX\x17\x00\x00\x00bams_srr/SRR1062460.bamq\x0bX\x17\x00\x00\x00bams_srr/SRR1062461.bamq\x0cX\x17\x00\x00\x00bams_srr/SRR1062462.bamq\rX\x17\x00\x00\x00bams_srr/SRR1062463.bamq\x0eX\x17\x00\x00\x00bams_srr/SRR1062464.bamq\x0fX\x17\x00\x00\x00bams_srr/SRR1062465.bamq\x10X\x17\x00\x00\x00bams_srr/SRR1062466.bamq\x11X\x17\x00\x00\x00bams_srr/SRR1062467.bamq\x12X\x17\x00\x00\x00bams_srr/SRR1062468.bamq\x13X\x17\x00\x00\x00bams_srr/SRR1062469.bamq\x14X\x17\x00\x00\x00bams_srr/SRR1062470.bamq\x15X\x17\x00\x00\x00bams_srr/SRR1062471.bamq\x16X\x17\x00\x00\x00bams_srr/SRR1062472.bamq\x17X\x17\x00\x00\x00bams_srr/SRR1062473.bamq\x18X\x17\x00\x00\x00bams_srr/SRR1062474.bamq\x19X\x17\x00\x00\x00bams_srr/SRR1062475.bamq\x1aX\x17\x00\x00\x00bams_srr/SRR1062476.bamq\x1bX\x17\x00\x00\x00bams_srr/SRR1062477.bamq\x1cX\x17\x00\x00\x00bams_srr/SRR1062478.bamq\x1dX\x17\x00\x00\x00bams_srr/SRR1062479.bamq\x1eX\x17\x00\x00\x00bams_srr/SRR1062480.bamq\x1fX\x17\x00\x00\x00bams_srr/SRR1062481.bamq X\x17\x00\x00\x00bams_srr/SRR1062482.bamq!X\x17\x00\x00\x00bams_srr/SRR1062483.bamq"X\x17\x00\x00\x00bams_srr/SRR1062484.bamq#X\x17\x00\x00\x00bams_srr/SRR1062485.bamq$X\x17\x00\x00\x00bams_srr/SRR1062486.bamq%X\x17\x00\x00\x00bams_srr/SRR1062487.bamq&X\x17\x00\x00\x00bams_srr/SRR1062488.bamq\'X\x17\x00\x00\x00bams_srr/SRR1062489.bamq(X\x17\x00\x00\x00bams_srr/SRR1062490.bamq)X\x17\x00\x00\x00bams_srr/SRR1062491.bamq*X\x17\x00\x00\x00bams_srr/SRR1062492.bamq+X\x17\x00\x00\x00bams_srr/SRR1062493.bamq,X\x17\x00\x00\x00bams_srr/SRR1062494.bamq-X\x17\x00\x00\x00bams_srr/SRR1062495.bamq.X\x17\x00\x00\x00bams_srr/SRR1062496.bamq/X\x17\x00\x00\x00bams_srr/SRR1062497.bamq0X\x17\x00\x00\x00bams_srr/SRR1062498.bamq1X\x17\x00\x00\x00bams_srr/SRR1062499.bamq2X\x17\x00\x00\x00bams_srr/SRR1062500.bamq3X\x17\x00\x00\x00bams_srr/SRR1062501.bamq4X\x17\x00\x00\x00bams_srr/SRR1062502.bamq5X\x17\x00\x00\x00bams_srr/SRR1062503.bamq6X\x17\x00\x00\x00bams_srr/SRR1062504.bamq7X\x17\x00\x00\x00bams_srr/SRR1062505.bamq8X\x17\x00\x00\x00bams_srr/SRR1062506.bamq9X\x17\x00\x00\x00bams_srr/SRR1062507.bamq:X\x17\x00\x00\x00bams_srr/SRR1062508.bamq;X\x17\x00\x00\x00bams_srr/SRR1062509.bamq<X\x17\x00\x00\x00bams_srr/SRR1062510.bamq=X\x17\x00\x00\x00bams_srr/SRR1062511.bamq>X\x17\x00\x00\x00bams_srr/SRR1062512.bamq?X\x17\x00\x00\x00bams_srr/SRR1062513.bamq@X\x17\x00\x00\x00bams_srr/SRR1062514.bamqAX\x17\x00\x00\x00bams_srr/SRR1062515.bamqBX\x17\x00\x00\x00bams_srr/SRR1062516.bamqCX\x17\x00\x00\x00bams_srr/SRR1062517.bamqDX\x17\x00\x00\x00bams_srr/SRR1062518.bamqEX\x17\x00\x00\x00bams_srr/SRR1062519.bamqFX\x17\x00\x00\x00bams_srr/SRR1062520.bamqGX\x17\x00\x00\x00bams_srr/SRR1062521.bamqHX\x17\x00\x00\x00bams_srr/SRR1062522.bamqIX\x17\x00\x00\x00bams_srr/SRR1062523.bamqJX\x17\x00\x00\x00bams_srr/SRR1062524.bamqKX\x17\x00\x00\x00bams_srr/SRR1062525.bamqLX\x17\x00\x00\x00bams_srr/SRR1062526.bamqMX\x17\x00\x00\x00bams_srr/SRR1062527.bamqNX\x17\x00\x00\x00bams_srr/SRR1062528.bamqOX\x17\x00\x00\x00bams_srr/SRR1062529.bamqPX\x17\x00\x00\x00bams_srr/SRR1062530.bamqQX\x17\x00\x00\x00bams_srr/SRR1062531.bamqRX\x17\x00\x00\x00bams_srr/SRR1062532.bamqSX\x17\x00\x00\x00bams_srr/SRR1062533.bamqTX\x17\x00\x00\x00bams_srr/SRR1062534.bamqUX\x17\x00\x00\x00bams_srr/SRR1062535.bamqVX\x17\x00\x00\x00bams_srr/SRR1062536.bamqWX\x17\x00\x00\x00bams_srr/SRR1062537.bamqXX\x17\x00\x00\x00bams_srr/SRR1062538.bamqYX\x17\x00\x00\x00bams_srr/SRR1062539.bamqZX\x17\x00\x00\x00bams_srr/SRR1062540.bamq[X\x17\x00\x00\x00bams_srr/SRR1062541.bamq\\X\x17\x00\x00\x00bams_srr/SRR1062542.bamq]X\x17\x00\x00\x00bams_srr/SRR1062543.bamq^X\x17\x00\x00\x00bams_srr/SRR1062544.bamq_X\x17\x00\x00\x00bams_srr/SRR1062545.bamq`X\x17\x00\x00\x00bams_srr/SRR1062546.bamqaX\x17\x00\x00\x00bams_srr/SRR1062547.bamqbX\x17\x00\x00\x00bams_srr/SRR1062548.bamqcX\x17\x00\x00\x00bams_srr/SRR1062549.bamqdX\x17\x00\x00\x00bams_srr/SRR1062550.bamqeX\x17\x00\x00\x00bams_srr/SRR1062551.bamqfX\x17\x00\x00\x00bams_srr/SRR1062552.bamqgX\x17\x00\x00\x00bams_srr/SRR1062553.bamqhX\x17\x00\x00\x00bams_srr/SRR1062554.bamqie}qjX\x06\x00\x00\x00_namesqk}qlsbX\t\x00\x00\x00wildcardsqmcsnakemake.io\nWildcards\nqn)\x81qoX\t\x00\x00\x00SRX399822qpa}qq(hk}qrX\x06\x00\x00\x00sampleqsK\x00N\x86qtsX\x06\x00\x00\x00samplequhpubX\x06\x00\x00\x00outputqvcsnakemake.io\nOutputFiles\nqw)\x81qxX\x12\x00\x00\x00bams/SRX399822.bamqya}qzhk}q{sbX\x03\x00\x00\x00logq|csnakemake.io\nLog\nq})\x81q~}q\x7fhk}q\x80sbX\x07\x00\x00\x00threadsq\x81K\x01X\t\x00\x00\x00resourcesq\x82csnakemake.io\nResources\nq\x83)\x81q\x84(K\x01K\x01e}q\x85(X\x06\x00\x00\x00_nodesq\x86K\x01X\x06\x00\x00\x00_coresq\x87K\x01hk}q\x88(h\x86K\x00N\x86q\x89h\x87K\x01N\x86q\x8auubX\x06\x00\x00\x00paramsq\x8bcsnakemake.io\nParams\nq\x8c)\x81q\x8dX\x04\x00\x00\x00/tmpq\x8ea}q\x8f(hk}q\x90X\x07\x00\x00\x00tmp_dirq\x91K\x00N\x86q\x92sh\x91h\x8eubX\x04\x00\x00\x00ruleq\x93X\n\x00\x00\x00merge_bamsq\x94X\x06\x00\x00\x00configq\x95}q\x96X\x0b\x00\x00\x00config_pathq\x97X\x1b\x00\x00\x00configs/GRCz10_SRP034750.pyq\x98sub.'); from snakemake.logging import logger; logger.printshellcmds = True
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