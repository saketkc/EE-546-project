######## Snakemake header ########
import sys
sys.path.append(
    "/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"
)
import pickle
snakemake = pickle.loads(
    b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00outputq\x03csnakemake.io\nOutputFiles\nq\x04)\x81q\x05X:\x00\x00\x00plots/metagene_lengthwise/SRX365487/25/3prime_combined.pngq\x06a}q\x07X\x06\x00\x00\x00_namesq\x08}q\tsbX\x05\x00\x00\x00inputq\ncsnakemake.io\nInputFiles\nq\x0b)\x81q\x0cX=\x00\x00\x00metagene_coverage_lengthwise/SRX365487/25/3prime_combined.tsvq\ra}q\x0eh\x08}q\x0fsbX\x06\x00\x00\x00paramsq\x10csnakemake.io\nParams\nq\x11)\x81q\x12}q\x13h\x08}q\x14sbX\x04\x00\x00\x00ruleq\x15X\x1f\x00\x00\x00plot_metagene_individual_lengthq\x16X\t\x00\x00\x00resourcesq\x17csnakemake.io\nResources\nq\x18)\x81q\x19(K\x01K\x01e}q\x1a(h\x08}q\x1b(X\x06\x00\x00\x00_coresq\x1cK\x00N\x86q\x1dX\x06\x00\x00\x00_nodesq\x1eK\x01N\x86q\x1fuh\x1eK\x01h\x1cK\x01ubX\x03\x00\x00\x00logq csnakemake.io\nLog\nq!)\x81q"}q#h\x08}q$sbX\x07\x00\x00\x00threadsq%K\x01X\x06\x00\x00\x00configq&}q\'X\x0b\x00\x00\x00config_pathq(X\x19\x00\x00\x00configs/mm10_SRP031501.pyq)sX\t\x00\x00\x00wildcardsq*csnakemake.io\nWildcards\nq+)\x81q,(X\t\x00\x00\x00SRX365487q-X\x02\x00\x00\x0025q.X\x08\x00\x00\x00combinedq/X\x06\x00\x00\x003primeq0e}q1(X\x06\x00\x00\x00sampleq2h-X\x0f\x00\x00\x00fragment_lengthq3h.X\x0b\x00\x00\x00orientationq4h0X\x06\x00\x00\x00strandq5h/h\x08}q6(X\x06\x00\x00\x00sampleq7K\x00N\x86q8X\x0f\x00\x00\x00fragment_lengthq9K\x01N\x86q:X\x06\x00\x00\x00strandq;K\x02N\x86q<X\x0b\x00\x00\x00orientationq=K\x03N\x86q>uubub.'
)
from snakemake.logging import logger
logger.printshellcmds = True
######## Original script #########
import os
from snakemake.shell import shell
if snakemake.wildcards.orientation == '5prime':
    RANGE = '-60:100'
else:
    RANGE = '-100:60'
if os.stat(str(snakemake.input)).st_size:
    shell(r'''riboraptor plot-metagene \
            --counts {snakemake.input} \
            --saveto {snakemake.output} \
            --positions {RANGE}''')
else:
    # Just touch the file
    shell(r'''touch {snakemake.output}''')
