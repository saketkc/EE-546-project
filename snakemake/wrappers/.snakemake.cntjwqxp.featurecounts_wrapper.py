
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00configq\x03}q\x04X\x0b\x00\x00\x00config_pathq\x05X\x19\x00\x00\x00configs/hg38_SRP055009.pyq\x06sX\x05\x00\x00\x00inputq\x07csnakemake.io\nInputFiles\nq\x08)\x81q\t(X\x12\x00\x00\x00bams/SRX876766.bamq\nX\x12\x00\x00\x00bams/SRX876750.bamq\x0bX\x12\x00\x00\x00bams/SRX876732.bamq\x0cX\x12\x00\x00\x00bams/SRX876756.bamq\rX\x12\x00\x00\x00bams/SRX876725.bamq\x0eX\x12\x00\x00\x00bams/SRX876728.bamq\x0fX\x12\x00\x00\x00bams/SRX876737.bamq\x10X\x12\x00\x00\x00bams/SRX876733.bamq\x11X\x12\x00\x00\x00bams/SRX876771.bamq\x12X\x12\x00\x00\x00bams/SRX876764.bamq\x13X\x12\x00\x00\x00bams/SRX876762.bamq\x14X\x12\x00\x00\x00bams/SRX876731.bamq\x15X\x12\x00\x00\x00bams/SRX876761.bamq\x16X\x12\x00\x00\x00bams/SRX876772.bamq\x17X\x12\x00\x00\x00bams/SRX876743.bamq\x18X\x12\x00\x00\x00bams/SRX876760.bamq\x19X\x12\x00\x00\x00bams/SRX876735.bamq\x1aX\x12\x00\x00\x00bams/SRX876758.bamq\x1bX\x12\x00\x00\x00bams/SRX876755.bamq\x1cX\x12\x00\x00\x00bams/SRX876768.bamq\x1dX\x12\x00\x00\x00bams/SRX876744.bamq\x1eX\x12\x00\x00\x00bams/SRX876763.bamq\x1fX\x12\x00\x00\x00bams/SRX876734.bamq X\x12\x00\x00\x00bams/SRX876770.bamq!X\x12\x00\x00\x00bams/SRX876726.bamq"X\x12\x00\x00\x00bams/SRX876751.bamq#X\x12\x00\x00\x00bams/SRX876765.bamq$X\x12\x00\x00\x00bams/SRX876769.bamq%X\x12\x00\x00\x00bams/SRX876740.bamq&X\x12\x00\x00\x00bams/SRX876747.bamq\'X\x12\x00\x00\x00bams/SRX876727.bamq(X\x12\x00\x00\x00bams/SRX876724.bamq)X\x12\x00\x00\x00bams/SRX876738.bamq*X\x12\x00\x00\x00bams/SRX876741.bamq+X\x12\x00\x00\x00bams/SRX876742.bamq,X\x12\x00\x00\x00bams/SRX876759.bamq-X\x12\x00\x00\x00bams/SRX876736.bamq.X\x12\x00\x00\x00bams/SRX876746.bamq/X\x12\x00\x00\x00bams/SRX876752.bamq0X\x12\x00\x00\x00bams/SRX876745.bamq1X\x12\x00\x00\x00bams/SRX876723.bamq2X\x12\x00\x00\x00bams/SRX876754.bamq3X\x12\x00\x00\x00bams/SRX876749.bamq4X\x12\x00\x00\x00bams/SRX876730.bamq5X\x12\x00\x00\x00bams/SRX876757.bamq6X\x12\x00\x00\x00bams/SRX876767.bamq7X\x12\x00\x00\x00bams/SRX876748.bamq8X\x12\x00\x00\x00hdf/SRX876766.hdf5q9X\x12\x00\x00\x00hdf/SRX876750.hdf5q:X\x12\x00\x00\x00hdf/SRX876732.hdf5q;X\x12\x00\x00\x00hdf/SRX876756.hdf5q<X\x12\x00\x00\x00hdf/SRX876725.hdf5q=X\x12\x00\x00\x00hdf/SRX876728.hdf5q>X\x12\x00\x00\x00hdf/SRX876737.hdf5q?X\x12\x00\x00\x00hdf/SRX876733.hdf5q@X\x12\x00\x00\x00hdf/SRX876771.hdf5qAX\x12\x00\x00\x00hdf/SRX876764.hdf5qBX\x12\x00\x00\x00hdf/SRX876762.hdf5qCX\x12\x00\x00\x00hdf/SRX876731.hdf5qDX\x12\x00\x00\x00hdf/SRX876761.hdf5qEX\x12\x00\x00\x00hdf/SRX876772.hdf5qFX\x12\x00\x00\x00hdf/SRX876743.hdf5qGX\x12\x00\x00\x00hdf/SRX876760.hdf5qHX\x12\x00\x00\x00hdf/SRX876735.hdf5qIX\x12\x00\x00\x00hdf/SRX876758.hdf5qJX\x12\x00\x00\x00hdf/SRX876755.hdf5qKX\x12\x00\x00\x00hdf/SRX876768.hdf5qLX\x12\x00\x00\x00hdf/SRX876744.hdf5qMX\x12\x00\x00\x00hdf/SRX876763.hdf5qNX\x12\x00\x00\x00hdf/SRX876734.hdf5qOX\x12\x00\x00\x00hdf/SRX876770.hdf5qPX\x12\x00\x00\x00hdf/SRX876726.hdf5qQX\x12\x00\x00\x00hdf/SRX876751.hdf5qRX\x12\x00\x00\x00hdf/SRX876765.hdf5qSX\x12\x00\x00\x00hdf/SRX876769.hdf5qTX\x12\x00\x00\x00hdf/SRX876740.hdf5qUX\x12\x00\x00\x00hdf/SRX876747.hdf5qVX\x12\x00\x00\x00hdf/SRX876727.hdf5qWX\x12\x00\x00\x00hdf/SRX876724.hdf5qXX\x12\x00\x00\x00hdf/SRX876738.hdf5qYX\x12\x00\x00\x00hdf/SRX876741.hdf5qZX\x12\x00\x00\x00hdf/SRX876742.hdf5q[X\x12\x00\x00\x00hdf/SRX876759.hdf5q\\X\x12\x00\x00\x00hdf/SRX876736.hdf5q]X\x12\x00\x00\x00hdf/SRX876746.hdf5q^X\x12\x00\x00\x00hdf/SRX876752.hdf5q_X\x12\x00\x00\x00hdf/SRX876745.hdf5q`X\x12\x00\x00\x00hdf/SRX876723.hdf5qaX\x12\x00\x00\x00hdf/SRX876754.hdf5qbX\x12\x00\x00\x00hdf/SRX876749.hdf5qcX\x12\x00\x00\x00hdf/SRX876730.hdf5qdX\x12\x00\x00\x00hdf/SRX876757.hdf5qeX\x12\x00\x00\x00hdf/SRX876767.hdf5qfX\x12\x00\x00\x00hdf/SRX876748.hdf5qge}qh(X\x04\x00\x00\x00bamsqicsnakemake.io\nNamedlist\nqj)\x81qk(h\nh\x0bh\x0ch\rh\x0eh\x0fh\x10h\x11h\x12h\x13h\x14h\x15h\x16h\x17h\x18h\x19h\x1ah\x1bh\x1ch\x1dh\x1eh\x1fh h!h"h#h$h%h&h\'h(h)h*h+h,h-h.h/h0h1h2h3h4h5h6h7h8e}qlX\x06\x00\x00\x00_namesqm}qnsbhm}qo(hiK\x00K/\x86qpX\x04\x00\x00\x00hdfsqqK/K^\x86qruhqhj)\x81qs(h9h:h;h<h=h>h?h@hAhBhChDhEhFhGhHhIhJhKhLhMhNhOhPhQhRhShThUhVhWhXhYhZh[h\\h]h^h_h`hahbhchdhehfhge}qthm}qusbubX\x03\x00\x00\x00logqvcsnakemake.io\nLog\nqw)\x81qx}qyhm}qzsbX\x06\x00\x00\x00outputq{csnakemake.io\nOutputFiles\nq|)\x81q}X\x19\x00\x00\x00featureCounts/fcounts.tsvq~a}q\x7fhm}q\x80sbX\x04\x00\x00\x00ruleq\x81X\r\x00\x00\x00featurecountsq\x82X\t\x00\x00\x00wildcardsq\x83csnakemake.io\nWildcards\nq\x84)\x81q\x85}q\x86hm}q\x87sbX\x06\x00\x00\x00paramsq\x88csnakemake.io\nParams\nq\x89)\x81q\x8aXN\x00\x00\x00/home/cmb-panasas2/skchoudh/genomes/hg38/annotation/gencode.v25.annotation.gtfq\x8ba}q\x8c(X\n\x00\x00\x00annotationq\x8dh\x8bhm}q\x8eh\x8dK\x00N\x86q\x8fsubX\x07\x00\x00\x00threadsq\x90K\x01X\t\x00\x00\x00resourcesq\x91csnakemake.io\nResources\nq\x92)\x81q\x93(K\x01K\x01e}q\x94(X\x06\x00\x00\x00_coresq\x95K\x01hm}q\x96(h\x95K\x01N\x86q\x97X\x06\x00\x00\x00_nodesq\x98K\x00N\x86q\x99uh\x98K\x01ubub.'); from snakemake.logging import logger; logger.printshellcmds = True
######## Original script #########
import os
import h5py
from snakemake.shell import shell
from collections import defaultdict
import pandas as pd

protocols = defaultdict(list)
for hdf, bam in zip(snakemake.input['hdfs'], snakemake.input['bams']):
    hdf = h5py.File(hdf, 'r')
    protocol = hdf.attrs['protocol']
    protocols[protocol].append(bam)
    hdf.close()
outputs = []
for protocol, bams in protocols.items():
    if protocol == 'forward':
        count_strat = '-s 1'
    elif protocol == 'unstranded':
        count_strat = '-s 2'
    else:
        count_strat = ''
    bams = sorted(bams)
    output = os.path.abspath(str(snakemake.output)) + '-' + protocol
    outputs.append(output)
    shell(
        r'''featureCounts {count_strat} -a {snakemake.params.annotation} -o {output} -t exon -g gene_id -Q 4 -T {snakemake.threads} {bams}'''
    )
df = pd.read_table(outputs[0], skiprows=[0])
if len(outputs) > 1:
    for f in outputs[1:]:
        temp_df = pd.read_table(f, skiprows=[0])
        temp_df = temp_df.drop(
            columns=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'])
        df = pd.concat([df, temp_df], axis=1)

df.to_csv(str(snakemake.output), sep='\t', index=False, header=True)
