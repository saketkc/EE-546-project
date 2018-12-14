
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00outputq\x03csnakemake.io\nOutputFiles\nq\x04)\x81q\x05X\x19\x00\x00\x00featureCounts/fcounts.tsvq\x06a}q\x07X\x06\x00\x00\x00_namesq\x08}q\tsbX\x04\x00\x00\x00ruleq\nX\r\x00\x00\x00featurecountsq\x0bX\x06\x00\x00\x00configq\x0c}q\rX\x0b\x00\x00\x00config_pathq\x0eX\x19\x00\x00\x00configs/hg38_SRP055009.pyq\x0fsX\x05\x00\x00\x00inputq\x10csnakemake.io\nInputFiles\nq\x11)\x81q\x12(X\x12\x00\x00\x00hdf/SRX876747.hdf5q\x13X\x12\x00\x00\x00hdf/SRX876723.hdf5q\x14X\x12\x00\x00\x00hdf/SRX876740.hdf5q\x15X\x12\x00\x00\x00hdf/SRX876731.hdf5q\x16X\x12\x00\x00\x00hdf/SRX876738.hdf5q\x17X\x12\x00\x00\x00hdf/SRX876771.hdf5q\x18X\x12\x00\x00\x00hdf/SRX876768.hdf5q\x19X\x12\x00\x00\x00hdf/SRX876746.hdf5q\x1aX\x12\x00\x00\x00hdf/SRX876733.hdf5q\x1bX\x12\x00\x00\x00hdf/SRX876726.hdf5q\x1cX\x12\x00\x00\x00hdf/SRX876732.hdf5q\x1dX\x12\x00\x00\x00hdf/SRX876742.hdf5q\x1eX\x12\x00\x00\x00hdf/SRX876735.hdf5q\x1fX\x12\x00\x00\x00hdf/SRX876748.hdf5q X\x12\x00\x00\x00hdf/SRX876756.hdf5q!X\x12\x00\x00\x00hdf/SRX876767.hdf5q"X\x12\x00\x00\x00hdf/SRX876770.hdf5q#X\x12\x00\x00\x00hdf/SRX876749.hdf5q$X\x12\x00\x00\x00hdf/SRX876762.hdf5q%X\x12\x00\x00\x00hdf/SRX876744.hdf5q&X\x12\x00\x00\x00hdf/SRX876764.hdf5q\'X\x12\x00\x00\x00hdf/SRX876751.hdf5q(X\x12\x00\x00\x00hdf/SRX876741.hdf5q)X\x12\x00\x00\x00hdf/SRX876769.hdf5q*X\x12\x00\x00\x00hdf/SRX876758.hdf5q+X\x12\x00\x00\x00hdf/SRX876765.hdf5q,X\x12\x00\x00\x00hdf/SRX876759.hdf5q-X\x12\x00\x00\x00hdf/SRX876763.hdf5q.X\x12\x00\x00\x00hdf/SRX876757.hdf5q/X\x12\x00\x00\x00hdf/SRX876760.hdf5q0X\x12\x00\x00\x00hdf/SRX876737.hdf5q1X\x12\x00\x00\x00hdf/SRX876725.hdf5q2X\x12\x00\x00\x00hdf/SRX876736.hdf5q3X\x12\x00\x00\x00hdf/SRX876761.hdf5q4X\x12\x00\x00\x00hdf/SRX876745.hdf5q5X\x12\x00\x00\x00hdf/SRX876743.hdf5q6X\x12\x00\x00\x00hdf/SRX876734.hdf5q7X\x12\x00\x00\x00hdf/SRX876730.hdf5q8X\x12\x00\x00\x00hdf/SRX876724.hdf5q9X\x12\x00\x00\x00hdf/SRX876755.hdf5q:X\x12\x00\x00\x00hdf/SRX876754.hdf5q;X\x12\x00\x00\x00hdf/SRX876772.hdf5q<X\x12\x00\x00\x00hdf/SRX876752.hdf5q=X\x12\x00\x00\x00hdf/SRX876750.hdf5q>X\x12\x00\x00\x00hdf/SRX876727.hdf5q?X\x12\x00\x00\x00hdf/SRX876766.hdf5q@X\x12\x00\x00\x00hdf/SRX876728.hdf5qAX\x12\x00\x00\x00bams/SRX876747.bamqBX\x12\x00\x00\x00bams/SRX876723.bamqCX\x12\x00\x00\x00bams/SRX876740.bamqDX\x12\x00\x00\x00bams/SRX876731.bamqEX\x12\x00\x00\x00bams/SRX876738.bamqFX\x12\x00\x00\x00bams/SRX876771.bamqGX\x12\x00\x00\x00bams/SRX876768.bamqHX\x12\x00\x00\x00bams/SRX876746.bamqIX\x12\x00\x00\x00bams/SRX876733.bamqJX\x12\x00\x00\x00bams/SRX876726.bamqKX\x12\x00\x00\x00bams/SRX876732.bamqLX\x12\x00\x00\x00bams/SRX876742.bamqMX\x12\x00\x00\x00bams/SRX876735.bamqNX\x12\x00\x00\x00bams/SRX876748.bamqOX\x12\x00\x00\x00bams/SRX876756.bamqPX\x12\x00\x00\x00bams/SRX876767.bamqQX\x12\x00\x00\x00bams/SRX876770.bamqRX\x12\x00\x00\x00bams/SRX876749.bamqSX\x12\x00\x00\x00bams/SRX876762.bamqTX\x12\x00\x00\x00bams/SRX876744.bamqUX\x12\x00\x00\x00bams/SRX876764.bamqVX\x12\x00\x00\x00bams/SRX876751.bamqWX\x12\x00\x00\x00bams/SRX876741.bamqXX\x12\x00\x00\x00bams/SRX876769.bamqYX\x12\x00\x00\x00bams/SRX876758.bamqZX\x12\x00\x00\x00bams/SRX876765.bamq[X\x12\x00\x00\x00bams/SRX876759.bamq\\X\x12\x00\x00\x00bams/SRX876763.bamq]X\x12\x00\x00\x00bams/SRX876757.bamq^X\x12\x00\x00\x00bams/SRX876760.bamq_X\x12\x00\x00\x00bams/SRX876737.bamq`X\x12\x00\x00\x00bams/SRX876725.bamqaX\x12\x00\x00\x00bams/SRX876736.bamqbX\x12\x00\x00\x00bams/SRX876761.bamqcX\x12\x00\x00\x00bams/SRX876745.bamqdX\x12\x00\x00\x00bams/SRX876743.bamqeX\x12\x00\x00\x00bams/SRX876734.bamqfX\x12\x00\x00\x00bams/SRX876730.bamqgX\x12\x00\x00\x00bams/SRX876724.bamqhX\x12\x00\x00\x00bams/SRX876755.bamqiX\x12\x00\x00\x00bams/SRX876754.bamqjX\x12\x00\x00\x00bams/SRX876772.bamqkX\x12\x00\x00\x00bams/SRX876752.bamqlX\x12\x00\x00\x00bams/SRX876750.bamqmX\x12\x00\x00\x00bams/SRX876727.bamqnX\x12\x00\x00\x00bams/SRX876766.bamqoX\x12\x00\x00\x00bams/SRX876728.bamqpe}qq(X\x04\x00\x00\x00hdfsqrcsnakemake.io\nNamedlist\nqs)\x81qt(h\x13h\x14h\x15h\x16h\x17h\x18h\x19h\x1ah\x1bh\x1ch\x1dh\x1eh\x1fh h!h"h#h$h%h&h\'h(h)h*h+h,h-h.h/h0h1h2h3h4h5h6h7h8h9h:h;h<h=h>h?h@hAe}quh\x08}qvsbX\x04\x00\x00\x00bamsqwhs)\x81qx(hBhChDhEhFhGhHhIhJhKhLhMhNhOhPhQhRhShThUhVhWhXhYhZh[h\\h]h^h_h`hahbhchdhehfhghhhihjhkhlhmhnhohpe}qyh\x08}qzsbh\x08}q{(hrK\x00K/\x86q|hwK/K^\x86q}uubX\x07\x00\x00\x00threadsq~K\x01X\x06\x00\x00\x00paramsq\x7fcsnakemake.io\nParams\nq\x80)\x81q\x81XN\x00\x00\x00/home/cmb-panasas2/skchoudh/genomes/hg38/annotation/gencode.v25.annotation.gtfq\x82a}q\x83(X\n\x00\x00\x00annotationq\x84h\x82h\x08}q\x85h\x84K\x00N\x86q\x86subX\t\x00\x00\x00resourcesq\x87csnakemake.io\nResources\nq\x88)\x81q\x89(K\x01K\x01e}q\x8a(X\x06\x00\x00\x00_coresq\x8bK\x01X\x06\x00\x00\x00_nodesq\x8cK\x01h\x08}q\x8d(h\x8bK\x00N\x86q\x8eh\x8cK\x01N\x86q\x8fuubX\t\x00\x00\x00wildcardsq\x90csnakemake.io\nWildcards\nq\x91)\x81q\x92}q\x93h\x08}q\x94sbX\x03\x00\x00\x00logq\x95csnakemake.io\nLog\nq\x96)\x81q\x97}q\x98h\x08}q\x99sbub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
