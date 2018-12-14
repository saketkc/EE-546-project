
######## Snakemake header ########
import sys; sys.path.append("/home/cmb-panasas2/skchoudh/software_frozen/anaconda27/envs/riboraptor/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x03\x00\x00\x00logq\x03csnakemake.io\nLog\nq\x04)\x81q\x05}q\x06X\x06\x00\x00\x00_namesq\x07}q\x08sbX\t\x00\x00\x00wildcardsq\tcsnakemake.io\nWildcards\nq\n)\x81q\x0b}q\x0ch\x07}q\rsbX\x06\x00\x00\x00configq\x0e}q\x0fX\x0b\x00\x00\x00config_pathq\x10X\x19\x00\x00\x00configs/hg38_SRP055009.pyq\x11sX\x06\x00\x00\x00outputq\x12csnakemake.io\nOutputFiles\nq\x13)\x81q\x14X\x19\x00\x00\x00featureCounts/fcounts.tsvq\x15a}q\x16h\x07}q\x17sbX\x04\x00\x00\x00ruleq\x18X\r\x00\x00\x00featurecountsq\x19X\t\x00\x00\x00resourcesq\x1acsnakemake.io\nResources\nq\x1b)\x81q\x1c(K\x01K\x01e}q\x1d(X\x06\x00\x00\x00_nodesq\x1eK\x01h\x07}q\x1f(X\x06\x00\x00\x00_coresq K\x00N\x86q!h\x1eK\x01N\x86q"uh K\x01ubX\x06\x00\x00\x00paramsq#csnakemake.io\nParams\nq$)\x81q%XN\x00\x00\x00/home/cmb-panasas2/skchoudh/genomes/hg38/annotation/gencode.v25.annotation.gtfq&a}q\'(X\n\x00\x00\x00annotationq(h&h\x07}q)h(K\x00N\x86q*subX\x05\x00\x00\x00inputq+csnakemake.io\nInputFiles\nq,)\x81q-(X\x12\x00\x00\x00bams/SRX876758.bamq.X\x12\x00\x00\x00bams/SRX876771.bamq/X\x12\x00\x00\x00bams/SRX876752.bamq0X\x12\x00\x00\x00bams/SRX876760.bamq1X\x12\x00\x00\x00bams/SRX876756.bamq2X\x12\x00\x00\x00bams/SRX876751.bamq3X\x12\x00\x00\x00bams/SRX876747.bamq4X\x12\x00\x00\x00bams/SRX876736.bamq5X\x12\x00\x00\x00bams/SRX876764.bamq6X\x12\x00\x00\x00bams/SRX876759.bamq7X\x12\x00\x00\x00bams/SRX876740.bamq8X\x12\x00\x00\x00bams/SRX876728.bamq9X\x12\x00\x00\x00bams/SRX876742.bamq:X\x12\x00\x00\x00bams/SRX876765.bamq;X\x12\x00\x00\x00bams/SRX876745.bamq<X\x12\x00\x00\x00bams/SRX876726.bamq=X\x12\x00\x00\x00bams/SRX876743.bamq>X\x12\x00\x00\x00bams/SRX876744.bamq?X\x12\x00\x00\x00bams/SRX876730.bamq@X\x12\x00\x00\x00bams/SRX876761.bamqAX\x12\x00\x00\x00bams/SRX876755.bamqBX\x12\x00\x00\x00bams/SRX876725.bamqCX\x12\x00\x00\x00bams/SRX876750.bamqDX\x12\x00\x00\x00bams/SRX876767.bamqEX\x12\x00\x00\x00bams/SRX876769.bamqFX\x12\x00\x00\x00bams/SRX876732.bamqGX\x12\x00\x00\x00bams/SRX876738.bamqHX\x12\x00\x00\x00bams/SRX876770.bamqIX\x12\x00\x00\x00bams/SRX876741.bamqJX\x12\x00\x00\x00bams/SRX876768.bamqKX\x12\x00\x00\x00bams/SRX876772.bamqLX\x12\x00\x00\x00bams/SRX876727.bamqMX\x12\x00\x00\x00bams/SRX876754.bamqNX\x12\x00\x00\x00bams/SRX876763.bamqOX\x12\x00\x00\x00bams/SRX876762.bamqPX\x12\x00\x00\x00bams/SRX876735.bamqQX\x12\x00\x00\x00bams/SRX876757.bamqRX\x12\x00\x00\x00bams/SRX876724.bamqSX\x12\x00\x00\x00bams/SRX876737.bamqTX\x12\x00\x00\x00bams/SRX876734.bamqUX\x12\x00\x00\x00bams/SRX876723.bamqVX\x12\x00\x00\x00bams/SRX876731.bamqWX\x12\x00\x00\x00bams/SRX876749.bamqXX\x12\x00\x00\x00bams/SRX876748.bamqYX\x12\x00\x00\x00bams/SRX876766.bamqZX\x12\x00\x00\x00bams/SRX876746.bamq[X\x12\x00\x00\x00bams/SRX876733.bamq\\X\x12\x00\x00\x00hdf/SRX876758.hdf5q]X\x12\x00\x00\x00hdf/SRX876771.hdf5q^X\x12\x00\x00\x00hdf/SRX876752.hdf5q_X\x12\x00\x00\x00hdf/SRX876760.hdf5q`X\x12\x00\x00\x00hdf/SRX876756.hdf5qaX\x12\x00\x00\x00hdf/SRX876751.hdf5qbX\x12\x00\x00\x00hdf/SRX876747.hdf5qcX\x12\x00\x00\x00hdf/SRX876736.hdf5qdX\x12\x00\x00\x00hdf/SRX876764.hdf5qeX\x12\x00\x00\x00hdf/SRX876759.hdf5qfX\x12\x00\x00\x00hdf/SRX876740.hdf5qgX\x12\x00\x00\x00hdf/SRX876728.hdf5qhX\x12\x00\x00\x00hdf/SRX876742.hdf5qiX\x12\x00\x00\x00hdf/SRX876765.hdf5qjX\x12\x00\x00\x00hdf/SRX876745.hdf5qkX\x12\x00\x00\x00hdf/SRX876726.hdf5qlX\x12\x00\x00\x00hdf/SRX876743.hdf5qmX\x12\x00\x00\x00hdf/SRX876744.hdf5qnX\x12\x00\x00\x00hdf/SRX876730.hdf5qoX\x12\x00\x00\x00hdf/SRX876761.hdf5qpX\x12\x00\x00\x00hdf/SRX876755.hdf5qqX\x12\x00\x00\x00hdf/SRX876725.hdf5qrX\x12\x00\x00\x00hdf/SRX876750.hdf5qsX\x12\x00\x00\x00hdf/SRX876767.hdf5qtX\x12\x00\x00\x00hdf/SRX876769.hdf5quX\x12\x00\x00\x00hdf/SRX876732.hdf5qvX\x12\x00\x00\x00hdf/SRX876738.hdf5qwX\x12\x00\x00\x00hdf/SRX876770.hdf5qxX\x12\x00\x00\x00hdf/SRX876741.hdf5qyX\x12\x00\x00\x00hdf/SRX876768.hdf5qzX\x12\x00\x00\x00hdf/SRX876772.hdf5q{X\x12\x00\x00\x00hdf/SRX876727.hdf5q|X\x12\x00\x00\x00hdf/SRX876754.hdf5q}X\x12\x00\x00\x00hdf/SRX876763.hdf5q~X\x12\x00\x00\x00hdf/SRX876762.hdf5q\x7fX\x12\x00\x00\x00hdf/SRX876735.hdf5q\x80X\x12\x00\x00\x00hdf/SRX876757.hdf5q\x81X\x12\x00\x00\x00hdf/SRX876724.hdf5q\x82X\x12\x00\x00\x00hdf/SRX876737.hdf5q\x83X\x12\x00\x00\x00hdf/SRX876734.hdf5q\x84X\x12\x00\x00\x00hdf/SRX876723.hdf5q\x85X\x12\x00\x00\x00hdf/SRX876731.hdf5q\x86X\x12\x00\x00\x00hdf/SRX876749.hdf5q\x87X\x12\x00\x00\x00hdf/SRX876748.hdf5q\x88X\x12\x00\x00\x00hdf/SRX876766.hdf5q\x89X\x12\x00\x00\x00hdf/SRX876746.hdf5q\x8aX\x12\x00\x00\x00hdf/SRX876733.hdf5q\x8be}q\x8c(h\x07}q\x8d(X\x04\x00\x00\x00bamsq\x8eK\x00K/\x86q\x8fX\x04\x00\x00\x00hdfsq\x90K/K^\x86q\x91uh\x8ecsnakemake.io\nNamedlist\nq\x92)\x81q\x93(h.h/h0h1h2h3h4h5h6h7h8h9h:h;h<h=h>h?h@hAhBhChDhEhFhGhHhIhJhKhLhMhNhOhPhQhRhShThUhVhWhXhYhZh[h\\e}q\x94h\x07}q\x95sbh\x90h\x92)\x81q\x96(h]h^h_h`hahbhchdhehfhghhhihjhkhlhmhnhohphqhrhshthuhvhwhxhyhzh{h|h}h~h\x7fh\x80h\x81h\x82h\x83h\x84h\x85h\x86h\x87h\x88h\x89h\x8ah\x8be}q\x97h\x07}q\x98sbubX\x07\x00\x00\x00threadsq\x99K\x01ub.'); from snakemake.logging import logger; logger.printshellcmds = True
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
