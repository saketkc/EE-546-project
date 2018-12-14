#!/usr/bin/env python

'''
FastQC checker for Galaxy biomedical data analysis platform

@author: Ilya Sytchev

Input: one or more files in fastq format
Output: sequencing quality report in text format

Requires FastQC 0.10.0 (http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/)

Partially based on:
fastqcwrapper (http://toolshed.g2.bx.psu.edu/repos/jjohnson/fastqc)
rgFastQC (https://bitbucket.org/galaxy/galaxy-dist/src/tip/tools/rgenetics/rgFastQC.py)

Tested with Python 2.6.1 and 2.7.2 on Mac OS 10.6.8
'''

import sys, os, optparse, tempfile, shutil, subprocess

def stop_err(msg, returncode=1):
    sys.stderr.write(msg)
    sys.exit(returncode)

def __main__():
    usage = "Usage: %prog -e fastqc_executable -o output_file fastq_file [fastq_file ... ]"
    version = "%prog 1.0.0"
    op = optparse.OptionParser(usage=usage, version=version)
    op.add_option('-e', '--executable', dest="executable", help="location of the FastQC program")
    op.add_option('-o', '--output', dest="outfile", help="location of the output file")
    (options, infiles) = op.parse_args()

    # check if location of the FastQC program was provided
    if options.executable == None:
        op.error("Missing location of FastQC")

    # check if FastQC program exists at the provided location
    if not os.path.isfile(options.executable):
        op.error("Cannot find FastQC at %s" % options.executable)

    # check if any input files were provided
    if infiles == None:
        op.error("Missing input files")

    # check if all input files exist
    for f in infiles:
        if not os.path.isfile(f):
            op.error("Cannot find input file %s" % f)
    
    # check if output file was provided
    if options.outfile == None:
        op.error("Missing output file name")
    
    # assemble FastQC command line
    cmd = []    # list is more secure than string for subprocess call 
    cmd.append(options.executable)
    tmpdir = tempfile.mkdtemp()    # create temp dir for FastQC output
    cmd.extend(['-o', tmpdir])
    cmd.extend(infiles)

    # prepare files for FastQC stdout and stderr
    tmp_stderr_name = tempfile.NamedTemporaryFile(dir=tmpdir, suffix='.err').name
    tmp_stderr = open(tmp_stderr_name, 'w')
    tmp_stdout_name = tempfile.NamedTemporaryFile(dir=tmpdir, suffix='.out').name
    tmp_stdout = open(tmp_stdout_name, 'w')
    # run FastQC
    try:
        subprocess.check_call(cmd, stderr=tmp_stderr.fileno(), stdout=tmp_stdout.fileno())
    except subprocess.CalledProcessError as e:
        stop_err("Error executing FastQC\n", e.returncode)
    finally:
        tmp_stderr.close()
        tmp_stdout.close()

    outfile = open(options.outfile, 'w')
    
    # parse all summary.txt files produced by FastQC and write results into the output file
    for f in infiles:
        filename = os.path.basename(f)
        (datasetname, extension) = os.path.splitext(filename)
        # Need to account for FastQC removing .fastq extension from input file names before using them to create output file and dir names
        # Alternative solution is to iterate over report directories instead of input file names
        if extension == '.fastq':
            summaryfilename = os.path.join(tmpdir, datasetname + '_fastqc', 'summary.txt')
        else:
            summaryfilename = os.path.join(tmpdir, filename + '_fastqc', 'summary.txt')
        outfile.write("%s results:\n" % datasetname)
        # if summary file exists, process and add results to the output file
        if os.path.isfile(summaryfilename):
            summaryfile = open(summaryfilename, 'r')
            for line in summaryfile:
                (result, test) = line.split('\t')[:2]
                outfile.write(result + '\t' + test + '\n') 
            summaryfile.close()
        else:
            outfile.write("FastQC summary report was not found at %s.\n" % summaryfilename)
        outfile.write("\n")
    
    outfile.close()

    # clean up temp dir, put in a try block so we don't fail on stale nfs handles
    try: 
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
    except:
        pass

if __name__ == '__main__': __main__()