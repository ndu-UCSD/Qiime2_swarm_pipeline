#!/usr/bin/env python
import argparse,os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--THREADS', type=str, default = '1', help = 'Number of threads')
    parser.add_argument('-r', '--reference', type=str)
    parser.add_argument('-i', '--inputfile', type=str)
    parser.add_argument('-o', '--outputfile', type=str)
    opts = parser.parse_args()

    THREADS = opts.THREADS
    reference = opts.reference
    input_file = opts.inputfile
    output_file = opts.outputfile

    cmd = 'source activate qiime2-2018.6 && qiime feature-classifier classify-sklearn --p-n-jobs ' +THREADS+' --i-classifier '+reference+' --i-reads '+input_file+' --o-classification '+output_file
    
    os.system(cmd)
    #print(cmd)
        
if __name__ == "__main__":
    main()