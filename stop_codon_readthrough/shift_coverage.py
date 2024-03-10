#!/usr/bin/env python
# coding: utf-8


# import packages
import numpy
import pandas



def loadCoverage(result):
    NAMES = ['chr','start','end','transcript','score','strand','tStart','tEnd','rgb','numE','lenE','startE',"index","count"]
    result = pandas.read_csv(result,sep='\t',header=None,names=NAMES)
    return(result)


def shift_index(df):
    if df["index"] <30 or df["count"] == 0:
        return df["index"]
    else:
        shift = df["index"] - numpy.array([int(i) for i in df["startE"].split(",")])

        for i in shift:
            if 0<=i <30:
                return int(i)

# Function to read in command line arguments
def getInputs():
    parser = argparse.ArgumentParser(
        description='Shift index of bedtools perbase coverage files.',
        epilog='Expects a comma-separated list of perbase coverage files.',
        add_help=True,
        allow_abbrev=True
    )
    parser.add_argument(
        '--replications',
        metavar='replications',
        help='Base names of the bedtools perbase coverage files.',
        required=True,
        type=str,
        dest='bed_file'
    )
    parser.add_argument(
        '--input',
        metavar='input',
        help='Path to input directory.',
        required=True,
        type=str,
        dest='stop_codon_entries'
    )
    parser.add_argument(
        '--out',
        metavar='Output directory',
        help='Path to output directory.',
        required=True,
        default=None,
        dest='out'
    )
    return parser.parse_args()

########### shift index for all runs

args = getInputs()

replications = args.replications.split(",")

replications = [name.split("_perbase.txt")[0] for name in replications]


for name in replications:
    print(args.input+name+'_perbase.txt')
    perbase = loadCoverage(args.input+name+'_perbase.txt')
    print("Loaded File!")
    perbase['index']=perbase.apply(shift_index, axis=1)
    perbase.to_csv(args.out+name+'_perbase_shifted.txt', header=None, index=False, sep='\t')
    print("done")
