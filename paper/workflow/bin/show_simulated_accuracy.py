#!/usr/bin/env python
import argparse
import pandas as pd
from common import fileopen
from compute_simulated_accuracy import summary_fields

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default="-")
    parser.add_argument("-o", "--output")
    parser.add_argument("-n", "--name")
    parser.add_argument("-c", "--caption")
    args = parser.parse_args()
    
    with fileopen(args.input, 'rt') as inp:
        table = pd.read_csv(inp, sep="\t", names=summary_fields)
    
    import pickle
    with fileopen(args.output + '.pickle', 'wb') as out:
        pickle.dump(table, out)

if __name__ == "__main__":
    main()