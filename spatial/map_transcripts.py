#!/usr/bin/env python3

import argparse
import csv
import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sparse
import scipy.io as sio
import subprocess



def main():
    # Parse input arguments.
    args = parse_args()

    # Check for existence of input file.
    if (not os.path.exists(args.baysor)):
        print("The specified Baysor output (%s) does not exist!" % args.baysor)
        sys.exit(0)

    # Check if output folder already exist.
    if (os.path.exists(args.out)):
        print("The specified output folder (%s) already exists!" % args.out)
        sys.exit(0)


    # Read 5 columns from transcripts Parquet file
    transcripts_df = pd.read_csv(args.baysor,
                                 usecols=["gene",
                                          "cell",
                                          "assignment_confidence"])

    # Find distinct set of features.
    features = np.unique(transcripts_df["gene"])

    # Create lookup dictionary
    feature_to_index = dict()
    for index, val in enumerate(features):
        feature_to_index[str(val)] = index

    # Find distinct set of cells. Discard the first entry which is 0 (non-cell)
    cells = np.unique(transcripts_df["cell"])[1:]


    # Create a cells x features data frame, initialized with 0
    matrix = pd.DataFrame(0, index=range(len(features)), columns=cells, dtype=np.int32)


    # Iterate through all transcripts
    for index, row in transcripts_df.iterrows():
        if index % args.rep_int == 0:
            print(index, "transcripts processed.")

        feature = str(row['gene'])
        cell = row['cell']
        conf = row['assignment_confidence']

        # Ignore transcript below user-specified cutoff
        if conf < args.conf_cutoff:
            continue

        # If cell is not 0 at this point, it means the transcript is associated with a cell
        if cell != 0:
            # Increment count in feature-cell matrix
            matrix.at[feature_to_index[feature], cell] += 1

    # Call a helper function to create Seurat and Scanpy compatible MTX output
    write_sparse_mtx(args, matrix, cells, features)



#--------------------------
# Helper functions

def parse_args():
    """Parses command-line options for main()."""
    summary = 'Map Xenium transcripts to Baysor segmentation result. \
               Generate Seurat/Scanpy-compatible feature-cell matrix.'

    parser = argparse.ArgumentParser(description=summary)
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-baysor',
                               required = True,
                               help="The path to the *segmentation.csv file produced " +
                                    "by Baysor.")
    requiredNamed.add_argument('-out',
                               required = True,
                               help="The name of output folder in which feature-cell " +
                                    "matrix is written.")
    parser.add_argument('-conf_cutoff',
                        default='0.9',
                        type=float,
                        help="Ignore transcripts with assignment confidence " +
                             " below this threshold. (default: 20.0)")
    parser.add_argument('-rep_int',
                        default='100000',
                        type=int,
                        help="Reporting interval. Will print message to stdout " +
                             "whenever the specified # of transcripts is processed. " +
                             "(default: 100000)")

    try:
        opts = parser.parse_args()
    except:
        sys.exit(0)

    return opts



def write_sparse_mtx(args, matrix, cells, features):
    """Write feature-cell matrix in Seurat/Scanpy-compatible MTX format"""

    # Create the matrix folder.
    os.mkdir(args.out)

    # Convert matrix to scipy's COO sparse matrix.
    sparse_mat = sparse.coo_matrix(matrix.values)

    # Write matrix in MTX format.
    sio.mmwrite(args.out + "/matrix.mtx", sparse_mat)

    # Write cells as barcodes.tsv. File name is chosen to ensure
    # compatibility with Seurat/Scanpy.
    with open(args.out + "/barcodes.tsv", 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        for cell in cells:
            writer.writerow(["cell_" + str(cell)])

    # Write features as features.tsv. Write 3 columns to ensure
    # compatibility with Seurat/Scanpy.
    with open(args.out + "/features.tsv", 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        for f in features:
            feature = str(f)
            if feature.startswith("NegControlProbe_") or feature.startswith("antisense_"):
                writer.writerow([feature, feature, "Negative Control Probe"])
            elif feature.startswith("NegControlCodeword_"):
                writer.writerow([feature, feature, "Negative Control Codeword"])
            elif feature.startswith("BLANK_"):
                writer.writerow([feature, feature, "Blank Codeword"])
            else:
                writer.writerow([feature, feature, "Gene Expression"])

    # Seurat expects all 3 files to be gzipped
    subprocess.run("gzip -f " + args.out + "/*", shell=True)



if __name__ == "__main__":
    main()
