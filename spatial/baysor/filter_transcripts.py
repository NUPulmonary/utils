#!/usr/bin/env python3

import argparse
import sys
import pandas as pd


def main():
    # Parse input arguments.
    args = parse_args()

    data_frame = pd.read_csv(args.transcript)

    # Filter transcripts. Ignore negative controls
    filtered_frame = data_frame[(data_frame["qv"] >= args.min_qv) &
                                (data_frame["x_location"] >= args.min_x) &
                                (data_frame["x_location"] <= args.max_x) &
                                (data_frame["y_location"] >= args.min_y) &
                                (data_frame["y_location"] <= args.max_y) &
                                (~data_frame["feature_name"].str.startswith("NegControlProbe_")) &
                                (~data_frame["feature_name"].str.startswith("antisense_")) &
                                (~data_frame["feature_name"].str.startswith("NegControlCodeword_")) &
                                (~data_frame["feature_name"].str.startswith("BLANK_"))]

    # Change cell_id of cell-free transcripts from -1 to 0
    neg_cell_row = filtered_frame["cell_id"] == 'UNASSIGNED'
    filtered_frame.loc[neg_cell_row,"cell_id"] = 0

    # Output filtered transcripts to CSV
    filtered_frame.to_csv("filtered_transcripts.csv",
                          index=False,
                          encoding = 'utf-8')


#--------------------------
# Helper functions

def parse_args():
    """Parses command-line options for main()."""
    summary = 'Filter transcripts from transcripts.csv based on Q-Score threshold \
               and upper bounds on x and y coordinates. Remove negative controls.'

    parser = argparse.ArgumentParser(description=summary)
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-transcript',
                               required = True,
                               help="The path to the transcripts.csv file produced " +
                                    "by Xenium.")
    parser.add_argument('-min_qv',
                        default='20.0',
                        type=float,
                        help="The minimum Q-Score to pass filtering. (default: 20.0)")
    parser.add_argument('-min_x',
                        default='0.0',
                        type=float,
                        help="Only keep transcripts whose x-coordinate is greater than specified limit. " +
                             "If no limit is specified, the default minimum value will be 0.0")
    parser.add_argument('-max_x',
                        default='24000.0',
                        type=float,
                        help="Only keep transcripts whose x-coordinate is less than specified limit. " +
                             "If no limit is specified, the default value will retain all " +
                             "transcripts since Xenium slide is <24000 microns in x and y. " +
                             "(default: 24000.0)")
    parser.add_argument('-min_y',
                        default='0.0',
                        type=float,
                        help="Only keep transcripts whose y-coordinate is greater than specified limit. " +
                             "If no limit is specified, the default minimum value will be 0.0")
    parser.add_argument('-max_y',
                        default='24000.0',
                        type=float,
                        help="Only keep transcripts whose y-coordinate is less than specified limit. " +
                             "If no limit is specified, the default value will retain all " +
                             "transcripts since Xenium slide is <24000 microns in x and y. " +
                             "(default: 24000.0)")

    try:
        opts = parser.parse_args()
    except:
        sys.exit(0)

    return opts



if __name__ == "__main__":
    main()

