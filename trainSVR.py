#!/usr/bin/env python3

import sys
from argparse import ArgumentParser
import pickle

import pandas as pd
from sklearn.svm import SVR

import generate_training_data


def get_args():
    parser = ArgumentParser(
        description="Takes training data and fits the SVR model, saving it as a pickled python object."
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--training_data",
        help="File containing extracted feature information for the training data",
        metavar="FILE",
        type=str,
    )
    group.add_argument(
        "--vcf_metadata",
        help="Two column tab-delimited file containing full paths to VCF files and respective contamination fractions",
        metavar="FILE",
        type=str,
    )
    parser.add_argument(
        "--extracted_features",
        help="Output filename to which the extracted features should be written",
        metavar="FILE",
        type=str,
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Name for serialized SVR model output",
        metavar="FILE",
        type=str,
    )
    parser.add_argument(
        "--loci",
        help="BED file containing loci over which to extract the features",
        metavar="FILE",
        type=str,
    )
    return parser.parse_args()


def train_and_save_model(training_data, output_name):
    contamination_fraction = training_data["contamination_fraction"]
    predictor_variables = training_data[
        [
            # "average_depth",
            "fraction_heterozygous_minABfracBelow20",
            "fraction_homozygous_minAB0",
            "heterozygosity",
            "weights1",
            "weights2",
            "mu1",
            "mu2",
        ]
    ]
    svrMod = SVR(gamma="scale", C=1, epsilon=0.1)
    svrMod.fit(predictor_variables, contamination_fraction)
    pickle.dump(svrMod, open(output_name, "wb"))


if __name__ == "__main__":
    args = get_args()
    if args.training_data is not None:
        train_and_save_model(
            pd.read_csv(args.training_data, sep="\t", header=0), args.out
        )
    elif args.vcf_metadata is not None:
        # First make sure that we have an output filename to write the extracted features to
        if args.extracted_features is None:
            sys.exit("Must supply --extracted_features if using --vcf_metadata")
        # Now extract the features from the training VCFs
        training_data_dict = generate_training_data.parse_input_list(args.vcf_metadata)
        # The following stores the extracted features in the dict extracted_features,
        # and also writes the extracted features to the file args.extracted_features
        extracted_features = generate_training_data.extract_features_from_dict(
            training_data_dict, args.extracted_features, args.loci
        )
        train_and_save_model(
            pd.read_csv(args.extracted_features, sep="\t", header=0), args.out
        )
    else:
        print("Must supply either --training_data or --vcf_metadata")
