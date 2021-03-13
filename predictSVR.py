#!/usr/bin/env python3

import pickle
from argparse import ArgumentParser

import pandas as pd

import extract_features


def get_args():
    parser = ArgumentParser(
        description="Starting with a VCF and a serialized trained SVR model, extract features and predict the level of contamination."
    )
    parser.add_argument(
        "--model",
        help="Serialized/pickled SVR model to use for prediction",
        metavar="FILE",
        type=str,
    )
    parser.add_argument(
        "--vcf",
        help="VCF file from which to extract features and predict contamination",
        metavar="FILE",
        type=str,
    )
    parser.add_argument("--index", help="Index file for VCF", metavar="FILE", type=str)
    parser.add_argument(
        "--loci",
        help="BED file containing loci over which to extract features",
        metavar="FILE",
        type=str,
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    svr_model = pickle.load(open(args.model, "rb"))
    extracted_feature_dict = extract_features.extract(args.vcf, args.index, args.loci)
    extracted_feature_df = pd.DataFrame(data=extracted_feature_dict, index=[0])
    ordered_extracted_feature_df = extracted_feature_df[
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

    predicted_contamination = svr_model.predict(ordered_extracted_feature_df)
    print(predicted_contamination)
