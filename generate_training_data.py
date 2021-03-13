#!/usr/bin/env python3

from argparse import ArgumentParser
import yaml

import extract_features


def get_args():
    parser = ArgumentParser(
        description="Takes a list of input VCFs with known contamination levels and outputs training data with which to train the SVR"
    )

    parser.add_argument(
        "--inputlist",
        type=str,
        help="""
             Two-column, unheadered file where the first column includes paths to VCF files
             and the second column contains the respective known contamination fractions.
             """,
    )
    parser.add_argument(
        "--loci",
        type=str,
        help="bed file containing the loci over which to extract the features",
    )
    parser.add_argument(
        "--feature_out",
        type=str,
        help="Output file that will contain the extracted features",
    )
    return parser.parse_args()


def parse_input_list(filelist):
    training_dict = {}
    # return(training_dict)
    with open(filelist, "r") as filesFH:
        files = filesFH.read().splitlines()
        for line in files:
            fields = line.split()
            training_dict[fields[0]] = {}
            training_dict[fields[0]]["contamination_fraction"] = fields[1]
    return training_dict


def extract_features_from_dict(sample_dict, outfile, loci):
    """
    Returns a dictionary of dictionaries, where the primary keys are the VCF file names.
    The values are dictionaries that contains the following values:
    average_depth
    homozygosity
    heterozygosity
    fraction_homozygous_minAB0
    fraction_heterozygous_minABfracBelow30
    fraction_heterozygous_minABfracBelow20
    fraction_heterozygous_minABfracBelow15
    fraction_heterozygous_minABfracBelow10
    mean_bad_het_MAF
    median_bad_het_MAF
    mean_bad_hom_MAF
    median_bad_hom_MAF
    distribution parameters: mu1, mu2, sd1, sd2, weights1, weights2
    """
    extracted_feature_dict = {}
    for sample in sample_dict:
        index_file = "{}.tbi".format(sample)
        extracted_feature_dict[sample] = extract_features.extract(
            sample, index_file, loci
        )

    with open(outfile, "w") as outFH:
        print(
            "base\tcontamination_fraction\taverage_depth\thomozygosity\theterozygosity\tfraction_homozygous_minAB0\tfraction_heterozygous_minABfracBelow30\tfraction_heterozygous_minABfracBelow20\tfraction_heterozygous_minABfracBelow15\tfraction_heterozygous_minABfracBelow10\tmean_bad_het_MAF\tmedian_bad_het_MAF\tmean_bad_hom_MAF\tmedian_bad_hom_MAF\tmu1\tmu2\tsd1\tsd2\tweights1\tweights2",
            file=outFH,
        )
        for sample in sample_dict:
            print(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    sample,
                    sample_dict[sample]["contamination_fraction"],
                    extracted_feature_dict[sample]["average_depth"],
                    extracted_feature_dict[sample]["homozygosity"],
                    extracted_feature_dict[sample]["heterozygosity"],
                    extracted_feature_dict[sample]["fraction_homozygous_minAB0"],
                    extracted_feature_dict[sample][
                        "fraction_heterozygous_minABfracBelow30"
                    ],
                    extracted_feature_dict[sample][
                        "fraction_heterozygous_minABfracBelow20"
                    ],
                    extracted_feature_dict[sample][
                        "fraction_heterozygous_minABfracBelow15"
                    ],
                    extracted_feature_dict[sample][
                        "fraction_heterozygous_minABfracBelow10"
                    ],
                    extracted_feature_dict[sample]["mean_bad_het_MAF"],
                    extracted_feature_dict[sample]["median_bad_het_MAF"],
                    extracted_feature_dict[sample]["mean_bad_hom_MAF"],
                    extracted_feature_dict[sample]["median_bad_hom_MAF"],
                    extracted_feature_dict[sample]["mu1"],
                    extracted_feature_dict[sample]["mu2"],
                    extracted_feature_dict[sample]["sd1"],
                    extracted_feature_dict[sample]["sd2"],
                    extracted_feature_dict[sample]["weights1"],
                    extracted_feature_dict[sample]["weights2"],
                ),
                file=outFH,
            )
    return extracted_feature_dict


if __name__ == "__main__":
    args = get_args()

    # Parse the list of files with associated contamination fractions
    input_dict = parse_input_list(args.inputlist)
    # print(yaml.dump(input_dict))
    # Go through each of the VCFs and generate the summary data needed for each
    training_data = extract_features_from_dict(input_dict, args.feature_out, args.loci)

    # print(yaml.dump(training_data))
