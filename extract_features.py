#!/usr/bin/env python3

from argparse import ArgumentParser
from statistics import mean
from statistics import median

import numpy as np
from pysam import VariantFile

import estimateMixtureModel


def get_args():
    # Pass in the VCF and list of sites to evaluate
    parser = ArgumentParser(
        description="Characterize the allelic balance of heterozygotes at a collection of intermediate frequency SNPs."
    )
    parser.add_argument(
        "--vcf",
        help="VCF (including path). Should be bgzipped and tabix indexed",
        metavar="FILE",
        type=str,
    )
    parser.add_argument("--index", type=str, help="Index file for VCF")
    parser.add_argument(
        "--tag", type=str, help="BED file of tag SNPs over which to characterize AB"
    )
    # parser.add_argument(
    #    "--mafs", type=str, help="Output file to write the het MAFs"
    # )
    # parser.add_argument(
    #    "--allhetmafs", type=str, help="Output file to write the het MAFs"
    # )
    # parser.add_argument(
    #    "--homozygousmafs", type=str, help="Output file to write the het MAFs"
    # )
    # parser.add_argument(
    #    "--allmafs", type=str, help="Output file to write the het MAFs"
    # )
    # parser.add_argument(
    #    "--out",
    #    help="Output file to which min AB fractions will be written.",
    #    metavar="FILE",
    #    type=str,
    # )
    return parser.parse_args()


def extract(vcf, index, tag):
    vcf_in = VariantFile(vcf, index_filename=index)
    tagBED = open(tag, "r")
    allmafs = []

    total_depth = 0
    badAB = 0
    min_AD_zero = 0
    min_AD_nonzero = 0
    numSites = 0
    homozygousSites = 0
    homozygous11 = 0
    homozygous00 = 0
    nonzeroAltInHomozygous00 = 0
    nonzeroREFInHomozygous11 = 0
    heterozygousSites = 0
    heterozygousMAFbelow10p = 0
    heterozygousMAFbelow15p = 0
    heterozygousMAFbelow20p = 0
    heterozygousMAFbelow30p = 0
    badHetFracList = []
    badHomFracList = []

    # Iterate over the bed file of positions to characterize
    for record in tagBED:
        # Convert line to the fields we need for pysam VariantFile fetch
        fields = record.rstrip("\n").split("\t")
        for rec in vcf_in.fetch(fields[0], int(fields[1]), int(fields[2])):
            numSites += 1
            total_depth += rec.samples[0]["DP"]
            for sample in rec.samples.items():
                # Need to check if the AD field exists first
                if sample[1].get("AD"):
                    adFields = sample[1].get("AD")
                    if min(adFields) == 0:
                        min_AD_zero += 1
                    if sum(adFields) == 0:
                        continue
                else:
                    continue
                if sample[1].get("GT")[0] == sample[1].get("GT")[1]:
                    if sample[1].get("GT")[0] == 1 and sample[1].get("GT")[1] == 1:
                        homozygous11 += 1
                        if adFields[0] != 0:
                            nonzeroREFInHomozygous11 += 1
                            # print(min(adFields) / sum(adFields), file=homozygousmafsFile)
                            # print(min(adFields) / sum(adFields), file=allmafsFile)
                            allmafs.append(min(adFields) / sum(adFields))
                            badHomFracList.append(min(adFields) / sum(adFields))
                    elif sample[1].get("GT")[0] == 0 and sample[1].get("GT")[1] == 0:
                        homozygous00 += 1
                        if adFields[1] != 0:
                            nonzeroAltInHomozygous00 += 1
                            # print(min(adFields) / sum(adFields), file=homozygousmafsFile)
                            # print(min(adFields) / sum(adFields), file=allmafsFile)
                            allmafs.append(min(adFields) / sum(adFields))
                            badHomFracList.append(min(adFields) / sum(adFields))
                    homozygousSites += 1
                else:
                    heterozygousSites += 1
                    # print(min(adFields)/sum(adFields), file=mafsFile)
                    # print(min(adFields)/sum(adFields), file=allmafsFile)
                    allmafs.append(min(adFields) / sum(adFields))
                    if min(adFields) / sum(adFields) < 0.3:
                        heterozygousMAFbelow30p += 1
                        badHetFracList.append(min(adFields) / sum(adFields))
                    if min(adFields) / sum(adFields) < 0.2:
                        heterozygousMAFbelow20p += 1
                    if min(adFields) / sum(adFields) < 0.15:
                        heterozygousMAFbelow15p += 1
                    if min(adFields) / sum(adFields) < 0.10:
                        heterozygousMAFbelow10p += 1

    average_depth = total_depth / numSites
    lowABhetFraction = heterozygousMAFbelow30p / heterozygousSites
    veryLowABhetFraction = heterozygousMAFbelow20p / heterozygousSites
    below15ABhetFraction = heterozygousMAFbelow15p / heterozygousSites
    below10ABhetFraction = heterozygousMAFbelow10p / heterozygousSites
    fracAB0 = min_AD_zero / numSites
    fracHomozygousAB0 = min_AD_zero / homozygousSites
    homozygosity = homozygousSites / (homozygousSites + heterozygousSites)
    heterozygosity = heterozygousSites / (homozygousSites + heterozygousSites)
    nonzeroUncalledAlleleSupport = nonzeroAltInHomozygous00 + nonzeroREFInHomozygous11

    qc_results = dict()
    qc_results["average_depth"] = average_depth
    qc_results["homozygosity"] = homozygosity
    qc_results["heterozygosity"] = heterozygosity
    qc_results["fraction_homozygous_minAB0"] = fracHomozygousAB0
    qc_results["fraction_heterozygous_minABfracBelow30"] = lowABhetFraction
    qc_results["fraction_heterozygous_minABfracBelow20"] = veryLowABhetFraction
    qc_results["fraction_heterozygous_minABfracBelow15"] = below15ABhetFraction
    qc_results["fraction_heterozygous_minABfracBelow10"] = below10ABhetFraction
    qc_results["mean_bad_het_MAF"] = mean(badHetFracList)
    qc_results["median_bad_het_MAF"] = median(badHetFracList)
    qc_results["mean_bad_hom_MAF"] = mean(badHomFracList)
    qc_results["median_bad_hom_MAF"] = median(badHomFracList)

    ln_mixture = estimateMixtureModel.get_ln_mixture_model_stats(allmafs)

    qc_results["mu1"] = float(ln_mixture["mu1"])
    qc_results["mu2"] = float(ln_mixture["mu2"])
    qc_results["sd1"] = float(ln_mixture["sd1"])
    qc_results["sd2"] = float(ln_mixture["sd2"])
    qc_results["weights1"] = ln_mixture["weights1"]
    qc_results["weights2"] = ln_mixture["weights2"]
    print("Finished {}".format(vcf), flush=True)
    return qc_results


if __name__ == "__main__":
    args = get_args()
    extract(args.vcf, args.index, args.tag)
