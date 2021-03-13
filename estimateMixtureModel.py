#!/usr/bin/env python3

from argparse import ArgumentParser

import pomegranate as pom
import numpy as np
import json


### Functions ###
def get_args():
    parser = ArgumentParser(
        description="Takes an array of less commonly observed allele fractions and returns a lognormal mixture model fitted to the allele fractions"
    )

    parser.add_argument(
        "--mafsfile",
        type=str,
        help="File that contains minor allele fractions, one per line, with no header",
    )
    return parser.parse_args()


def logmu_to_mu(logmu, logsd):
    # Arithmetic mean
    return np.exp(logmu + 0.5 * logsd ** 2.0)


def logsd_to_sd(logmu, logsd):
    # sqrt(Artithmetic variance)
    return np.sqrt((np.exp(logsd ** 2.0) - 1.0) * np.exp(2.0 * logmu + logsd ** 2.0))


def get_ln_mixture_model_stats(mafs):
    # mafs is passed in as a list, so we'll turn it into a numpy array after removing 0s
    mafs = [i for i in mafs if float(i) > 0.0]
    mafs = np.array(mafs)
    # Initialize the mixture model
    d1 = pom.LogNormalDistribution(0, 1)
    d2 = pom.LogNormalDistribution(0.5, 1)

    mixtureMod = pom.GeneralMixtureModel([d1, d2], weights=[0.5, 0.5])

    # Fit the mixture model
    fittedMod = mixtureMod.fit(mafs)
    # Extract the information
    fittedModJSON = json.loads(fittedMod.to_json())
    mulog1 = fittedMod.distributions[0].parameters[0]
    sdlog1 = fittedMod.distributions[0].parameters[1]

    mulog2 = fittedMod.distributions[1].parameters[0]
    sdlog2 = fittedMod.distributions[1].parameters[1]

    results = dict()
    # Translate back into regular space
    results["mu1"] = logmu_to_mu(mulog1, sdlog1)
    results["mu2"] = logmu_to_mu(mulog2, sdlog2)
    results["sd1"] = logsd_to_sd(mulog1, sdlog1)
    results["sd2"] = logsd_to_sd(mulog2, sdlog2)
    results["mu1log"] = mulog1
    results["mu2log"] = mulog2
    results["sd1log"] = sdlog1
    results["sd2log"] = sdlog2
    results["weights1"] = fittedModJSON["weights"][0]
    results["weights2"] = fittedModJSON["weights"][1]
    return results


if __name__ == "__main__":
    args = get_args()
    with open(args.mafsfile, "r") as mafsFH:
        mafslist = mafsFH.read().splitlines()
        mixture_mod_stats = get_ln_mixture_model_stats(mafslist)
        print(mixture_mod_stats)
