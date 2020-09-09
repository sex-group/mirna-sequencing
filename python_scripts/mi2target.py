#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Checks the most relevant miRNAs targets and its functions.

Prepared for Vitis Annotation
"""


import argparse
import logging

from argparse import RawTextHelpFormatter
from multiprocessing import Pool

import pandas as pd


__author__ = "Miguel Ramos"
__copyright__ = "Copyright 2020, Miguel Ramos"
__credits__ = ["Miguel Ramos"]
__license__ = "MIT"
__version__ = "1"
__maintainer__ = "Miguel Ramos"
__email__ = "miguelramos22@gmail.com"
__status__ = "Production"


# Logger configurations
log_format = logging.Formatter(("%(asctime)s [%(levelname)s] "
                                "%(name)s@%(lineno)s: %(message)s"))

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler("log.txt")
file_handler.setFormatter(log_format)
file_handler.setLevel(logging.INFO)
logger.addHandler(file_handler)

console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
console_handler.setFormatter(log_format)
logger.addHandler(console_handler)


# Functions -- Alphabetically
def access_mirna(mirna):
    """Infers the miRNA target"""

    logger.info(f"{mirna}: starting...")

    # Search the mirna on psrnatarget
    condition = input_psrnatarget["miRNA_Acc."].str.contains(mirna,
                                                             na=False,
                                                             case=False
                                                             )

    target_subset = input_psrnatarget.loc[condition]
    target_subset = target_subset[["Target_Acc.","Expectation"]]
    target_subset = target_subset.drop_duplicates()

    targets_to_gene = target_subset["Target_Acc."].str.split(".")
    target_subset["Target_Acc."] = [target[0] for target in targets_to_gene]

    target_subset = target_subset.groupby(["Target_Acc."], as_index=False)\
                                 .min()

    targets_list = target_subset["Target_Acc."].tolist()

    # We are not processing this with multiprocessing, because this function is
    #     expected to be a child function of a multiprocessing one!
    target_information = list(map(vit_accession2name_function, targets_list))
    target_information = pd.DataFrame(target_information)

    target_subset = target_subset.set_index("Target_Acc.")
    target_information = target_information.set_index("Target_Acc.")

    target_complete = pd.concat([target_subset, target_information],
                                axis=1,
                                sort=False
                                )

    target_complete = target_complete.sort_values(by = "Expectation")

    target_complete.to_csv(f"m2t_{mirna}.tsv",sep="\t",index_label="Target id")
    logger.info(f"{mirna}: completed")


def mirnas_to_list(mirna_str):
    """Transforms the provided miRNA list from arguments to a python list"""

    logger.info("Loading miRNA list provided as input...")
    mirna_list = [mirna.strip() for mirna in mirna_str.split(',')]
    logger.info(f"miRNAs provided as input: {len(mirna_list)}")

    return mirna_list


def open_annotation(annotation_path):
    """Open the annotation file and load contact to pandas data frame"""

    logger.info("Loading annotation content...")
    annotation = pd.read_csv(annotation_path,
                             delimiter="\t",
                             header=0,
                             )
    logger.info(f"The provided file contains {annotation.shape[0]} lines")

    return annotation


def open_psrnatarget(psrnatarget_path):
    """Open the psrnatarget file and load contact to pandas data frame

    This function expects the first line to be a comment one
    """

    logger.info("Loading psrnatarget content...")
    psrnatarget = pd.read_csv(psrnatarget_path,
                              delimiter="\t",
                              header=0,
                              skiprows=1
                              )
    logger.info(f"The provided file contains {psrnatarget.shape[0]} lines")

    return psrnatarget


def vit_accession2name_function(vit_accession):
    """For a specific vit_accession retrieve name and function from annotation
    """

    condition = input_annotation["GeneID"] == vit_accession
    gene_information = input_annotation.loc[condition]

    vit_information = {
        "Target_Acc.": vit_accession,
        "Gene Name": gene_information["GeneName"].iloc[0],
        "Gene Function": gene_information["GeneFunction"].iloc[0]
        }

    return vit_information


# Run functions on direct call
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--mirnas",
        help="miRNA lists to search for, separated by coma",
        required=True
    )

    parser.add_argument(
        "--psrnatarget",
        help="path to psrnatarget output file",
        required=True
    )

    parser.add_argument(
        "--annotation",
        help="path to Vitis annotation",
        required=True
    )

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s (version {__version__} | {__status__})"
    )

    args = parser.parse_args()


    # Say hi!
    logger.info("mi2target is starting... Hello!")
    logger.info("These were the provided inputs:")
    logger.info(f"--mirnas {args.mirnas}")
    logger.info(f"--psrnatarget {args.psrnatarget}")
    logger.info(f"--annotation {args.annotation}")


    # Run functions from main
    # Load and prepare inputs
    input_mirnas = mirnas_to_list(args.mirnas)
    input_psrnatarget = open_psrnatarget(args.psrnatarget)
    input_annotation = open_annotation(args.annotation)

    # Loop over miRNAs using multiprocessing
    with Pool() as p:
        p.map(access_mirna,input_mirnas)
    #list(map(access_mirna, input_mirnas))
