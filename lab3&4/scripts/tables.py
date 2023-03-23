import gzip
from collections import defaultdict
from os import listdir
from os.path import isfile, join
from typing import Any, Tuple, List

from tabulate import tabulate
import numpy as np
import pandas as pd
from Bio import SeqIO
from numpy import ndarray
from scipy import stats


def get_aa(mypath: str) -> tuple[list[tuple[list[int], ndarray, ndarray | float, str]], defaultdict[Any, int]]:
    """
    Get amino acid counts and protein lengths from files in `mypath`
    :param mypath: path to the files
    :return: List of tuples of form: (file, aa_dict, lens)
    """
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

    orgs = []
    for file in onlyfiles:
        lens = []
        aa_dict = defaultdict(float)
        if file.endswith('gz'):
            with gzip.open(f"{mypath}/{file}", "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.seq:
                        lens.append(len(record.seq))
                        for letter in record.seq:
                            aa_dict[letter] += 1
        elif file.endswith('fasta'):
            print(file)
            for record in SeqIO.parse(f"{mypath}/{file}", "fasta"):
                if record.seq:
                    lens.append(len(record.seq))
                    for letter in record.seq:
                        aa_dict[letter] += 1
        aa_dict["Avg. protein length"] = np.mean(lens).item()
        orgs.append((file, aa_dict, lens))

    return orgs


def get_king(orgs, king):
    """
    Agregate info from files in directory to one
    :param orgs:
    :param king:
    :return:
    """
    aa = defaultdict(float)
    tot_len = 0
    count_seq = 0
    for org in orgs:
        file, aa_dict, lens = org
        for k, v in aa_dict.items():
            if k != "Avg. protein length":
                aa[k] += v
        tot_len += sum(lens)
        count_seq += len(lens)
    aa["Avg. protein length"] = tot_len / count_seq

    return king, aa, [tot_len]


def make_table(orgs):
    """
    Create table
    :param orgs:
    :return:
    """
    orgs_c = []
    for org in orgs:
        file, aa_dict, lens = org
        aa_dict = {k: v / sum(lens) if k != "Avg. protein length" else v for k, v in aa_dict.items()}
        orgs_c.append((file, aa_dict))
    table = pd.DataFrame([x[1] for x in orgs_c], index=[x[0] for x in orgs_c])
    print(tabulate(table.T.dropna(), headers='keys', tablefmt='psql', floatfmt='.3f'))


if __name__ == '__main__':
    # make_table(get_aa("../data/misc"))
    kings = ['Archaea', "Eukaryota", "Viruses", "Bacteria"]
    make_table([get_king(get_aa(f"../data/{king}"), king) for king in kings])
