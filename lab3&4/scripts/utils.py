import gzip
import os
from collections import defaultdict
from os import listdir
from os.path import isfile, join
from typing import Tuple, List, Any

import numpy as np
from Bio import SeqIO


def get_length_org(mypath: str, recursive: bool = False, statistic: str = 'mean') -> list:
    """
    Calculates lengths of proteins of organisms in given directory
    :param mypath: path to directory with proteomes
    :param recursive: Should the function get files from subdirectories
    :param statistic: what statistic should be used mean or median
    :return: list of tuples of form (filename, mean_length, sem)
    """
    assert statistic in {'mean', "median"}, "Statistic should be one of mean or median"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    if recursive:
        myfiles = []
        for root, subdirs, files in os.walk(mypath):
            for name in files:
                myfiles.append(join(root, name))
        onlyfiles = myfiles
        print(myfiles)
    orgs = []
    for file in onlyfiles:
        lens = []
        if file.endswith('gz'):
            with gzip.open(f"{file}" if recursive else f"{mypath}/{file}", "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.seq:
                        lens.append(len(record.seq))
        elif file.endswith('fasta'):
            for record in SeqIO.parse(f"{file}" if recursive else f"{mypath}/{file}", "fasta"):
                if record.seq:
                    lens.append(len(record.seq))

        orgs.append((lens, file.split("/")[-1]))
    if statistic == 'mean':
        ret = [(org[0], np.mean(org[0]), np.std(org[0])/np.sqrt(len(org[0])), org[1]) for org in orgs]
    else:
        ret = [(org[0], np.median(org[0]), np.std(org[0]) / np.sqrt(len(org[0])), org[1]) for org in orgs]

    return ret


def get_length_king(orgs: list, statistic: str = 'mean') -> Tuple:
    """
    Calculates mean length of protein in given list of organisms
    :param orgs: list of form as returned by `get_length_org`
    :param statistic: what statistic should be used mean or median
    :return: tuple of mean and sem
    """
    assert statistic in {'mean', "median"}, "Statistic should be one of mean or median"
    ret = [x[0] for x in orgs]
    ret = [item for sublist in ret for item in sublist]

    return np.mean(ret) if statistic == 'mean' else np.median(ret), np.std(ret)/np.sqrt(len(ret))


def get_aa(mypath: str, file_flag: bool = False, statistic: str = 'mean') -> Tuple[List[Tuple[str, defaultdict[Any, float], List[int]]], List[Any]]:
    """
    Get amino acid counts and protein lengths from files in `mypath`
    :param mypath: path to the files
    :param file_flag: Is the path a path to the file?
    :param statistic: what statistic should be used mean or median
    :return: List of tuples of form: (file, aa_dict, lens), and list of n terminus aminoacids
    """
    assert statistic in {'mean', "median"}, "Statistic should be one of mean or median"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))] if not file_flag else [mypath.split('/')[-1]]
    if file_flag:
        mypath = "/".join(mypath.split('/')[:-1])

    orgs = []
    n_term = []
    for file in onlyfiles:
        lens = []
        aa_dict = defaultdict(float)
        if file.endswith('gz'):
            with gzip.open(f"{mypath}/{file}", "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.seq:
                        lens.append(len(record.seq))
                        n_term.append(record.seq[0])
                        for letter in record.seq:
                            aa_dict[letter] += 1
        elif file.endswith('fasta'):
            for record in SeqIO.parse(f"{mypath}/{file}", "fasta"):
                if record.seq:
                    lens.append(len(record.seq))
                    n_term.append(record.seq[0])
                    for letter in record.seq:
                        aa_dict[letter] += 1
        if statistic == 'mean':
            aa_dict["Avg. protein length"] = np.mean(lens).item()
        else:
            aa_dict["Avg. protein length"] = np.median(lens).item()
        orgs.append((file, aa_dict, lens))

    return orgs, n_term

