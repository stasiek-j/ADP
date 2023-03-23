import gzip
from os import listdir
from os.path import isfile, join
from matplotlib import pyplot as plt
import numpy as np
from Bio import SeqIO
from scipy import stats
from typing import Tuple


def get_length_org(mypath: str) -> list:
    """
    Calculates lengths of proteins of organisms in given directory
    :param mypath: path to directory with proteomes
    :return: list of tuples of form (filename, mean_length, sem)
    """
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

    orgs = []
    for file in onlyfiles:
        lens = []
        if file.endswith('gz'):
            with gzip.open(f"{mypath}/{file}", "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.seq:
                        lens.append(len(record.seq))
        elif file.endswith('fasta'):
            print(file)
            for record in SeqIO.parse(f"{mypath}/{file}", "fasta"):
                if record.seq:
                    lens.append(len(record.seq))

        orgs.append((lens, file))

    ret = [(org[0], np.mean(org[0]), stats.sem(org[0]), org[1]) for org in orgs]

    return ret


def get_length_king(orgs: list) -> Tuple:
    """
    Calculates mean length of protein in given list of organisms
    :param orgs: list of form as returned by `get_length_org`
    :return: tuple of mean and sem
    """
    ret = [x[0] for x in orgs]
    ret = [item for sublist in ret for item in sublist]

    return np.mean(ret), stats.sem(ret)


def plot_kingdoms(path: str):
    """
    Plots mean length of proteins in given directory
    :param path: path to directory for which mean protein lengths should be calculated
    :return: Nothing
    """
    data = {key: get_length_org(path + key) for key in ['Archaea', "Eukaryota", "Viruses", "Bacteria"]}
    kings = {key: get_length_king(value) for key, value in data.items()}
    print([x[1] for x in kings.values()])
    plt.bar(range(len(kings)), [x[0] for x in kings.values()],
            yerr=[x[1] for x in kings.values()], tick_label=[key for key, _ in kings.items()])
    plt.show()


def plot_databases(path: str):
    """
    Plots mean length of proteins for each file in given directory
    :param path: path to directory for which mean protein lengths should be calculated
    :return: Nothing
    """
    data = get_length_org(path)
    print(data[0])
    plt.bar(range(len(data)), [x[1] for x in data],
            yerr=[x[2] for x in data], tick_label=[x[3].split('.')[0] for x in data])
    plt.show()


if __name__ == '__main__':
    data_path = '../data/'
    plot_databases(data_path + 'misc')
    plot_databases(data_path + 'databases')
    plot_kingdoms(data_path)

