from argparse import ArgumentParser

import numpy as np
from matplotlib import pyplot as plt
from utils import get_aa, get_length_org, get_length_king
import pandas as pd


def plot_kingdoms(path: str, swiss: bool = False, stat='mean', save=None):
    """
    Plots mean length of proteins in given directory
    :param swiss:
    :param path: path to directory for which mean protein lengths should be calculated
    :return: Nothing
    """
    kings = ['Archaea', "Eukaryota", "Viruses", "Bacteria", "databases/SwissProt"] if swiss else ['Archaea',
                                                                                                  "Eukaryota",
                                                                                                  "Viruses", "Bacteria"]
    data = {key.split('/')[-1]: get_length_org(path + key, statistic=stat) for key in kings}
    kings = {key: get_length_king(value, statistic=stat) for key, value in data.items()}
    plt.bar(range(len(kings)), [x[0] for x in kings.values()],
            yerr=[x[1] for x in kings.values()], tick_label=[key for key, _ in kings.items()])
    plt.title(f"Avg. proetin lengths using {stat}")
    if save:
        plt.savefig(save)
    plt.show()


def plot_databases(path: str, stat='mean', save=None):
    """
    Plots mean length of proteins for each file in given directory
    :param path: path to directory for which mean protein lengths should be calculated
    :return: Nothing
    """
    data = get_length_org(path, recursive=True, statistic=stat)
    plt.bar(range(len(data)), [x[1] for x in data],
            yerr=[x[2] for x in data], tick_label=[".".join(x[3].split('.')[:-1]) for x in data])
    plt.title(f"Average protein length using {stat}")
    if save:
        plt.savefig(save)
    plt.show()


def pretty_plot_lengths(path: str, boxplot: bool = False, stat='mean', colors=None, save=None):

    if colors is None:
        colors = {"E.coli": ['r', "bacteria"], "human": ['b', "Vertebraes"], "yeast": ['g', "Fungi"],
                  "A. thaliana": ['y', "Plants"],
                  "D. melanogaster": ['k', "Invertebraes"], "C. elegans": ['k', "Invertebraes"],
                  "Mouse": ['b', "Vertebraes"],
                  "D. rerio": ['b', "Vertebraes"], "B. subtilis": ['r', "bacteria"]}
    plt.figure(figsize=(6, 10))
    if not boxplot:
        data = get_length_org(path, statistic=stat)
        plt.bar(range(len(data)), [x[1] for x in data],
                yerr=[x[2] for x in data],
                label=[colors[".".join(x[3].split('.')[:-1])][1] for x in data],
                color=[colors[".".join(x[3].split('.')[:-1])][0] for x in data],
                tick_label=[".".join(x[3].split('.')[:-1]) for x in data]
                )
        plt.xticks(rotation=30)
    if boxplot:
        aa_tuple, _ = get_aa(path, statistic=stat)
        data = [(file, lens) for (file, _, lens) in aa_tuple]
        plt.boxplot([x[1] for x in data], bootstrap=1000)
        plt.xticks(ticks=range(1, len(data) + 1), labels=[".".join(x[0].split('.')[:-1]) for x in data], rotation=30)
    plt.xlabel("Organism")
    plt.ylabel("Avg. protein length")
    plt.title(f"Averege protein length using {stat}")
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    if save:
        plt.savefig(save)
    plt.show()


def pretty_plot_lengths_kings(path: str, boxplot: bool = False, stat='mean', swiss=False, save=None):

    kings = ['Archaea', "Eukaryota", "Viruses", "Bacteria", "databases/SwissProt"] if swiss else ['Archaea',
                                                                                                  "Eukaryota",
                                                                                                  "Viruses", "Bacteria"]

    data = {key.split('/')[-1]: get_length_org(path + key, statistic=stat) for key in kings}
    kings = {key: get_length_king(value, statistic=stat) for key, value in data.items()}

    plt.bar(range(len(kings)), [x[0] for x in kings.values()],
            yerr=[x[1] for x in kings.values()], tick_label=[key for key, _ in kings.items()],
            label=[key for key, _ in kings.items()]
            )

    plt.xlabel("Organism")
    plt.ylabel("Avg. protein length")
    plt.title(f"Averege protein length in kingdoms using {stat}")
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    if save:
        plt.savefig(save)
    plt.show()


def plot_aa(path: str, stat: str = 'mean', save=None):
    aa_tuple, _ = get_aa(path, statistic=stat)
    to_plot = ['yeast.fasta', 'human.fasta', 'E.coli.fasta']

    aa_tuple = [x for x in aa_tuple if x[0] in to_plot]
    df = pd.DataFrame([x[1] for x in aa_tuple],
                      index=[x[0] for x in aa_tuple]).T.dropna().drop("Avg. protein length")

    df = df / df.sum(axis=0)
    df.plot.bar(xlabel="Aminoacid", ylabel="Frequency")
    plt.xticks(rotation=90)
    plt.title("Aminoacid frequency plot for selected organisms")
    if save:
        plt.savefig(save)
    plt.show()


def plot_histogram(path: str, threshold: int = 3000, stat='mean', save=None):
    data = get_aa(path, statistic=stat)[0]
    for i, (file, _, lens) in enumerate(data):
        plt.hist(lens, bins=range(0, threshold, threshold // 10))
        plt.title(f"{file} histogram")
        plt.xlabel("Protein length")
        plt.ylabel("# of proteins")
        if save:
            plt.savefig(save+file.split(".")[0])
        plt.show()


if __name__ == '__main__':
    parser = ArgumentParser(prog="plots",
                            description="generate plots")
    parser.add_argument('data_path', action='store', help="where should data be sotred")
    parser.add_argument("--fa", help="first a", action="store_true")
    parser.add_argument("--fc", help="first c", action="store_true")
    parser.add_argument("--fbc", help="first for b and c", action="store_true")
    parser.add_argument("--sa", help="second a", action="store_true")
    parser.add_argument("--sc", help="second c", action="store_true")
    parser.add_argument("--sd", help='second d', action="store_true")
    parser.add_argument("--save_to")
    args = parser.parse_args()
    data_path = args.data_path

    save_to = args.save_to if args.save_to else None
    if args.fbc:
        plot_databases(data_path + 'databases', save=save_to+"fbc")

    if args.fc:
        plot_kingdoms(data_path, True, save=save_to+'fc')

    if args.sa:
        pretty_plot_lengths(data_path + 'misc', save=save_to+'sa1')
        plot_aa(data_path + 'misc', save=save_to + 'sa2')

    if args.sc:
        pretty_plot_lengths_kings(data_path, save=save_to + 'sc')

    if args.sd:
        plot_histogram(data_path + 'misc', save=save_to + 'sd1')
        pretty_plot_lengths(data_path + 'misc', True, save=save_to + 'sd2')
        pretty_plot_lengths(data_path + 'misc', stat='median', save=save_to + 'sd3')
