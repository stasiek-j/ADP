from collections import defaultdict
from tabulate import tabulate
import pandas as pd
from utils import get_aa


def _get_king(orgs, king):
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
    make_table(get_aa("../data/misc")[0])
    kings = ['Archaea', "Eukaryota", "Viruses", "Bacteria"]
    make_table([_get_king(get_aa(f"../data/{king}")[0], king) for king in kings])
    n_term = get_aa(f"../data/Eukaryota")[1]
    n_term_freq = {k: n_term.count(k)/len(n_term) for k in set(n_term)}
    print(max(n_term_freq.items(), key=lambda x: x[1]))
    make_table(get_aa('../data/PDB/')[0])

