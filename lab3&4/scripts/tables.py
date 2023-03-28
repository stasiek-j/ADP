from argparse import ArgumentParser
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


def make_table(orgs, format='psql', save=None):
    """
    Create table
    :param orgs: list of form (file, aa_dict, lens)
    :param format: format in which table should be displayed & optionally saved
    :param save: where to save table, if None then it's not saved, only printed, else only saved
    :return:
    """
    orgs_c = []
    for org in orgs:
        file, aa_dict, lens = org
        aa_dict = {k: v / sum(lens) if k != "Avg. protein length" else v for k, v in aa_dict.items()}
        orgs_c.append((file, aa_dict))
    table = pd.DataFrame([x[1] for x in orgs_c], index=[x[0] for x in orgs_c])
    table = tabulate(table.T.dropna(), headers='keys', tablefmt=format, floatfmt='.3f')
    if not save:
        print(table)
    if save:
        with open(save, 'w+') as f:
            f.write(table)


if __name__ == '__main__':
    parser = ArgumentParser(prog="tables",
                            description="generate tables")
    parser.add_argument('data_path', action='store', help="where should data be sotred")
    parser.add_argument("--fa", help="first a and second a it's the same result", action="store_true")
    parser.add_argument("--fc", help="first c", action="store_true")
    parser.add_argument("--fnt", help="first for b and c", action="store_true")
    parser.add_argument("--sb", help="second b", action="store_true")
    parser.add_argument("--format", help="format of the tables", default='psql')
    parser.add_argument("--save_to", help="Where to store tables")
    args = parser.parse_args()
    data_path = args.data_path
    FORMAT = args.format


    if args.fa:
        make_table(get_aa("%smisc" % data_path)[0], FORMAT, args.save_to + 'fa')
    if args.fc:
        kings = ['Archaea', "Eukaryota", "Viruses", "Bacteria"]
        make_table([_get_king(get_aa(f"%s{king}" % data_path)[0], king) for king in kings], FORMAT, args.save_to + 'fc')

    if args.fnt:
        n_term = get_aa(f"%sEukaryota" % data_path)[1]
        n_term_freq = {k: n_term.count(k)/len(n_term) for k in set(n_term)}
        print(f" Most frequent pair: {max(n_term_freq.items(), key=lambda x: x[1])}")

    if args.sb:
        make_table(get_aa(f'{data_path}databases/PDB/')[0], FORMAT, args.save_to + 'sb')

