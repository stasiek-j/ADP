import re
from argparse import ArgumentParser
from collections import OrderedDict


def read_gp(arg) -> list:
    ret = []
    with open(arg['input'], 'r') as inp:
        file = [x.strip() for x in inp.readlines()]
        # print(file)
        size = len(file)
        idx_list = [idx + 1 for idx, val in
                    enumerate(file) if val == '//']

        res = [file[i: j] for i, j in
               zip([0] + idx_list, idx_list +
                   ([size] if idx_list[-1] != size else []))]
        res = [x for x in res if x != ['']]

    for sublist in res:
        data = OrderedDict()
        data['seq'] = ''
        seq = False
        for line in sublist:
            if 'organism' in arg:
                if line.startswith("ORGANISM"):
                    data['organism'] = " ".join(line.split()[1:3])  # Weź dwa słowa pojawiające się po 'ORGANISM'
            if 'id' in arg:
                if arg['id'] == "GI":  # Pierwsze co się pojawi po 'GI'
                    if line.startswith("VERSION"):
                        matches = re.findall(r"GI:\d+", line)
                        if matches:
                            data["gi"] = matches[0][2:]

                else:
                    if line.startswith("LOCUS"):   # Pierwsze co jest po 'LOCUS'
                        data['locus'] = line.split()[1]

            if arg['additional']:
                ...
            if arg['gene_name']:  # Pierwsze co jest po 'DEFINITION'
                if line.startswith("DEFINITION"):
                    data['gene_name'] = " ".join([x for x in line.split()[1:] if x != "PREDICTED:"])  # Weź dwa słowa pojawiające się po 'ORGANISM'
            # Wszystko co jest po 'ORIGIN' do końca
            if line.startswith("ORIGIN"):
                seq = True
            if seq:
                if line != "//":
                    data['seq'] += ''.join(line.split()[1:])
        data['seq'] = data['seq'].upper()
        ret.append(data)
    return ret


def write_fasta(path: str, data_list: list, sep: str) -> None:
    with open(path, "w+") as f:
        for data in data_list:
            f.write(f">{sep.join([x for k, x in data.items() if k != 'seq'])}\n")
            f.write(data['seq'])
            f.write('\n')


def change_format(data: list, form: str) -> list:
    if form == "M.Musculeu":
        for d in data:
            d["organism"] = f"{d['organism'].split()[0][0]}. {d['organism'].split()[1].lower()}"
    if form == "Musmus":
        for d in data:
            d["organism"] = f"{d['organism'].split()[0][0:3]}{d['organism'].split()[1][0:3].lower()}"
    return data


if __name__ == '__main__':
    parser = ArgumentParser(prog="gp2fasta",
                            description="convert gp to fasta format")
    parser.add_argument("input", action="store", help="Path to the input file")
    parser.add_argument("--output", action="store", default="output.fasta", help="Path to the output file")
    parser.add_argument("--separator", "-s", action="store", default="-", help="Separator to be used in the header")
    parser.add_argument("--organism", "-o", action="store", choices=["Mus musculus", "M.musculeu", "Musmus"],
                        help="Format of organism name to be used in the header")
    parser.add_argument("--id", "-i", action="store", choices=["GI", "LOC"],
                        help="Type of id to be used in the header")
    parser.add_argument("--additional", "-a", action="store_true", help="Should additional info be provided.")
    parser.add_argument("--gene_name", "-g", action="store_true", help="Should the gene name be added in the header.")

    args = parser.parse_args()
    args = vars(args)

    data = read_gp(args)
    data = change_format(data, args["organism"])
    write_fasta(args['output'], data, args['separator'])


