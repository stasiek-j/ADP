import re
from argparse import ArgumentParser


def read_gp(arg) -> list:
    ret = []
    with open(arg['input'], 'r') as inp:
        file = [x.strip() for x in inp.readlines()]
        ### Part from geeksforgeeks, slightly modified##
        size = len(file)
        idx_list = [idx + 1 for idx, val in
                    enumerate(file) if val == '//']

        res = [file[i: j] for i, j in
               zip([0] + idx_list, idx_list +
                   ([size] if idx_list[-1] != size else []))]
        ################################
        res = [x for x in res if x != ['']]

    for sublist in res:
        data = dict()
        data['seq'] = ''
        seq = False
        defin = False
        for line in sublist:
            if 'organism' in arg:
                if line.startswith("ORGANISM"):
                    data['organism'] = " ".join(line.split()[1:3])  # Weź dwa słowa pojawiające się po 'ORGANISM'
            if 'id' in arg:
                if arg['id'] == "GI":  # Pierwsze co się pojawi po 'GI'
                    if line.startswith("VERSION"):
                        matches = re.findall(r"GI:\d+", line)
                        if matches:
                            data["id"] = matches[0][2:]
                else:
                    if line.startswith("LOCUS"):  # Pierwsze co jest po 'LOCUS'
                        data['id'] = line.split()[1]
            if arg['definition']:  # Pierwsze co jest po 'DEFINITION'
                if defin:
                    if not line.startswith("ACCESSION"):
                        data['definition'] += " ".join([x for x in line.split()[1:] if x != "PREDICTED:"])
                    else:
                        defin = False
                if line.startswith("DEFINITION"):
                    if arg['additional']:
                        data['additional'] = ''.join(['P' if line.find("PREDICTED") != -1 else '',
                                                      's' if line.find('similar') != -1 else '',
                                                      'h' if line.find("hypothetical") != -1 else '',
                                                      'u' if line.find("unnamed protein product") != -1 else '',
                                                      'n' if line.find("novel") != -1 else '',
                                                      'p' if line.find("putative") != -1 else '',
                                                      'o' if line.find("open reading frame") != -1 else ''])
                    data['definition'] = " ".join([x for x in line.split()[1:] if x != "PREDICTED:"])
                    defin = True
            # Sekwencją jest wszystko co jest po 'ORIGIN' do końca
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
            f.write(f">{sep.join([x for k, x in sorted(data.items(), key=lambda x: key(x[0])) if k != 'seq'])}\n")
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


def key(el):
    r = {'organism': 0, 'id': 1, 'definition': 2, 'additional': 3, "sequence": 4, 'seq': 5}
    return r[el]


if __name__ == '__main__':
    parser = ArgumentParser(prog="gp2fasta",
                            description="convert gp to fasta format")
    parser.add_argument("input", action="store", help="Path to the input file")
    parser.add_argument("--output", action="store", default="output.fas", help="Path to the output file")
    parser.add_argument("--separator", "-s", action="store", default="-", help="Separator to be used in the header")
    parser.add_argument("--organism", "-o", action="store", choices=["Mus musculus", "M.musculeu", "Musmus"],
                        help="Format of organism name to be used in the header")
    parser.add_argument("--id", "-i", action="store", choices=["GI", "LOC"],
                        help="Type of id to be used in the header")
    parser.add_argument("--additional", "-a", action="store_true", help="Should additional info be provided.")
    parser.add_argument("--definition", "-d", action="store_true", help="Should the definition be added in the header.")

    args = parser.parse_args()
    args = vars(args)

    data = read_gp(args)
    data = change_format(data, args["organism"])

    write_fasta(args['output'], data, args['separator'])
