import gzip
from os import listdir
from os.path import isfile, join
from matplotlib import pyplot as plt
from Bio import SeqIO
from scipy import stats


def get_length_org(mypath):
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    orgs = []
    for file in onlyfiles:
        lens = []
        if file.endswith('fasta.gz'):
            with gzip.open(file, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    lens.append(len(record.seq))
        orgs.append(lens)

    ret = [(org, stats.sem(org)) for org in orgs]

    return ret
