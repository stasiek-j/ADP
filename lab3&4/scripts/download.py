import json
import os

import requests as r
from bs4 import BeautifulSoup
from random import shuffle
import wget
from tqdm import tqdm
from argparse import ArgumentParser


def download_kingdom(king: str, path: str, count: int = 100):
    ftp_url = f"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/{king}/"

    html_page = r.get(ftp_url).text
    soup = BeautifulSoup(html_page, 'html.parser')
    proteomes = [link.get('href') for link in soup.findAll('a') if link.get('href').startswith('UP')]
    shuffle(proteomes)
    added = 0
    i = 0
    pbar_while = tqdm(desc=f"Downloading {king}", total=count)
    while added < count:
        prot = proteomes[i]
        html_page = r.get(ftp_url + prot).text
        soup = BeautifulSoup(html_page, 'html.parser')
        links = [link.get('href') for link in soup.findAll('a') if link.get('href').startswith('UP') if
                 link.get('href').find('fasta') != -1 and link.get('href').find('DNA.fasta') == -1]
        if len(links) == 1:
            added += 1
            try:
                wget.download(ftp_url + prot + links[0], path + links[0])
            except FileNotFoundError:
                os.mkdir(path)
                wget.download(ftp_url + prot + links[0], path + links[0])

        i += 1
        pbar_while.update()
    return king


def download_query(tax_id: str, path: str, organism: str):
    rest_link = f"https://rest.uniprot.org/proteomes/stream?query=reference:true+taxonomy_id:{tax_id}&top_node&format=list"
    proteome_id = r.get(rest_link).text.split()[0]
    rest_link_prot = f"https://rest.uniprot.org/uniprotkb/stream?query=proteome:{proteome_id}&format=fasta"
    with r.get(rest_link_prot, stream=True) as req:
        req.raise_for_status()
        try:
            with open(path + organism + '.fasta', 'wb') as f:
                for chunk in req.iter_content(chunk_size=2 ** 20):
                    f.write(chunk)
        except FileNotFoundError:
            os.mkdir(path)
            with open(path + organism + '.fasta', 'wb') as f:
                for chunk in req.iter_content(chunk_size=2 ** 20):
                    f.write(chunk)
    return proteome_id


def download_swiss(path: str):
    url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    try:
        wget.download(url, f"{path}swissprot.fasta.gz")
    except FileNotFoundError:
        os.makedirs(path)
        wget.download(url, f"{path}swissprot.fasta.gz")


def download_pdb(path: str):
    url = "https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"
    try:
        wget.download(url, f"{path}pdb.fasta.gz")
    except FileNotFoundError:
        os.makedirs(path)
        wget.download(url, f"{path}pdb.fasta.gz")


if __name__ == '__main__':
    parser = ArgumentParser(prog="download",
                            description="download proteomes")
    parser.add_argument('data_path', action='store', help="where should data be sotred")
    parser.add_argument('--queries', '-q', action='store',
                        help="path to JSON file containing queries to be downloaded",
                        type=str)
    parser.add_argument('--king', '-k', action='store_true', help="should proteomes of organisms in kingdoms given "
                                                                  "in exercise description be downloaded")
    parser.add_argument('--counts', '-c', action='store', type=int, help="how many organisms from each kingdom "
                                                                         "should be downloaded", default=100)
    parser.add_argument('--pdb', '-p', action='store_true', help='should pdb be downloaded')
    parser.add_argument('--swiss', '-s', action='store_true', help='should swissprot be downloaded')

    args = parser.parse_args()

    path = args.data_path
    os.makedirs(path, exist_ok=True)

    if args.queries:
        with open(args.queries, 'r') as f:
            organisms = json.load(f)
        pbar = tqdm(organisms.items())
        for query in pbar:
            pbar.set_description(f"downloaded proteome: {download_query(query[1], f'{path}/misc/', query[0])}")

    if args.king:
        kingdoms = ['Archaea', "Eukaryota", "Viruses", "Bacteria"]
        for kingdom in kingdoms:
            download_kingdom(kingdom, f'{path}/{kingdom}/', args.counts)

    if args.pdb:
        print("Downloading PDB")
        download_pdb(f'{path}/databases/PDB/')

    if args.swiss:
        print("Downloading SwissProt")
        download_swiss(f'{path}/databases/SwissProt/')
