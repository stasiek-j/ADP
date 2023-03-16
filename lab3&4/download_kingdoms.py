import os

import requests as r
from bs4 import BeautifulSoup
from random import shuffle
import wget


def download_kingdom(king: str, path: str, count: int=100):
    ftp_url = f"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/{king}/"

    html_page = r.get(ftp_url).text
    soup = BeautifulSoup(html_page, 'html.parser')
    proteomes = [link.get('href') for link in soup.findAll('a') if link.get('href').startswith('UP')]
    shuffle(proteomes)
    added = 0
    i = 0
    print(f"Downloading {king}:", end=' ')
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
        if not i % 10:
            print('.', end='')
    print(f"\nDownloaded {added} {king}'s proteomes.")


if __name__ == '__main__':
    download_kingdom("Archaea", 'Archaea/')
    download_kingdom('Eukaryota', 'Eukaryota/')
    download_kingdom('Bacteria', 'Bacteria/')
    download_kingdom('Viruses', 'Viruses/')

