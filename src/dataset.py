"""
Download NCBI dataset of mitochondrial DNA.
Stored as a csv file in 'dat/data.csv', containing species name, NCBI assembly identifier and DNA sequence.

Onno Eberhard, Jan 2022
"""

from zipfile import ZipFile

import ncbi.datasets
import pandas as pd
import requests
from tqdm.auto import tqdm, trange


# Initialize NCBI API connection
api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())

# Search the NCBI database for all listed eukaryotes
genome_summary = api_instance.assembly_descriptors_by_taxon(taxon='eukaryotes', page_size=1000)
eukaryotes = genome_summary.assemblies
page_token = genome_summary.next_page_token
pbar = tqdm(total=genome_summary.total_count)
pbar.update(len(eukaryotes))
while page_token:
    genome_summary = api_instance.assembly_descriptors_by_taxon(
        taxon='eukaryotes', page_size=1000, page_token=page_token)
    a = genome_summary.assemblies
    page_token = genome_summary.next_page_token
    pbar.update(len(a))
    eukaryotes += a
pbar.close()

# Filter out all species without mitochondrium sequence data and build a DataFrame
# of species name and assembly identifier
species = pd.DataFrame([{'id': x['assembly']['assembly_accession'], 'name': x['assembly']['org']['title']}
                        for x in eukaryotes for y in x['assembly']['chromosomes'] if y['name'] == 'MT'])

# Download assemly descriptions (as zip-files), which includes URLs to DNA sequences
for i in trange(0, len(species), 100):
    api_response = api_instance.download_assembly_package(
        species.iloc[i : i + 100].id.to_list(), chromosomes=['MT'], hydrated="DATA_REPORT_ONLY")

    # Write zip file with assembly descriptions
    with open('tmp/assemblies.zip', 'wb') as f:
        f.write(api_response.read())

    # Extract zip file
    with ZipFile('tmp/assemblies.zip', 'r') as f:
        f.extractall('tmp/assemblies')
    
    # Download mitochondrium dna sequences
    with open("tmp/assemblies/ncbi_dataset/fetch.txt") as f:
        for line in f:
            ls = line.split()
            
            # For those assemblies that actually contain a file with the mitochondrium sequence,
            # download and store DNA
            if ls[-1].endswith('chrMT.fna'):
                r = requests.get(ls[0])
                dna = ''.join(r.text.split('\n')[1:])
                index = species[species.id == ls[-1].split('/')[1]].index[0]
                species.loc[index, 'dna'] = dna

# Filter out species with no DNA information and remove duplicates
species = species[~species.dna.isna()]
species = species[species.dna.str.len() > 0]
species = species.drop_duplicates('name')

# Save data as csv file
species.to_csv('dat/data.csv', index=False)
