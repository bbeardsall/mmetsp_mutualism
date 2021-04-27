import pandas as pd
from Bio import SeqIO
import sys

data = pd.read_csv(snakemake.input['genome_csv'])
data.set_index("genome_filename", inplace=True)

# make full organism name with assembly accession for second column
data['full_organism_name'] = data['assembly_accession'] + " " + data['organism_name'] + data['infraspecific_name'].replace("strain=", " ", regex=True).fillna('')

# load snakemake variables
genomeName = snakemake.wildcards['genomeName']
genomeFile = genomeName + ".fna.gz"
maskedFile = snakemake.input['masked_file']
outputFile = snakemake.output[0]

seqIds = []

with open(maskedFile, 'rt') as genome:
    records = SeqIO.parse(genome, 'fasta')

    for record in records:
        seqIds.append(record.id)

outData = pd.DataFrame(list(
    zip(seqIds, [data.loc[genomeFile]['taxon_id']]*len(seqIds), [data.loc[genomeFile]['full_organism_name']]*len(seqIds))
    ),
    columns=['seq', 'taxid', 'name'])

outData.to_csv(outputFile, sep = '\t', header=False, index=False)
    

