import pandas as pd
from Bio import SeqIO
import sys
import gzip

data = pd.read_csv(snakemake.input['genome_csv'])
data.set_index("genome_filename", inplace=True)

original_file = snakemake.input['original_file']

genomeFile = snakemake.wildcards['genomeName'] + '.fna.gz'

taxId = str(data.loc[genomeFile]['taxon_id'])

appended_file = snakemake.output['appended_file']

count = 0

with gzip.open(original_file, 'rt') as original, open(appended_file, 'w') as appended:
    records = SeqIO.parse(original, 'fasta')
    
    prefix = '|kraken:taxid|' + taxId
    
    for record in records:          
        record.id = record.id + prefix
        SeqIO.write(record, appended, 'fasta')
        count += 1
    

