IDS = glob_wildcards("data/genomes/{genomeName}.fna.gz")

rule all:
  input:
    #"data/taxIdAppendedGenomes/taxIdAppended_GCA_001652855.1_Pro_Min_D-127.v1.0_genomic.fna.gz"
    expand("data/taxIdAppendedGenomes/taxIdAppended_{genomeName}.fasta", genomeName = IDS.genomeName)

rule addKrakenTaxidGenome:
  input:
    genome_csv = "data/mmetsp_ncbi_genome_info.csv",
    original_file = "data/genomes/{genomeName}.fna.gz"
  output:
    appended_file = "data/taxIdAppendedGenomes/taxIdAppended_{genomeName}.fasta"
  script:
    "scripts/addKrakenTaxid.py"
  