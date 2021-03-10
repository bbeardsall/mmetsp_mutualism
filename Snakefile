IDS = glob_wildcards("data/genomes/{genomeName}.fna.gz")
appendedIDS = glob_wildcards("data/taxIdAppendedGenomes/{appended_genomeName}.fasta")

rule all:
  input:
    #"data/taxIdAppendedGenomes/taxIdAppended_GCA_001652855.1_Pro_Min_D-127.v1.0_genomic.fna.gz"
    #expand("data/taxIdAppendedGenomes/taxIdAppended_{genomeName}.fasta", genomeName = IDS.genomeName)
    expand("data/kraken2_add_seq/{appended_genomeName}.out", appended_genomeName = appendedIDS.appended_genomeName)

rule addKrakenTaxidGenome:
  input:
    genome_csv = "data/mmetsp_ncbi_genome_info.csv",
    original_file = "data/genomes/{genomeName}.fna.gz"
  output:
    appended_file = "data/taxIdAppendedGenomes/taxIdAppended_{genomeName}.fasta"
  script:
    "scripts/addKrakenTaxid.py"
    
rule addFastaToDb:
  input:
    appended_file = "data/taxIdAppendedGenomes/{appended_genomeName}.fasta",
    database = "data/mmetsp_genome_db"
  output:
    "data/kraken2_add_seq/{appended_genomeName}.out"
  shell:
    "kraken2-build --add-to-library {input.appended_file} --db {input.database} && touch {output}"
    
  