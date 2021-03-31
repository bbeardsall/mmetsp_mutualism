import pandas as pd

genome_data = pd.read_csv("data/mmetsp_ncbi_genome_info.csv")
#/iplant/home/shared/imicrobe/projects/104/transcriptomes/MMETSP005/MMETSP0053.fastq.tar
#fastq_files = genome_data['fastq_file']

_sample_ids = genome_data['sample_id'].tolist()
_sample_names = genome_data['sample_name_main'].tolist()
_genus = genome_data['genus'].tolist()
_species = genome_data['species'].tolist()
_strain = genome_data['strain'].tolist()

IDS = glob_wildcards("data/genomes/{genomeName}.fna.gz")
appendedIDS = glob_wildcards("data/taxIdAppendedGenomes/{appended_genomeName}.fasta")

rule all:
  input:
    "results/krakenOutputs/2475-MMETSP0053-Prorocentrum-minimum-CCMP1329.kraken"
    #expand("data/tar/{sample}/{sample_name_main}.fastq.tar", zip, sample = _sample_ids, sample_name_main = _sample_names)

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
    
rule downloadMmetspTar:
    output:
        "data/tar/{sample}/{sample_name_main}.fastq.tar"
    shell:
        "iget -PT /iplant/home/shared/imicrobe/projects/104/samples/{wildcards.sample}/{wildcards.sample_name_main}.fastq.tar data/tar/{wildcards.sample}"
    
checkpoint untarFastq:
  input:
    "data/tar/{sample}/{sample_name_main}.fastq.tar"
  output:
    directory("data/MMETSP_fastq/{sample}-{sample_name_main}")
  shell:
    "mkdir -p data/MMETSP_fastq/{wildcards.sample}-{wildcards.sample_name_main} && tar xfv {input} -C data/MMETSP_fastq/{wildcards.sample}-{wildcards.sample_name_main}"

def aggregate_untarFastq(wildcards):
  checkpoint_output = checkpoints.untarFastq.get(**wildcards).output[0]
  return expand("results/krakenOutputs/{sample}-{sample_name_main}-{fullSpeciesName}.kraken",
  sample = wildcards.sample,
  sample_name_main = wildcards.sample_name_main,
  i = [1,2],
  fullSpeciesName = glob_wildcards(os.path.join(checkpoint_output, "{sample_name_main}-{fullSpeciesName}.1.fastq.gz")).fullSpeciesName),
  
rule kraken2:
    input:
        fastq1 = "data/MMETSP_fastq/{sample}/{sample_name_main}-{fullSpeciesName}.1.fastq.gz",
        fastq2 = "data/MMETSP_fastq/{sample}/{sample_name_main}-{fullSpeciesName}.2.fastq.gz"
    output:
        krakenOut = "results/krakenOutputs/{sample}-{sample_name_main}-{fullSpeciesName}.kraken",
        reportOut = "results/krakenReports/{sample}-{sample_name_main}-{fullSpeciesName}.report"
    shell:
        "kraken2 --use-names --threads 32 --db data/mmetsp_genome_db --report {output.reportOut} --paired {input.fastq1} {input.fastq2} > {output.krakenOut}"

rule finished:
  input: aggregate_untarFastq
  output: "results/finished.txt"
  shell: "touch {output}"

rule getMmetspFastq:
  output:
    "data/MMETSP_fastq/{sample}/{sample_name_main}-{genus}-{species}-{strain}.1.fastq",
    "data/MMETSP_fastq/{sample}/{sample_name_main}-{genus}-{species}-{strain}.2.fastq"
  shell:
    "iget -PT /iplant/home/shared/imicrobe/projects/104/samples/{wildcards.sample}/{wildcards.sample_name_main}-{wildcards.genus}-{wildcards.species}-{wildcards.strain}.1.fastq data/MMETSP_fastq/{wildcards.sample} && "
    "iget -PT /iplant/home/shared/imicrobe/projects/104/samples/{wildcards.sample}/{wildcards.sample_name_main}-{wildcards.genus}-{wildcards.species}-{wildcards.strain}.2.fastq data/MMETSP_fastq/{wildcards.sample}"
    
