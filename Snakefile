import pandas as pd

genome_data = pd.read_csv("data/mmetsp_ncbi_genome_info.csv")
#/iplant/home/shared/imicrobe/projects/104/transcriptomes/MMETSP005/MMETSP0053.fastq.tar
#fastq_files = genome_data['fastq_file']
_genome_names = genome_data['genome_name'].tolist()

_sample_ids = genome_data['sample_id'].tolist()
#_sample_ids = [2475]
_sample_names = genome_data['sample_name_main'].tolist()
#_sample_names = ["MMETSP0053"]
# _genus = genome_data['genus'].tolist()
# _species = genome_data['species'].tolist()
# _strain = genome_data['strain'].tolist()

# IDS = glob_wildcards("data/genomes/{genomeName}.fna.gz")
# appendedIDS = glob_wildcards("data/taxIdAppendedGenomes/{appended_genomeName}.fasta")

rule all:
  input:
    #"results/krakenOutputs/2475_MMETSP0053_Prorocentrum-minimum-CCMP1329.kraken"
    # expand("data/tar/{sample}/{sample_name_main}.fastq.tar", zip, sample = _sample_ids, sample_name_main = _sample_names)
    #expand("results/krakenReports/{sample_id}_-_{sample_name_main}_-_{fullSpeciesName}.report", zip, sample_id = _sample_ids, sample_name_main = _sample_names)
    expand("results/done_krona/{sample_id}_-_{sample_name_main}.done", zip, sample_id = _sample_ids, sample_name_main = _sample_names)
    #expand("data/krakenUniq_mmetsp_genome_db/library/added/{genomeName}.fna.map", genomeName = _genome_names)

### Preparing kraken2 database ###

rule addKrakenTaxidGenome:
  input:
    genome_csv = "data/mmetsp_ncbi_genome_info.csv",
    original_file = "data/genomes/{genomeName}.fna.gz"
  output:
    appended_file = "data/taxIdAppendedGenomes/taxIdAppended_-_{genomeName}.fasta"
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

### Preparing krakenuniq database ###
rule unzipGenome:
  input:
    "data/genomes/{genomeName}.fna.gz"
  output:
    temp("data/unzipped_genomes/{genomeName}.fna")
  shell:
    "gzip -dc {input} > {output}"

rule dustmaskGenome:
  input:
    "data/unzipped_genomes/{genomeName}.fna"
  output:
    "data/krakenUniq_mmetsp_genome_db/library/added/{genomeName}-dustmasked.fna"
  shell:
    "dustmasker -infmt fasta -in {input} -level 20 -outfmt fasta | sed '/^>/! s/[^AGCT]/N/g' > {output}"

rule mapGenomeSequences:
  input:
    masked_file = "data/krakenUniq_mmetsp_genome_db/library/added/{genomeName}-dustmasked.fna",
    genome_csv = "data/mmetsp_ncbi_genome_info.csv"
  output:
    "data/krakenUniq_mmetsp_genome_db/library/added/{genomeName}.fna.map"
  script:
    "scripts/map_krakenuniq_taxonomy.py"

### Downloading MMETSP and running Kraken ###   

# download tar of fastq files from cyverse
rule downloadMmetspTar:
    output:
        "data/tar/{sample_id}/{sample_name_main}.fastq.tar"
    shell:
        "iget -PT /iplant/home/shared/imicrobe/projects/104/samples/{wildcards.sample_id}/{wildcards.sample_name_main}.fastq.tar data/tar/{wildcards.sample_id}"

# extract fastq files
checkpoint untarFastq:
  input:
    "data/tar/{sample_id}/{sample_name_main}.fastq.tar"
  output:
    directory("data/MMETSP_fastq/{sample_id}_-_{sample_name_main}")
  shell:
    "mkdir -p data/MMETSP_fastq/{wildcards.sample_id}_-_{wildcards.sample_name_main} && tar xfv {input} -C data/MMETSP_fastq/{wildcards.sample_id}_-_{wildcards.sample_name_main}"

# checkpoint function to get names of fastq
def aggregate_untarFastq(wildcards):
  checkpoint_output = checkpoints.untarFastq.get(**wildcards).output[0]

  # includes MMETSP name along with species name
  fullSpeciesNames, Index = glob_wildcards(os.path.join(checkpoint_output, "{fullSpeciesName}.{i}.fastq.gz"))
  
  # kraken_filenames = expand("results/krakenOutputs/{sample_id}_-_{sample_name_main}_-_{fullSpeciesName}.kraken", zip,
  # sample_id = wildcards.sample_id,
  # sample_name_main = wildcards.sample_name_main,
  # fullSpeciesName = fullSpeciesNames
  # )
  bracken_filenames = expand("results/kronaGraph/{sample_id}_-_{sample_name_main}_-_{fullSpeciesName}.kraken.krona.html", zip,
  sample_id = wildcards.sample_id,
  sample_name_main = wildcards.sample_name_main,
  fullSpeciesName = fullSpeciesNames
  )

  return bracken_filenames

# run Kraken 2 to classify
rule kraken2:
    input:
        fastq1 = "data/MMETSP_fastq/{sample_id}_-_{sample_name_main}/{fullSpeciesName}.1.fastq.gz",
        fastq2 = "data/MMETSP_fastq/{sample_id}_-_{sample_name_main}/{fullSpeciesName}.2.fastq.gz"
    output:
        krakenOut = "results/krakenOutputs/{sample_id}_-_{sample_name_main}_-_{fullSpeciesName}.kraken",
        reportOut = "results/krakenReports/{sample_id}_-_{sample_name_main}_-_{fullSpeciesName}.report"
    resources:
      mem_mb=50000
    shell:
        "kraken2 --threads 32 --report-minimizer-data --db data/mmetsp_genome_db --report {output.reportOut} --paired {input.fastq1} {input.fastq2} > {output.krakenOut}"

# Bracken to reestimate species abundance
rule bracken:
  input:
    krakenReport = "results/krakenReports/{sample_id}_-_{sample_name_main}_-_{fullSpeciesName}.report"
  output:
    brackenOut = "results/bracken/{sample_id}_-_{sample_name_main}_-_{fullSpeciesName}.bracken"
  shell:
    "bracken -d data/mmetsp_genome_db -i {input.krakenReport} -o {output.brackenOut} -r 50"

rule krakenKrona:
  input:
    kraken = "results/krakenOutputs/{sample_id}_-_{sample_name_main}_-_{fullSpeciesName}.kraken"
  output:
    krakenKrona = "results/krakenKrona/{sample_id}_-_{sample_name_main}_-_{fullSpeciesName}.kraken.krona"
  shell:
    "cat {input.kraken} | cut -f 3 > {output.krakenKrona}"

rule kronaGraph:
  input:
    krakenKrona = "results/krakenKrona/{sample_id}_-_{sample_name_main}_-_{fullSpeciesName}.kraken.krona"
  output:
    kronaGraph = "results/kronaGraph/{sample_id}_-_{sample_name_main}_-_{fullSpeciesName}.kraken.krona.html"
  shell:
    "ktImportTaxonomy {input.krakenKrona} -o {output.kronaGraph} -t 1"

rule finished:
  input: aggregate_untarFastq
  output: "results/done_krona/{sample_id}_-_{sample_name_main}.done"
  shell: "touch {output}"