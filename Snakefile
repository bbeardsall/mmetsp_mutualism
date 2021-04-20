import pandas as pd

genome_data = pd.read_csv("data/mmetsp_ncbi_genome_info.csv")
#/iplant/home/shared/imicrobe/projects/104/transcriptomes/MMETSP005/MMETSP0053.fastq.tar
#fastq_files = genome_data['fastq_file']

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
    #expand("results/krakenOutputs/{sample_id}_{sample_name_main}_{fullSpeciesName}.kraken", zip, sample = _sample_ids, sample_name_main = _sample_names)
    expand("results/done_bracken/{sample_id}_-_{sample_name_main}.done", zip, sample_id = _sample_ids, sample_name_main = _sample_names)
    #aggregate_untarFastq

### Preparing database ###

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
  bracken_filenames = expand("results/bracken/{sample_id}_-_{sample_name_main}_-_{fullSpeciesName}.bracken", zip,
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
        "kraken2 --use-names --threads 32 --db data/mmetsp_genome_db --report {output.reportOut} --paired {input.fastq1} {input.fastq2} > {output.krakenOut}"

# Bracken to reestimate species abundance
rule bracken:
  input:
    krakenReport = "results/krakenReports/{sample_id}_-_{sample_name_main}_-_{fullSpeciesName}.report"
  output:
    brackenOut = "results/bracken/{sample_id}_-_{sample_name_main}_-_{fullSpeciesName}.bracken"
  shell:
    "bracken -d data/mmetsp_genome_db -i {input.krakenReport} -o {output.brackenOut} -r 50"

rule finished:
  input: aggregate_untarFastq
  output: "results/done_bracken/{sample_id}_-_{sample_name_main}.done"
  shell: "touch {output}"

