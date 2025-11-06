####################################
#### FastQC and MultiQC Reports ####
####################################

# Unzip trimmed files as fastqc 0.12.0 cannot properly read compressed files. 
# Tried with 0.12.1 and it works!
rule unzip_files:
  input:
    zipped_r1="results/fastq_trimmed/{batch}/{sample}_R1.fastq.gz",
    zipped_r2="results/fastq_trimmed/{batch}/{sample}_R2.fastq.gz"
  output: 
    unzipped_r1="results/fastq_trimmed_unzip/{batch}/{sample}_R1.fastq",
    unzipped_r2="results/fastq_trimmed_unzip/{batch}/{sample}_R2.fastq"
  shell:
    "gunzip -c {input.zipped_r1} > {output.unzipped_r1} && gunzip -c {input.zipped_r2} > {output.unzipped_r2}"

# Run FastQC per each trimmed sequence file
# Attempting zipped files since using fastqc/0.12.1
rule fastqc:
  input:
    fastq_r1="results/fastq_trimmed/{batch}/{sample}_R1.fastq.gz",
    fastq_r2="results/fastq_trimmed/{batch}/{sample}_R2.fastq.gz"
  output:
    html_report_r1="results/fastqc/{batch}/{sample}_R1_fastqc.html",
    zip_report_r1="results/fastqc/{batch}/{sample}_R1_fastqc.zip",
    html_report_r2="results/fastqc/{batch}/{sample}_R2_fastqc.html",
    zip_report_r2="results/fastqc/{batch}/{sample}_R2_fastqc.zip"
  log:
    path="results/logs/fastQC/{batch}/{sample}.log"
  envmodules:
    "fastqc/0.12.1"
  shell:
    "fastqc {input.fastq_r1} {input.fastq_r2} --outdir results/fastqc/{wildcards.batch} 2> {log.path}"

# A little different from above because using the raw fastq paths as well as deriving the filenames to get expected fastqc outputs
# Call in rule all with: expected_fastqc_outputs (this is a predetermined list from above)
rule fastqc_raw:
  input:
    fastq_r1=lambda wildcards: get_fastq_paths(wildcards.batch, wildcards.sample, "R1"),
    fastq_r2=lambda wildcards: get_fastq_paths(wildcards.batch, wildcards.sample, "R2")
  output:
    html_report_r1="results/fastqc_raw/{batch}/{sample}_R1_fastqc.html",
    zip_report_r1="results/fastqc_raw/{batch}/{sample}_R1_fastqc.zip",
    html_report_r2="results/fastqc_raw/{batch}/{sample}_R2_fastqc.html",
    zip_report_r2="results/fastqc_raw/{batch}/{sample}_R2_fastqc.zip"
  params:
    # Compute the default names produced by fastqc from the full input filenames.
    r1_default=lambda wildcards: os.path.basename(get_fastq_paths(wc.batch, wc.sample, "R1")).replace(".fastq.gz", ""),
    r2_default=lambda wildcards: os.path.basename(get_fastq_paths(wc.batch, wc.sample, "R2")).replace(".fastq.gz", "")
  log:
    path="results/logs/fastQC_raw/{batch}/{sample}.log"
  envmodules:
    "fastqc/0.12.1"
  shell:
    """
    # Run FastQC: it will produce files named like:
    #   {params.r1_default}_fastqc.html and {params.r1_default}_fastqc.zip
    #   for read 1, and similarly for read 2.
    fastqc {input.fastq_r1} {input.fastq_r2} --outdir results/fastqc_raw/{wildcards.batch} 2> {log.path}
    
    # Rename the outputs to use your config sample name (the shortened version)
    mv results/fastqc_raw/{wildcards.batch}/{params.r1_default}_fastqc.html {output.html_report_r1}
    mv results/fastqc_raw/{wildcards.batch}/{params.r1_default}_fastqc.zip  {output.zip_report_r1}
    mv results/fastqc_raw/{wildcards.batch}/{params.r2_default}_fastqc.html {output.html_report_r2}
    mv results/fastqc_raw/{wildcards.batch}/{params.r2_default}_fastqc.zip  {output.zip_report_r2}
    """

# Create an aggregated FastQC report using MultiQC.
# Note that we create separate MultiQC reports for batch 1 to 3; and .fastqc files must exist prior to calling it on snakemake
rule multiqc_trimmed:
  input:
    fastqc_dir="results/fastqc/{batch}"
  output:
    html_report="results/multiqc/multiqc_report_trimmed_{batch}.html"
  log:
    path="results/logs/multiqc/{batch}.log"
  params:
    fastqc_dir="results/fastqc/{batch}"
  shell:
    "multiqc -n {output.html_report} {input.fastqc_dir} 2> {log.path}"

rule multiqc_raw:
  input:
    fastqc_dir="results/fastqc_raw/{batch}"
  output:
    html_report="results/multiqc_raw/multiqc_report_raw_{batch}.html"
  log:
    path="results/logs/multiqc_raw/{batch}.log"
  params:
    fastqc_dir="results/fastqc_raw/{batch}"
  shell:
    "multiqc -n {output.html_report} {input.fastqc_dir} 2> {log.path}"

