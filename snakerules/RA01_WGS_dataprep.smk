##########################
#### DATA PREPARATION ####
##########################

# Trim adapter ends off each sequence file using fastp
rule trim_reads:
  input:
    r1=lambda wildcards: get_fastq_paths(wildcards.batch, wildcards.sample, "R1"),
    r2=lambda wildcards: get_fastq_paths(wildcards.batch, wildcards.sample, "R2")
  output:
    r1="results/fastq_trimmed/{batch}/{sample}_R1_trimmed.fastq.gz",
    r2="results/fastq_trimmed/{batch}/{sample}_R2_trimmed.fastq.gz"
  log:
    "results/logs/trim_reads/{batch}/{sample}.log"
  threads: 4
  envmodules:
    "apptainer/1.3.5"
  params:
    container="/home/socamero/containers/fastp_0.20.1.sif"
  shell:
    """
    set -euo pipefail

    mkdir -p "$(dirname {output.r1})" "$(dirname {log})" /home/socamero/containers

    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"

    if [[ ! -d "$HOST_SCRATCH" ]]; then
        echo "Host scratch path not found: $HOST_SCRATCH" >&2
        exit 1
    fi

    if [[ ! -f "{params.container}" ]]; then
        apptainer pull "{params.container}" docker://quay.io/biocontainers/fastp:0.20.1--h8b12597_0
    fi

    apptainer exec --cleanenv \
      -B "$PWD:$PWD","$HOST_SCRATCH:/links/scratch" \
      --pwd "$PWD" \
      "{params.container}" \
      fastp \
        -w {threads} \
        -i {input.r1} \
        -I {input.r2} \
        -o {output.r1} \
        -O {output.r2} \
      > {log} 2>&1
    """

# Creating faidx index for reference genome
rule faidx_reference:
  input:
    "data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    "data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  log:
    "results/logs/refgen/lupinehap{hap}_faidx.log"
  envmodules:
    "apptainer/1.3.5"
  params:
    container=CONTAINERS["samtools"]
  shell:
    """
    set -euo pipefail
    mkdir -p "$(dirname {log})"

    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"

    apptainer exec --cleanenv \
      -B "$PWD:$PWD","$HOST_SCRATCH:/links/scratch" \
      --pwd "$PWD" \
      {params.container} \
      samtools faidx {input} \
      2> {log}
    """

# Rules for indexing reference genomes (haplotypes 1 and 2)
rule index_reference:
  input:
    "data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    multiext("data/reference/hap{hap}/lupinehap{hap}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
  log:
    "results/logs/refgen/lupinehap{hap}_bwa_index.log"
  envmodules:
    "apptainer/1.3.5"
  params:
    container=CONTAINERS["bwa"]
  shell:
    """
    set -euo pipefail
    mkdir -p "$(dirname {log})"

    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"

    apptainer exec --cleanenv \
      -B "$PWD:$PWD","$HOST_SCRATCH:/links/scratch" \
      --pwd "$PWD" \
      {params.container} \
      bwa index -a bwtsw {input} \
      2> {log}
    """

# Rules for creating dictionary files
rule create_dict:
  input:
    "data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    "data/reference/hap{hap}/lupinehap{hap}.dict"
  log:
    "results/logs/refgen/hap{hap}_dict.log"
  envmodules:
    "apptainer/1.3.5"
  params:
    container=CONTAINERS["picard"]
  shell:
    """
    set -euo pipefail
    mkdir -p "$(dirname {log})"

    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"

    apptainer exec --cleanenv \
      -B "$PWD:$PWD","$HOST_SCRATCH:/links/scratch" \
      --pwd "$PWD" \
      {params.container} \
      picard CreateSequenceDictionary \
        R={input} \
        O={output} \
      2> {log}
    """


# NOTE: Prior to mapping, some people like to merge fastqs from the same individual/library and remove PCR duplicates prior to mapping using SuperDeduper from HTStream (last used 2020)
# This might be helpful in reducing heterozygote excess, however I have opted NOT to do this as it may be outdated.

# bwa version '0.7.17' for RepAdapt
# Mapping/Aligning reads to reference haplotypes
rule map_reads:
  input:
    r1="results/fastq_trimmed/{batch}/{sample}_R1_trimmed.fastq.gz",
    r2="results/fastq_trimmed/{batch}/{sample}_R2_trimmed.fastq.gz",
    genome="data/reference/hap{hap}/lupinehap{hap}.fasta",
    idx=multiext("data/reference/hap{hap}/lupinehap{hap}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
  output:
    bam="results/bam_raw/hap{hap}/{batch}/{sample}_hap{hap}.bam",
    bai="results/bam_raw/hap{hap}/{batch}/{sample}_hap{hap}.bam.bai"
  log:
    "results/logs/map_reads/hap{hap}/{batch}/{sample}_hap{hap}.log"
  threads: 4
  envmodules:
    "apptainer/1.3.5"
  params:
    bwa_container=CONTAINERS["bwa"],
    samtools_container=CONTAINERS["samtools"]
  shell:
    """
    set -euo pipefail
    mkdir -p "$(dirname {output.bam})" "$(dirname {log})"

    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"
    BIND="$PWD:$PWD,$HOST_SCRATCH:/links/scratch"

    apptainer exec --cleanenv -B "$BIND" --pwd "$PWD" {params.bwa_container} \
      bwa mem -t {threads} {input.genome} {input.r1} {input.r2} \
      2> {log} \
    | apptainer exec --cleanenv -B "$BIND" --pwd "$PWD" {params.samtools_container} \
      samtools view -Sb -q 10 - \
    | apptainer exec --cleanenv -B "$BIND" --pwd "$PWD" {params.samtools_container} \
      samtools sort --threads {threads} -o {output.bam} -

    apptainer exec --cleanenv -B "$BIND" --pwd "$PWD" {params.samtools_container} \
      samtools index {output.bam} {output.bai} \
      2>> {log}
    """

# We add metadata associated with each .bam file since there are differences in sequencing platforms, etc. 
# This is implemented for RepAdapt data but NOT for past L. perennis analyses
# Batch_1 and Batch_2 were created from the same libraries, so no information is lost when merging!
rule add_read_groups:
  input:
    bam="results/bam_raw/hap{hap}/{batch}/{sample}_hap{hap}.bam"
  output:
    bam="results/bam_raw/hap{hap}/{batch}_RG/{sample}_hap{hap}_RG.bam"
  params:
    rgid=lambda wildcards: wildcards.sample,
    rglb=lambda wildcards: f"{wildcards.sample}_LB",
    rgpl="ILLUMINA",
    rgpu=lambda wildcards: wildcards.batch,
    rgsm=lambda wildcards: wildcards.sample,
    container=CONTAINERS["picard"]
  log:
    "results/logs/add_read_groups/hap{hap}/{batch}/{sample}_hap{hap}_RG.log"
  envmodules:
    "apptainer/1.3.5"
  shell:
    """
    set -euo pipefail
    mkdir -p "$(dirname {output.bam})" "$(dirname {log})"

    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"

    apptainer exec --cleanenv \
      -B "$PWD:$PWD","$HOST_SCRATCH:/links/scratch" \
      --pwd "$PWD" \
      {params.container} \
      picard AddOrReplaceReadGroups \
        I={input.bam} \
        O={output.bam} \
        RGID={params.rgid} \
        RGLB={params.rglb} \
        RGPL={params.rgpl} \
        RGPU={params.rgpu} \
        RGSM={params.rgsm} \
      2> {log}
    """
# ID = unique ID | LB = library | PL = platform | PU = platform unit | SM = sample


# To call 'raw' .bam files, we use the get_fastq_paths function in previous rules because the naming conventions are different between batch_1/batch_2 and batch_3. 
# In rule all, you can call:
    #[f"results/bam_raw/hap2/{batch}_RG/{sample}_hap2_RG.bam"
     #for batch in get_batches()
     #for sample in list(get_samples(batch).keys())],

##############################
#### DATA QUALITY CONTROL ####
##############################

# Steps in QC :
  # 1. Remove PCR and optical duplicates
  # 2. Identify paralogous regions causing mapping problems
      # - Requires ngsParalog and 1st run of ANGSD
      # - ANGSD: SNP call and input SNPs (aka polymorphic sites) into ngsParalog
      # - ngsParalog: probablistic call of paralogous regions
      # - After ngsParalog, calculate p-values based on chi-sq df = 1 with Benjamini-Hochberg correction
      # NOTE: This script is NOT shown, but can be provided. It was previously run on another analysis and so the output SNPs are provided.


## STEP 1: MERGE REPLICATES & REMOVE PCR & OPTICAL DUPLICATES 

# merge replicates if found, otherwise rename and move to hap{hap}/merged/
# NEW: added rm results/bam_raw/hap{hap}/merged/*rep2*.bam because the file remained in the folder. 
# NOTE: merge_replicates must be called before any downstream rules!
rule merge_replicates:
  input:
    lambda wildcards: find_replicates(wildcards.sample_prefix, wildcards.hap)
  output:
    "results/bam_raw/hap{hap}/merged/{sample_prefix}_hap{hap}_RG.bam"
  log:
    "results/logs/merge_replicates/hap{hap}/{sample_prefix}_hap{hap}.log"
  envmodules:
    "apptainer/1.3.5"
  params:
    container=CONTAINERS["samtools"]
  shell:
    """
    set -euo pipefail
    mkdir -p "$(dirname {output})" "$(dirname {log})"

    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"
    BIND="$PWD:$PWD,$HOST_SCRATCH:/links/scratch"

    echo "Input files: {input}" > {log}
    echo "Output file: {output}" >> {log}

    if [ -z "$(echo {input} | tr -d '[:space:]')" ]; then
      echo "No files found for {wildcards.sample_prefix} in hap{wildcards.hap}" >> {log}
      exit 1
    elif [ $(echo {input} | wc -w) -eq 1 ]; then
      echo "Single file found; copying." >> {log}
      cp {input} {output}
    else
      echo "Multiple files found; merging." >> {log}
      apptainer exec --cleanenv \
        -B "$BIND" \
        --pwd "$PWD" \
        {params.container} \
        samtools merge -f -o {output} {input} \
        2>> {log}
    fi

    rm -f results/bam_raw/hap{wildcards.hap}/merged/*rep2*.bam
    """


# Marking and removing PCR duplicates + index
# Use PICARD tools v2.26.3
rule mark_remove_duplicates:
  input:
    raw_bam="results/bam_raw/hap{hap}/merged/{sample_prefix}_hap{hap}_RG.bam"
  output:
    bam="results/bam_mkdup/hap{hap}/{sample_prefix}_hap{hap}_mkdup.bam",
    bai="results/bam_mkdup/hap{hap}/{sample_prefix}_hap{hap}_mkdup.bai",
    metrics="results/qc/mkdup_metrics/{sample_prefix}_hap{hap}.metrics"
  log:
    "results/logs/mark_remove_duplicates/hap{hap}/{sample_prefix}_hap{hap}.log"
  envmodules:
    "apptainer/1.3.5"
  params:
    container=CONTAINERS["picard"]
  shell:
    """
    set -euo pipefail
    mkdir -p "$(dirname {output.bam})" "$(dirname {output.metrics})" "$(dirname {log})"

    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"

    apptainer exec --cleanenv \
      -B "$PWD:$PWD","$HOST_SCRATCH:/links/scratch" \
      --pwd "$PWD" \
      {params.container} \
      picard MarkDuplicates \
        INPUT={input.raw_bam} \
        OUTPUT={output.bam} \
        METRICS_FILE={output.metrics} \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true \
      2> {log}

    if [[ -f "{output.bam}.bai" && ! -f "{output.bai}" ]]; then
      mv "{output.bam}.bai" "{output.bai}"
    fi
    """


# After downstream analyses (admixture), it appears that LCTGP-19 and SWCP-19 are mislabeled and need to swap names.
# The first and second raw reads (fastq files) for both SWCP-19 and LCTGP-19 still match. It was just mislabelled at the sequencing step!
# NOTE: The metadata in the RG group is incorrect still and hasn't been updated!
rule rename_specific_bam_files:
  input:
    lctgp_bam="results/bam_mkdup/hap{hap}/LCTGP-19_hap{hap}_mkdup.bam",
    lctgp_bai="results/bam_mkdup/hap{hap}/LCTGP-19_hap{hap}_mkdup.bai",
    swcp_bam="results/bam_mkdup/hap{hap}/SWCP-19_hap{hap}_mkdup.bam",
    swcp_bai="results/bam_mkdup/hap{hap}/SWCP-19_hap{hap}_mkdup.bai"
  output:
    checkpoint="results/checkpoints/hap{hap}/rename_specific_files_checkpoint.txt"
  log:
    "results/logs/rename_specific_files/hap{hap}/rename_specific_files.log"
  shell:
    """
    # Temporary files to avoid overwriting
    tmp_lctgp_bam="results/bam_mkdup/hap{wildcards.hap}/LCTGP-19_temp.bam"
    tmp_lctgp_bai="results/bam_mkdup/hap{wildcards.hap}/LCTGP-19_temp.bai"
    tmp_swcp_bam="results/bam_mkdup/hap{wildcards.hap}/SWCP-19_temp.bam"
    tmp_swcp_bai="results/bam_mkdup/hap{wildcards.hap}/SWCP-19_temp.bai"

    # Copy the files to temporary locations
    cp {input.lctgp_bam} $tmp_lctgp_bam
    cp {input.lctgp_bai} $tmp_lctgp_bai
    cp {input.swcp_bam} $tmp_swcp_bam
    cp {input.swcp_bai} $tmp_swcp_bai

    # Rename LCTGP-19 to SWCP-19 and vice versa using the temporary files
    mv $tmp_lctgp_bam results/bam_mkdup/hap{wildcards.hap}/SWCP-19_hap{wildcards.hap}_mkdup.bam
    mv $tmp_lctgp_bai results/bam_mkdup/hap{wildcards.hap}/SWCP-19_hap{wildcards.hap}_mkdup.bai
    mv $tmp_swcp_bam results/bam_mkdup/hap{wildcards.hap}/LCTGP-19_hap{wildcards.hap}_mkdup.bam
    mv $tmp_swcp_bai results/bam_mkdup/hap{wildcards.hap}/LCTGP-19_hap{wildcards.hap}_mkdup.bai

    # Remove the temporary files
    rm -f $tmp_lctgp_bam $tmp_lctgp_bai $tmp_swcp_bam $tmp_swcp_bai

    # Create a checkpoint file to indicate renaming is complete
    echo "LCTGP-19 and SWCP-19 file names successfully swapped and renamed!" > {output.checkpoint}
    """

# Create bam list per population for entry into realign indels 
# Set for haplotype 2
rule generate_mkdup_bam_list_per_population:
  input:
    expand("results/bam_mkdup/hap2/{sample_prefix}_hap2_mkdup.bam", sample_prefix=sample_prefixes),
    checkpoint="results/checkpoints/hap2/rename_specific_files_checkpoint.txt"
  output:
    "data/lists/hap2/{population}_mkdup_hap2.list"
  wildcard_constraints:
    population="|".join(POPULATIONS)
  run:
    bam_files = input
    output_file = output[0]
    population = wildcards.population

    with open(output_file, "w") as output:
        for bam_file in bam_files:
            if population in bam_file and f"_hap2_" in bam_file:
                output.write(f"{bam_file}\n")

# Create interval list of indels
# NOTE: We turn the badmate filter OFF! (-drf BadMate) because our genome assembly is scaffold-based and so there may be reads mapped on different areas than their mate. 
# NOTE: We establish targets by a per-sample basis, rather than identifying realignment targets by population. This differs from the original approach!
# Create interval list of indels
# We use apptainer, since I've had issues running gatk/3.8 with simply module load on the DRAC clusters
# MUST PRE-INSTALL gatk/3.8 container in login node
# apptainer pull gatk3.sif docker://broadinstitute/gatk3:3.8-1
rule indel_list:
  input:
    bam_list = "data/lists/hap{hap}/{population}_mkdup_hap{hap}.list",
    reference = "data/reference/hap{hap}/lupinehap{hap}.fasta",
  output:
    intervals = "data/lists/hap{hap}/{population}_hap{hap}_indels.intervals",
  log:
    "results/logs/indel_list/hap{hap}/{population}_hap{hap}_indel_list.log"
  threads: 4
  envmodules:
    "apptainer/1.3.5"
  params:
    container=CONTAINERS["gatk3"]
  shell:
    """
    set -euo pipefail
    mkdir -p "$(dirname {output.intervals})" "$(dirname {log})"

    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"

    env -u JAVA_TOOL_OPTIONS \
    apptainer exec --cleanenv \
      -B "$PWD:$PWD","$HOST_SCRATCH:/links/scratch" \
      --pwd "$PWD" \
      {params.container} \
      java -Xms2g -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R {input.reference} \
        -I {input.bam_list} \
        -o {output.intervals} \
        -drf BadMate \
        -nt {threads} \
      2> {log}
    """

rule realign_indels:
  input:
    bam = "results/bam_mkdup/hap{hap}/{sample_prefix}_hap{hap}_mkdup.bam",
    ref = "data/reference/hap{hap}/lupinehap{hap}.fasta",
    intervals = lambda wildcards: "data/lists/hap" + wildcards.hap + "/" + next(pop for pop in POPULATIONS if pop in wildcards.sample_prefix) + f"_hap{wildcards.hap}_indels.intervals",
  output:
    realigned_bam = "results/bam_realign/hap{hap}/{sample_prefix}_hap{hap}_realign.bam"
  log:
    "results/logs/realign_indels/hap{hap}/{sample_prefix}_hap{hap}_realign_indels.log"
  envmodules:
    "apptainer/1.3.5"
  params:
    container=CONTAINERS["gatk3"]
  shell:
    """
    set -euo pipefail
    mkdir -p "$(dirname {output.realigned_bam})" "$(dirname {log})"

    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"

    env -u JAVA_TOOL_OPTIONS \
    apptainer exec --cleanenv \
      -B "$PWD:$PWD","$HOST_SCRATCH:/links/scratch" \
      --pwd "$PWD" \
      {params.container} \
      java -Xms2g -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R {input.ref} \
        -I {input.bam} \
        -targetIntervals {input.intervals} \
        --consensusDeterminationModel USE_READS \
        -o {output.realigned_bam} \
        -drf BadMate \
      &> {log}
    """


# Rule to create .txt file of BAM files for each haplotype
# Only applies to hap2
rule generate_bam_list_per_haplotype:
  input:
    expand("results/bam_realign/hap2/{sample_prefix}_hap2_realign.bam", sample_prefix=sample_prefixes)
  output:
    "data/lists/hap2/all_realign_hap2.txt"
  run:
    bam_files = input
    output_file = output[0]

    with open(output_file, "w") as output:
        for bam_file in bam_files:
            if f"_hap2_" in bam_file:
                output.write(f"{bam_file}\n")

# Merge bams into batches 1+2 and 3 because of different sequencing platforms
# Only applies to hap2
rule generate_bam_list_per_round:
  input:
    expand("results/bam_realign/hap2/{sample_prefix}_hap2_realign.bam", sample_prefix=sample_prefixes)
  output:
    first_round="data/lists/hap2/first_round_hap2.txt",
    second_round="data/lists/hap2/second_round_hap2.txt"
  run:
    first_round_bams = []
    second_round_bams = []

    # Define sequencing rounds
    first_round_populations = {"HPW", "IDNP-MW", "LCTGP", "MFNP", "PPP", "RLPLV", "SWCP"}
    
    # Extract population name dynamically
    def extract_population(sample_prefix):
        for pop in POPULATIONS:
            if sample_prefix.startswith(pop):
                return pop  # Return the correct population name
        return None  # If no match is found, which should not happen

    for bam_file in input:
        sample_prefix = bam_file.split("/")[-1].split("_hap")[0]  # Extract sample prefix
        population = extract_population(sample_prefix)  # Get correct population

        if population in first_round_populations:
            first_round_bams.append(bam_file)
        else:
            second_round_bams.append(bam_file)

    # Write BAM lists for each sequencing round
    with open(output.first_round, "w") as f:
        for bam in first_round_bams:
            f.write(f"{bam}\n")

    with open(output.second_round, "w") as f:
        for bam in second_round_bams:
            f.write(f"{bam}\n")


