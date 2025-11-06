##########################
#### DATA PREPARATION ####
##########################

# Trim adapter ends off each sequence file using Trimmomatic
rule trim_reads:
  input:
    r1=lambda wildcards: get_fastq_paths(wildcards.batch, wildcards.sample, "R1"),
    r2=lambda wildcards: get_fastq_paths(wildcards.batch, wildcards.sample, "R2")
  output:
    r1="results/fastq_trimmed/{batch}/{sample}_R1.fastq.gz",
    r1_unp="results/fastq_trimmed/{batch}/{sample}_R1_unpaired.fastq.gz",
    r2="results/fastq_trimmed/{batch}/{sample}_R2.fastq.gz",
    r2_unp="results/fastq_trimmed/{batch}/{sample}_R2_unpaired.fastq.gz"
  log:
    "results/logs/trim_reads/{batch}/{sample}.log"
  envmodules:
    "trimmomatic/0.39"
  params:
    adapters="$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa"
  shell:
    "java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE {input.r1} {input.r2} "
    "{output.r1} {output.r1_unp} {output.r2} {output.r2_unp} "
    "ILLUMINACLIP:{params.adapters}:2:30:10:1 " #previously had ':True' but this failed even though there's documentation of this online
    "LEADING:3 "
    "TRAILING:3 "
    "SLIDINGWINDOW:4:15 " # removes low quality bases
    "MINLEN:36 2> {log}"

# Creating faidx index for reference genome
rule faidx_reference:
  input:
    "data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    "data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  log:
    "results/logs/refgen/lupinehap{hap}_faidx.log"
  envmodules:
    "samtools/1.20"
  shell:
    "samtools faidx {input} 2> {log}"

# Rules for indexing reference genomes (haplotypes 1 and 2)
rule index_reference:
  input:
    "data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    multiext("data/reference/hap{hap}/lupinehap{hap}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
  log:
    "results/logs/refgen/lupinehap{hap}_bwa_index.log"
  envmodules:
    "bwa/0.7.18"
  shell:
    "bwa index {input} 2> {log}"

# Rules for creating dictionary files
rule create_dict:
  input:
    "data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    "data/reference/hap{hap}/lupinehap{hap}.dict"
  log:
    "results/logs/refgen/hap{hap}_dict.log"
  envmodules:
    "samtools/1.20"
  shell:
    "samtools dict {input} > {output} 2> {log}"

# Removing reads smaller than 70bp so 'bwa mem' works better
# Note for a read to be kept, it must be greater than 70bp in BOTH read 1 and 2.
# There are some tools available in envmodules 'seqkit', but I found none that were compatible
rule filter_short_reads_with_seqkit:
  input:
    r1="results/fastq_trimmed/{batch}/{sample}_R1.fastq.gz",
    r2="results/fastq_trimmed/{batch}/{sample}_R2.fastq.gz"
  output:
    r1_filtered="results/fastq_filtered/{batch}/{sample}_R1_filtered.fastq.gz",
    r2_filtered="results/fastq_filtered/{batch}/{sample}_R2_filtered.fastq.gz"
  envmodules:
    "StdEnv/2020",
    "seqkit/2.3.1"
  log:
    "results/logs/filter_fastq_70bp/{batch}/{sample}_filter_fastq.log"
  shell:
    """
    seqkit seq -m 70 {input.r1} | gzip > {output.r1_filtered}
    seqkit seq -m 70 {input.r2} | gzip > {output.r2_filtered}
    """

# Pair filtered reads back together
# Note: naming convention stays the same! Just different folder
rule pair_filtered_reads:
  input:
    r1_filtered="results/fastq_filtered/{batch}/{sample}_R1_filtered.fastq.gz",
    r2_filtered="results/fastq_filtered/{batch}/{sample}_R2_filtered.fastq.gz"
  output:
    r1_paired="results/fastq_paired/{batch}/{sample}_R1_filtered.fastq.gz",
    r2_paired="results/fastq_paired/{batch}/{sample}_R2_filtered.fastq.gz"
  envmodules:
    "StdEnv/2020",
    "seqkit/2.3.1"
  log:
    "results/logs/pair_fastq_70bp/{batch}/{sample}_pair_fastq.log"
  shell:
    """
    seqkit pair \
    -1 {input.r1_filtered} \
    -2 {input.r2_filtered} \
    -O results/fastq_paired/{wildcards.batch}
    """

# NOTE: Prior to mapping, some people like to merge fastqs from the same individual/library and remove PCR duplicates prior to mapping using SuperDeduper from HTStream (last used 2020)
# This might be helpful in reducing heterozygote excess, however I have opted NOT to do this as it may be outdated.

# Mapping/Aligning reads to reference haplotypes
rule map_reads:
  input:
    r1="results/fastq_paired/{batch}/{sample}_R1_filtered.fastq.gz",
    r2="results/fastq_paired/{batch}/{sample}_R2_filtered.fastq.gz",
    genome="data/reference/hap{hap}/lupinehap{hap}.fasta",
    idx=multiext("data/reference/hap{hap}/lupinehap{hap}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
  output:
    "results/bam_raw/hap{hap}/{batch}/{sample}_hap{hap}.bam"
  log:
    "results/logs/map_reads/hap{hap}/{batch}/{sample}_hap{hap}.log"
  envmodules:
    "bwa/0.7.18",
    "samtools/1.20"
  threads: 12 
  params:
    RG="-R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA' "
  shell:
    """
    bwa mem {params.RG} -t {threads} {input.genome} {input.r1} {input.r2} |\
    samtools view -u |\
    samtools sort - > {output}) 2> {log}
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
    # Use the sample name from the wildcard for all read group fields
    rgid=lambda wildcards: wildcards.sample,
    rglb=lambda wildcards: f"{wildcards.sample}_LB",
    rgpl="ILLUMINA",
    rgpu=lambda wildcards: wildcards.batch,
    rgsm=lambda wildcards: wildcards.sample
  log:
    "results/logs/add_read_groups/hap{hap}/{batch}/{sample}_hap{hap}_RG.log"
  envmodules:
    "picard/3.1.0"
  shell:
    """
    java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
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
    "results/bam_raw/hap{hap}/merged/{sample_prefix}_hap{hap}.bam"
  log:
    "results/logs/merge_replicates/hap{hap}/{sample_prefix}_hap{hap}.log"
  envmodules:
    "samtools/1.17"
  shell:
    """
    module load StdEnv/2020
    module load samtools/1.17
    echo "Input files: {input}" >> {log}
    echo "Output file: {output}" >> {log}
    
    # If no input files are found, exit with an error.
    if [ -z "$(echo {input} | tr -d '[:space:]')" ]; then
      echo "No files found for {wildcards.sample_prefix} in hap{wildcards.hap}" >> {log}
      exit 1
    # If there's only one file, just copy it.
    elif [ $(echo {input} | wc -w) -eq 1 ]; then
      echo "Single file found for {wildcards.sample_prefix} in hap{wildcards.hap}, copying to merged folder." >> {log}
      cp {input} {output}
    else
      echo "Multiple files found for {wildcards.sample_prefix} in hap{wildcards.hap}. Merging..." >> {log}
      samtools merge -f -o {output} {input} 2>> {log}
    fi
    sync

    rm -f results/bam_raw/hap{wildcards.hap}/merged/*rep2*.bam
    """


# Marking and removing PCR duplicates + index
# Note that: /bam_mkdup/ is where marked duplicates are marked AND removed.
# marked duplicates here are now removed (updated March 26, 2025)
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
    "gatk/4.4.0.0"
  shell:
    """
    gatk MarkDuplicates \
    --CREATE_INDEX \
    -I {input.raw_bam} \
    -O {output.bam} \
    -M {output.metrics} \
    --REMOVE_DUPLICATES true \
    2> {log}
    """

# Clip overlapping reads
# Previous we loaded packages nixpkgs/16.09 intel/2018.3 to run bamutil/1.0.14 but these might be depreciated with new changes to the Supercomputer
# We use samtools/1.18 because it is compatible with the output of bamutil/1.0.14
rule clip_overlapping_reads:
  input:
    bam="results/bam_mkdup/hap{hap}/{sample_prefix}_hap{hap}_mkdup.bam"
  output:
    clipped_bam="results/bam_clipped/hap{hap}/{sample_prefix}_hap{hap}_clipped.bam",
    clipped_index="results/bam_clipped/hap{hap}/{sample_prefix}_hap{hap}_clipped.bai"
  log:
    "results/logs/clip_overlap/hap{hap}/{sample_prefix}_clip_overlapping_reads.log"
  shell:
    """
    module --force purge
    module load StdEnv/2020
    module load bamutil/1.0.14
    bam clipOverlap --in {input.bam} --out {output.clipped_bam} 2> {log} || true
    module --force purge
    module load StdEnv/2023 samtools/1.18
    samtools index -b {output.clipped_bam} -o {output.clipped_index} --threads 4
    module --force purge
    """

# After downstream analyses (admixture), it appears that LCTGP-19 and SWCP-19 are mislabeled and need to swap names.
# The first and second raw reads (fastq files) for both SWCP-19 and LCTGP-19 still match. It was just mislabelled at the sequencing step!
# NOTE: The metadata in the RG group is incorrect still and hasn't been updated!
rule rename_specific_bam_files:
  input:
    lctgp_bam="results/bam_clipped/hap{hap}/LCTGP-19_hap{hap}_clipped.bam",
    lctgp_bai="results/bam_clipped/hap{hap}/LCTGP-19_hap{hap}_clipped.bai",
    swcp_bam="results/bam_clipped/hap{hap}/SWCP-19_hap{hap}_clipped.bam",
    swcp_bai="results/bam_clipped/hap{hap}/SWCP-19_hap{hap}_clipped.bai"
  output:
    checkpoint="results/checkpoints/hap{hap}/rename_specific_files_checkpoint.txt"
  log:
    "results/logs/rename_specific_files/hap{hap}/rename_specific_files.log"
  shell:
    """
    # Temporary files to avoid overwriting
    tmp_lctgp_bam="results/bam_clipped/hap{wildcards.hap}/LCTGP-19_temp.bam"
    tmp_lctgp_bai="results/bam_clipped/hap{wildcards.hap}/LCTGP-19_temp.bai"
    tmp_swcp_bam="results/bam_clipped/hap{wildcards.hap}/SWCP-19_temp.bam"
    tmp_swcp_bai="results/bam_clipped/hap{wildcards.hap}/SWCP-19_temp.bai"

    # Copy the files to temporary locations
    cp {input.lctgp_bam} $tmp_lctgp_bam
    cp {input.lctgp_bai} $tmp_lctgp_bai
    cp {input.swcp_bam} $tmp_swcp_bam
    cp {input.swcp_bai} $tmp_swcp_bai

    # Rename LCTGP-19 to SWCP-19 and vice versa using the temporary files
    mv $tmp_lctgp_bam results/bam_clipped/hap{wildcards.hap}/SWCP-19_hap{wildcards.hap}_clipped.bam
    mv $tmp_lctgp_bai results/bam_clipped/hap{wildcards.hap}/SWCP-19_hap{wildcards.hap}_clipped.bai
    mv $tmp_swcp_bam results/bam_clipped/hap{wildcards.hap}/LCTGP-19_hap{wildcards.hap}_clipped.bam
    mv $tmp_swcp_bai results/bam_clipped/hap{wildcards.hap}/LCTGP-19_hap{wildcards.hap}_clipped.bai

    # Remove the temporary files
    rm -f $tmp_lctgp_bam $tmp_lctgp_bai $tmp_swcp_bam $tmp_swcp_bai

    # Create a checkpoint file to indicate renaming is complete
    echo "LCTGP-19 and SWCP-19 file names successfully swapped and renamed!" > {output.checkpoint}
    """

# Create bam list per population for entry into realign indels 
# Set for haplotype 2
rule generate_clipped_bam_list_per_population:
  input:
    expand("results/bam_clipped/hap2/{sample_prefix}_hap2_clipped.bam", sample_prefix=sample_prefixes),
    checkpoint="results/checkpoints/hap2/rename_specific_files_checkpoint.txt"
  output:
    "data/lists/hap2/{population}_clipped_hap2.list"
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
    bam_list = "data/lists/hap{hap}/{population}_clipped_hap{hap}.list",
    reference = "data/reference/hap{hap}/lupinehap{hap}.fasta",
  output:
    intervals = "data/lists/hap{hap}/{population}_hap{hap}_indels.intervals",
  log:
    "results/logs/indel_list/hap{hap}/{population}_hap{hap}_indel_list.log",
  threads: 4
  envmodules:
    "apptainer/1.3.5",
  shell:
    """
    set -euo pipefail
    mkdir -p "$(dirname {output.intervals})" "$(dirname {log})"

    # Resolve your host scratch root (e.g., /home/socamero/links/scratch)
    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"

    if [[ ! -d "$HOST_SCRATCH" ]]; then
        echo "Host scratch path not found: $HOST_SCRATCH" >&2
        exit 1
    fi

    # Bind current project dir (for relative paths) AND map host scratch -> /links/scratch inside container
    env -u JAVA_TOOL_OPTIONS \
    apptainer exec --cleanenv \
      -B "$PWD:$PWD","$HOST_SCRATCH:/links/scratch" \
      --pwd "$PWD" \
      /home/socamero/gatk3.sif \
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
    bam = "results/bam_clipped/hap{hap}/{sample_prefix}_hap{hap}_clipped.bam",
    ref = "data/reference/hap{hap}/lupinehap{hap}.fasta",
    intervals = lambda wildcards: "data/lists/hap" + wildcards.hap + "/" + next(pop for pop in POPULATIONS if pop in wildcards.sample_prefix) + f"_hap{wildcards.hap}_indels.intervals",
  output:
    realigned_bam = "results/bam_realign/hap{hap}/{sample_prefix}_hap{hap}_realign.bam",
  log:
    "results/logs/realign_indels/hap{hap}/{sample_prefix}_hap{hap}_realign_indels.log",
  envmodules:
    "apptainer/1.3.5",
  shell:
    """
    set -euo pipefail
    mkdir -p "$(dirname {input.intervals})" "$(dirname {log})"

    # Resolve your host scratch root (e.g., /home/socamero/links/scratch)
    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"

    if [[ ! -d "$HOST_SCRATCH" ]]; then
        echo "Host scratch path not found: $HOST_SCRATCH" >&2
        exit 1
    fi

    # Bind your current project dir (for relative paths) + /links/scratch (for any absolute paths in lists)
    env -u JAVA_TOOL_OPTIONS \
    apptainer exec --cleanenv \
      -B "$PWD:$PWD","$HOST_SCRATCH:/links/scratch" \
      --pwd "$PWD" \
      /home/socamero/gatk3.sif \
      java -Xms2g -Xmx16g -jar /usr/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R {input.ref} \
        -I {input.bam} \
        -targetIntervals {input.intervals} \
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


