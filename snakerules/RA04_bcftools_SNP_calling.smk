####################################
#### SNP calling using bcftools ####
####################################

BCFTOOLS_HAP2_UNFILTERED_VCF = "results/bcftools/hap2/all_scaffolds_unfiltered.vcf.gz"
BCFTOOLS_HAP2_FINAL_VCF = "results/bcftools/hap2/rawg0138_Lupinus_perennis_SoSchoenHargreaves.vcf.gz"

# NOTE:  Ignore a ploidy map because L. perennis is diploid and the default setting is diploid
# SNP‐calling rule, one job per scaffold:
rule bcftools_snp_call_by_scaffold:
  input:
    ref     = "data/reference/hap2/lupinehap2.fasta",
    bamlist = "data/lists/hap2/all_realign_hap2.txt"
  output:
    vcf = "results/bcftools/hap2/vcf_by_scaffold/{hap2scaffolds}.vcf.gz",
    tbi = "results/bcftools/hap2/vcf_by_scaffold/{hap2scaffolds}.vcf.gz.tbi"
  log:
    "results/logs/bcftools/SNP_call/{hap2scaffolds}.log"
  params:
    scaffolds = lambda wc: wc.hap2scaffolds,
    container = CONTAINERS["bcftools"]
  threads: 4
  envmodules:
    "apptainer/1.3.5"
  shell:
    r"""
    set -euo pipefail
    mkdir -p "$(dirname {output.vcf})" "$(dirname {log})"

    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"
    BIND="$PWD:$PWD,$HOST_SCRATCH:/links/scratch"

    apptainer exec --cleanenv -B "$BIND" --pwd "$PWD" {params.container} \
      bcftools mpileup -Ou \
        -f {input.ref} \
        --bam-list {input.bamlist} \
        -q 5 \
        -r {params.scaffolds} \
        -I \
        -a FMT/AD,FMT/DP \
    | apptainer exec --cleanenv -B "$BIND" --pwd "$PWD" {params.container} \
      bcftools call \
        -G - \
        -f GQ \
        -mv \
        -Oz \
        -o {output.vcf} \
        --threads {threads} \
      2> {log}

    apptainer exec --cleanenv -B "$BIND" --pwd "$PWD" {params.container} \
      tabix -f -p vcf {output.vcf} \
      2>> {log}
    """

rule concat_scaffold_vcfs:
  input:
    ref_list = "results/scaffolds/hap2_all_scaffolds.txt",
    vcfs = expand("results/bcftools/hap2/vcf_by_scaffold/{hap2scaffolds}.vcf.gz", hap2scaffolds=HAP2SCAFFOLDS),
    tbi  = expand("results/bcftools/hap2/vcf_by_scaffold/{hap2scaffolds}.vcf.gz.tbi", hap2scaffolds=HAP2SCAFFOLDS)
  output:
    vcf = BCFTOOLS_HAP2_UNFILTERED_VCF,
    tbi = BCFTOOLS_HAP2_UNFILTERED_VCF + ".tbi"
  threads: 8
  envmodules:
    "apptainer/1.3.5"
  params:
    container = CONTAINERS["bcftools"]
  log:
    "results/logs/bcftools/concat/all_scaffolds_concat.log"
  shell:
    r"""
    set -euo pipefail
    mkdir -p "$(dirname {output.vcf})" "$(dirname {log})"

    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"
    BIND="$PWD:$PWD,$HOST_SCRATCH:/links/scratch"

    FILELIST="$(dirname {log})/vcf_filelist.txt"

    sed -e 's|^|results/bcftools/hap2/vcf_by_scaffold/|' \
        -e 's|$|.vcf.gz|' \
        {input.ref_list} > "$FILELIST"

    while IFS= read -r f; do
      if [[ ! -s "$f" ]]; then
        echo "Missing or empty VCF: $f" >&2
        exit 1
      fi
    done < "$FILELIST"

    apptainer exec --cleanenv -B "$BIND" --pwd "$PWD" {params.container} \
      bcftools concat \
        --threads {threads} \
        -f "$FILELIST" \
        -Oz \
        -o {output.vcf} \
      > {log} 2>&1

    apptainer exec --cleanenv -B "$BIND" --pwd "$PWD" {params.container} \
      tabix -f -p vcf {output.vcf} \
      2>> {log}
    """


rule filter_bcftools_vcf:
  input:
    vcf = BCFTOOLS_HAP2_UNFILTERED_VCF,
    tbi = BCFTOOLS_HAP2_UNFILTERED_VCF + ".tbi"
  output:
    vcf = BCFTOOLS_HAP2_FINAL_VCF,
    tbi = BCFTOOLS_HAP2_FINAL_VCF + ".tbi"
  log:
    "results/logs/bcftools/filter/all_scaffolds_filtered.log"
  threads: 2
  envmodules:
    "apptainer/1.3.5"
  params:
    container = CONTAINERS["bcftools"]
  shell:
    r"""
    set -euo pipefail
    mkdir -p "$(dirname {output.vcf})" "$(dirname {log})"

    HOST_SCRATCH="$(readlink -f /home/socamero/links/scratch)"
    BIND="$PWD:$PWD,$HOST_SCRATCH:/links/scratch"

    apptainer exec --cleanenv -B "$BIND" --pwd "$PWD" {params.container} \
      bcftools filter \
        -e 'AC=AN || MQ < 30' \
        -Oz \
        -o {output.vcf} \
        {input.vcf} \
      > {log} 2>&1

    apptainer exec --cleanenv -B "$BIND" --pwd "$PWD" {params.container} \
      tabix -f -p vcf {output.vcf} \
      2>> {log}
    """
