####################################
#### SNP calling using bcftools ####
####################################

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
    scaffolds = lambda wc: wc.hap2scaffolds
  threads: 4
  envmodules:
    "bcftools/1.22"
  shell:
    r"""
    bcftools mpileup -Ou \
      -f {input.ref} \
      -b {input.bamlist} \
      -q 5 -r {params.scaffolds} -I \
      -a FMT/AD,FMT/DP \
    | bcftools call \
      -G - -f GQ -mv -Oz -o {output.vcf} \
      --threads {threads} 2> {log}
    tabix -f {output.vcf}
    """

rule concat_scaffold_vcfs:
  input:
    ref_list = "results/scaffolds/hap2_all_scaffolds.txt",
    vcfs = expand("results/bcftools/hap2/vcf_by_scaffold/{hap2scaffolds}.vcf.gz", hap2scaffolds=HAP2SCAFFOLDS),
    tbi  = expand("results/bcftools/hap2/vcf_by_scaffold/{hap2scaffolds}.vcf.gz.tbi", hap2scaffolds=HAP2SCAFFOLDS)
  output:
    vcf = "results/bcftools/hap2/all_scaffolds.vcf.gz",
    tbi = "results/bcftools/hap2/all_scaffolds.vcf.gz.tbi"
  threads: 8
  envmodules:
    "bcftools/1.22"
  log:
    "results/logs/bcftools/concat/all_scaffolds_concat.log"
  shell:
    r"""
    set -euo pipefail
    mkdir -p "$(dirname {output.vcf})" "$(dirname {log})"

    # Build filelist from scaffold names (no shell vars needed)
    FILELIST="$(dirname {log})/vcf_filelist.txt"
    sed -e 's|^|results/bcftools/hap2/vcf_by_scaffold/|' \
        -e 's|$|.vcf.gz|' \
        {input.ref_list} > "$FILELIST"

    # sanity check
    while IFS= read -r f; do
      if [[ ! -s "$f" ]]; then
        echo "Missing or empty VCF: $f" >&2
        exit 1
      fi
    done < "$FILELIST"

    bcftools concat --threads {threads} -f "$FILELIST" -Oz -o {output.vcf} > {log} 2>&1
    tabix -f {output.vcf} 2>> {log}
    """