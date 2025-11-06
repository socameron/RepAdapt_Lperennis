####################################
#### ESTIMATE DEPTH OF COVERAGE ####
####################################

## ESTIMATE DEPTH USING RepAdapt Scripts (credit Gabriele Nocchi)
# See: https://github.com/RepAdapt/snp_calling_simple_bcftools_slurm/blob/main/08_cnv_2.sh

# The gene annotation is based on haplotype 2, so we use haplotype2
rule prepare_depth_files:
  input:
    fai="data/reference/hap2/lupinehap2.fasta.fai",
    fasta="data/reference/hap2/lupinehap2.fasta",
    gff="data/annotation/MCG3698_Lupinus_perennis.annotation.gff"
  output:
    genome_bed="data/reference/hap2/lupinus_perennis_genome.bed",
    windows_bed="data/reference/hap2/lupinus_perennis_windows.bed",
    windows_list="data/reference/hap2/lupinus_perennis_windows.list",
    genes_bed="data/reference/hap2/lupinus_perennis_genes.bed",
    genes_list="data/reference/hap2/lupinus_perennis_genes.list"
  shell:
    """
    # Create a genome file for bedtools (chromosome and length)
    awk '{{print $1"\\t"$2}}' {input.fai} > {output.genome_bed}
    
    # Create a BED file of 5000 bp windows using the FAI file
    awk -v w=5000 '{{chr = $1; chr_len = $2;
      for (start = 0; start < chr_len; start += w) {{
          end = ((start + w) < chr_len ? (start + w) : chr_len);
          print chr "\\t" start "\\t" end;
      }}
    }}' {input.fai} > {output.windows_bed}
    
    # Create a sorted list of window locations
    awk -F "\\t" '{{print $1":"$2"-"$3}}' {output.windows_bed} | sort -k1,1 > {output.windows_list}
    
    # Create a BED file for each gene using the GFF file
    awk '$3 == "gene" {{print $1"\\t"$4"\\t"$5}}' {input.gff} | uniq > {output.genes_bed}
    
    # Sort the gene BED file based on the order of the reference (from the FAI)
    cut -f1 {input.fai} | while read chr; do 
        awk -v chr="$chr" '$1 == chr {{print $0}}' {output.genes_bed} | sort -k2,2n; 
    done > genes.sorted.bed
    mv genes.sorted.bed {output.genes_bed}
    
    # Create a sorted list of gene locations
    awk -F "\\t" '{{print $1":"$2"-"$3}}' {output.genes_bed} | sort -k1,1 > {output.genes_list}
    """

rule estimate_depth_RepAdapt:
  input:
    fasta="data/reference/hap2/lupinehap2.fasta",
    bam="results/bam_realign/hap2/{sample_prefix}_hap2_realign.bam"
  output:
    temp_depth="results/depths/RepAdapt_temp/{sample_prefix}.depth"
  log:
    "results/logs/create_temp_depth/{sample_prefix}.log"
  envmodules:
    "samtools/1.20"
  shell:
    """
    samtools depth --reference {input.fasta} -aa {input.bam} > {output.temp_depth}
    """

rule estimate_depth_RepAdapt_stats:
  input:
    temp_depth="results/depths/RepAdapt_temp/{sample_prefix}.depth",
    genome_bed="data/reference/hap2/lupinus_perennis_genome.bed",
    windows_bed="data/reference/hap2/lupinus_perennis_windows.bed",
    windows_list="data/reference/hap2/lupinus_perennis_windows.list",
    genes_bed="data/reference/hap2/lupinus_perennis_genes.bed",
    genes_list="data/reference/hap2/lupinus_perennis_genes.list"
  output:
    wg="results/depths/RepAdapt_method/{sample_prefix}-wg.txt",
    genes_sorted="results/depths/RepAdapt_method/{sample_prefix}-genes.sorted.tsv",
    windows_sorted="results/depths/RepAdapt_method/{sample_prefix}-windows.sorted.tsv"
  params:
    temp_window="results/depths/RepAdapt_temp/{sample_prefix}-windows.tsv",
    temp_genes="results/depths/RepAdapt_temp/{sample_prefix}-genes.tsv"
  log:
    "results/logs/estimate_depth_RepAdapt/{sample_prefix}.log"
  envmodules:
    "samtools/1.20",
    "bedtools/2.31.0"
  shell:
    """
    set +o pipefail #we force the rule because we tested it line by line and it works, just somehow not in one rule

    # Gene depth analysis: compute the mean depth per gene
    awk '{{print $1"\\t"$2"\\t"$2"\\t"$3}}' {input.temp_depth} | bedtools map -a {input.genes_bed} -b stdin -c 4 -o mean -null 0 -g {input.genome_bed} | awk -F "\\t" '{{print $1":"$2"-"$3"\\t"$4}}' | sort -k1,1 > {params.temp_genes} || true

    join -a 1 -e 0 -o '1.1 2.2' -t $'\\t' {input.genes_list} {params.temp_genes} > {output.genes_sorted} || true

    # Window depth analysis: compute the mean depth per window
    awk '{{print $1"\\t"$2"\\t"$2"\\t"$3}}' {input.temp_depth} | bedtools map -a {input.windows_bed} -b stdin -c 4 -o mean -null 0 -g {input.genome_bed} | awk -F "\\t" '{{print $1":"$2"-"$3"\\t"$4}}' | sort -k1,1 > {params.temp_window} || true

    join -a 1 -e 0 -o '1.1 2.2' -t $'\\t' {input.windows_list} {params.temp_window} > {output.windows_sorted} || true

    # Overall genome depth (average depth across all positions)
    awk '{{sum += $3; count++}} END {{if (count > 0) print sum/count; else print "No data"}}' {input.temp_depth} > {output.wg}

    # Cleanup temporary files
    rm -f {input.temp_depth} {params.temp_genes} {params.temp_window}
    """


rule combine_depth_RepAdapt:
  input:
    wg=expand("results/depths/RepAdapt_method/{sample_prefix}-wg.txt", sample_prefix=sample_prefixes),
    genes=expand("results/depths/RepAdapt_method/{sample_prefix}-genes.sorted.tsv", sample_prefix=sample_prefixes),
    windows=expand("results/depths/RepAdapt_method/{sample_prefix}-windows.sorted.tsv", sample_prefix=sample_prefixes)
  output:
    combined_windows="results/depths/RepAdapt_method/combined_windows.tsv",
    combined_genes="results/depths/RepAdapt_method/combined_genes.tsv",
    combined_wg="results/depths/RepAdapt_method/combined_wg.tsv"
  params:
    depth_header="results/depths/RepAdapt_temp/depthheader.txt",
    samples_list="results/depths/RepAdapt_temp/samples.txt",
    genes_temp="results/depths/RepAdapt_temp/combined-genes.temp",
    windows_temp="results/depths/RepAdapt_temp/combined-windows.temp",
    # Create a newline-separated string of sample names from your Python variable.
    sample_names="\n".join(sample_prefixes)
  log:
    "results/logs/estimate_depth_RepAdapt/combine_depths.log"
  shell:
    """
    # Ensure temporary directory exists.
    mkdir -p results/depths/RepAdapt_temp

    # Create a file with the sample names using the Python-provided list.
    echo -e "{params.sample_names}" > {params.samples_list}

    # Create a header using the sample names.
    # Replace newlines with tabs for the header.
    echo -e "location\\t$(cat {params.samples_list} | tr '\n' '\\t')" > {params.depth_header}

    # Combine window depth results:
    while read samp; do 
      cut -f2 results/depths/RepAdapt_method/${{samp}}-windows.sorted.tsv > results/depths/RepAdapt_method/${{samp}}-windows.depthcol; 
    done < {params.samples_list}
    first_sample=$(head -n1 {params.samples_list})
    paste results/depths/RepAdapt_method/${{first_sample}}-windows.sorted.tsv $(tail -n +2 {params.samples_list} | sed 's/.*/results\/depths\/RepAdapt_method\/&-windows.depthcol/') > {params.windows_temp}
    cat {params.depth_header} {params.windows_temp} > {output.combined_windows}

    # Combine gene depth results:
    while read samp; do 
      cut -f2 results/depths/RepAdapt_method/${{samp}}-genes.sorted.tsv > results/depths/RepAdapt_method/${{samp}}-genes.depthcol; 
    done < {params.samples_list}
    first_sample=$(head -n1 {params.samples_list})
    paste results/depths/RepAdapt_method/${{first_sample}}-genes.sorted.tsv $(tail -n +2 {params.samples_list} | sed 's/.*/results\/depths\/RepAdapt_method\/&-genes.depthcol/') > {params.genes_temp}
    cat {params.depth_header} {params.genes_temp} > {output.combined_genes}

    # Combine whole-genome depth results:
    while read samp; do 
      echo -e "${{samp}}\t$(cat results/depths/RepAdapt_method/${{samp}}-wg.txt)"; 
    done < {params.samples_list} > {output.combined_wg}
    """


