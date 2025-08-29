nextflow.enable.dsl = 2

params.mode = 'p'
params.prefix = null
params.cluster_map = null
params.outdir = 'results'

plink_files = Channel.fromPath("inputs/${params.prefix}.{bed,bim,fam}")
cluster_map = Channel.fromPath(params.cluster_map)

process make_keep_lists{
    input:
    path cluster_map_file

    output:
    path "keeps/*.txt"

    script:
    """
    mkdir -p keeps
    awk '{print \$3}' ${cluster_map_file} | sort -u > clusters.txt
    while read C; do
      awk -v c="\$C" '\$3==c {print \$1, \$2}' ${cluster_map_file} > keeps/\${C}.txt
    done < clusters.txt
    """
}

process plink_subset_per_cluster{
    tag { keepfile.getBaseName() }
    publishDir "${params.outdir}/intermediates", mode: 'copy'

    input:
    path keepfile

    output:
    path "*.bed"
    path "*.bim"
    path "*.fam"

    script:
    """
    mkdir -p results/intermediates

    plink \\
      --bfile ${projectDir}/inputs/${params.prefix} \\
      --keep ${keepfile} \\
      --make-bed \\
      --out ${keepfile.getBaseName()}
    """

}

workflow {
    log.info "Mode: ${params.mode}"
    log.info "Prefix: ${params.prefix}"
    log.info "Cluster map: ${params.cluster_map}"
    log.info "Output directory: ${params.outdir}"
    plink_files.view()
    cluster_map.view()
    keep_ch = make_keep_lists(cluster_map).flatten()
    plink_subset_per_cluster(keep_ch)
}