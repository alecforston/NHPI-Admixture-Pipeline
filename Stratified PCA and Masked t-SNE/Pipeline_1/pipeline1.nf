nextflow.enable.dsl = 2

params.prefix         = params.prefix         ?: 'inputs/prefix'
params.clusters       = params.clusters       ?: 'inputs/clusters.txt'
params.regroup_script = params.regroup_script ?: 'regroup.py'
params.outdir         = params.outdir         ?: 'results'

workflow {
    Channel.fromPath(params.clusters).set { CLUSTER_MAP }
    Channel.fromPath(params.regroup_script).set { REGROUP_SCRIPT }

    BASE_FILES = Channel.of([
        file("${params.prefix}.bed"),
        file("${params.prefix}.bim"),
        file("${params.prefix}.fam")
    ])

    MAKE_KEEP_FILES(CLUSTER_MAP)

    MAKE_KEEP_FILES.out.keep_files
        .flatten()
        .combine(BASE_FILES)
        .map { keep, bed, bim, fam ->
            tuple(keep.baseName.replace('.txt',''), keep, bed, bim, fam)
        }
        | PLINK_SUBSET

    PLINK_SUBSET.out | GCTA_PCA
    GCTA_PCA.out.map { cluster, eigenvec -> eigenvec }
              .collect()
              .set { eigenvecs_list }

    REGROUP(eigenvecs_list, REGROUP_SCRIPT)
}

process MAKE_KEEP_FILES {
    label 'keep'
    input:
    path cluster_map
    output:
    path "keep/*.txt", emit: keep_files
    script:
    """
    mkdir -p keep
    awk '{ print \$1, \$2 >> ("keep/" \$3 ".txt") }' ${cluster_map}
    """
}

process PLINK_SUBSET {
    label 'plink'
    publishDir "${params.outdir}/plink", mode: 'copy'
    input:
    tuple val(cluster), path(keep), path(bed), path(bim), path(fam)
    output:
    tuple val(cluster), path("${cluster}.bed"), path("${cluster}.bim"), path("${cluster}.fam")
    script:
    """
    plink --bed ${bed} --bim ${bim} --fam ${fam} \
          --keep ${keep} --make-bed --out ${cluster}
    """
}

process GCTA_PCA {
    label 'gcta'
    publishDir "${params.outdir}/gcta", mode: 'copy'
    input:
    tuple val(cluster), path(bed), path(bim), path(fam)
    output:
    tuple val(cluster), path("${cluster}.eigenvec")
    script:
    """
    gcta --bfile ${cluster} --make-grm --out ${cluster}
    gcta --grm ${cluster} --pca 10 --out ${cluster}
    """
}

process REGROUP {
    label 'regroup'
    publishDir "${params.outdir}/regroup", mode: 'copy'
    input:
    path eigenvecs
    path regroup_script
    output:
    path "regroup.log"
    script:
    """
    ARGS=""
    for f in ${eigenvecs}; do
      base=\$(basename "\$f" .eigenvec)
      ARGS="\$ARGS \$base \$f"
    done
    python3 ${regroup_script} \$ARGS | tee regroup.log
    """
}