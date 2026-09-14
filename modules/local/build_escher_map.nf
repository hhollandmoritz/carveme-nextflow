process BUILD_ESCHER_MAP {

    tag "${sample_id}:${medium}"

    input:
    tuple val(sample_id), val(medium), path(model)
    path escher_map

    output:
    tuple val(sample_id),
          val(medium),
          path("${sample_id}.${medium}.escher.html"),
          emit: html

    path "${sample_id}.${medium}.escher.log",
         emit: logs

    script:
    def prefix = "${sample_id}.${medium}"

    """
    set -o pipefail

    build_escher_map.py \
        --model '${model}' \
        --map '${escher_map}' \
        --output '${prefix}.escher.html' \
        --reaction-scale-preset '${params.escher_reaction_scale_preset}' \
        2>&1 | tee '${prefix}.escher.log'
    """

    stub:
    def prefix = "${sample_id}.${medium}"

    """
    touch '${prefix}.escher.html'
    touch '${prefix}.escher.log'
    """
}