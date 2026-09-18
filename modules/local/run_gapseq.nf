process RUN_GAPSEQ {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(faa)

    output:
    tuple val(sample_id),
          path("${sample_id}.xml"),
          emit: models

    tuple val(sample_id),
          path("${sample_id}.RDS"),
          emit: rds_models

    tuple val(sample_id),
          path("${sample_id}-all-Reactions.tbl"),
          path("${sample_id}-all-Pathways.tbl"),
          path("${sample_id}-Transporter.tbl"),
          emit: evidence

    tuple val(sample_id),
          path("${sample_id}.gapseq.log"),
          emit: logs

    script:
    def extra_args = params.gapseq_args?.toString() ?: ""

    def database_arg = params.gapseq_db
        ? "-D '${params.gapseq_db}' -O"
        : ""

    def medium_arg = params.gapseq_medium
        ? "-m '${params.gapseq_medium}'"
        : ""

    """
    set -o pipefail

    gapseq doall \
        -A '${params.gapseq_aligner}' \
        -K ${task.cpus} \
        -t '${params.gapseq_taxonomy}' \
        -f . \
        ${database_arg} \
        ${medium_arg} \
        ${extra_args} \
        '${faa}' \
        2>&1 | tee '${sample_id}.gapseq.log'
    """

    stub:
    """
    touch '${sample_id}.xml'
    touch '${sample_id}.RDS'

    touch '${sample_id}-all-Reactions.tbl'
    touch '${sample_id}-all-Pathways.tbl'
    touch '${sample_id}-Transporter.tbl'

    echo 'Stub gapseq run: ${sample_id}' \
        > '${sample_id}.gapseq.log'
    """
}