process RUN_CARVEME_DEFAULT {

    tag "${sample_id}:${medium}"

    input:
    tuple val(sample_id),
          path(faa),
          val(medium)

    output:
    tuple val(sample_id),
          val(medium),
          path("${sample_id}.${medium}.sbml"),
          emit: models

    tuple val(sample_id),
          val(medium),
          path("${sample_id}.${medium}.log"),
          emit: logs

    script:
    def model = "${sample_id}.${medium}.sbml"
    def log   = "${sample_id}.${medium}.log"

    """
    set -o pipefail

    carve '${faa}' \
        --output '${model}' \
        --gapfill '${medium}' \
        --init '${medium}' \
        --solver '${params.solver}' \
        --verbose \
        2>&1 | tee '${log}'
    """
    stub:
    """
    touch '${sample_id}.${medium}.sbml'
    echo 'Stub built-in medium run: ${sample_id} ${medium}' \
        > '${sample_id}.${medium}.log'
    """
}

process RUN_CARVEME_CUSTOM {

    tag "${sample_id}:${medium}"

    input:
    tuple val(sample_id),
          path(faa),
          val(medium),
          path(mediadb)

    output:
    tuple val(sample_id),
          val(medium),
          path("${sample_id}.${medium}.sbml"),
          emit: models

    tuple val(sample_id),
          val(medium),
          path("${sample_id}.${medium}.log"),
          emit: logs

    script:
    def model = "${sample_id}.${medium}.sbml"
    def log   = "${sample_id}.${medium}.log"

    """
    set -o pipefail

    carve '${faa}' \
        --output '${model}' \
        --gapfill '${medium}' \
        --init '${medium}' \
        --mediadb '${mediadb}' \
        --solver '${params.solver}' \
        --verbose \
        2>&1 | tee '${log}'
    """

    stub:
    """
    touch '${sample_id}.${medium}.sbml'
    echo 'Stub custom medium run: ${sample_id} ${medium} ${mediadb}' \
        > '${sample_id}.${medium}.log'
    """
}

process RUN_CARVEME_NO_GAPFILL {

    tag "${sample_id}:no-gapfill"

    input:
    tuple val(sample_id), path(faa)

    output:
    tuple val(sample_id),
          val('NO_GAPFILL'),
          path("${sample_id}.no-gapfill.sbml"),
          emit: models

    tuple val(sample_id),
          val('NO_GAPFILL'),
          path("${sample_id}.no-gapfill.log"),
          emit: logs

    script:
    def model = "${sample_id}.no-gapfill.sbml"
    def log   = "${sample_id}.no-gapfill.log"

    """
    set -o pipefail

    carve '${faa}' \
        --output '${model}' \
        --solver '${params.solver}' \
        --verbose \
        2>&1 | tee '${log}'
    """

    stub:
    """
    touch '${sample_id}.no-gapfill.sbml'
    echo 'Stub non-gapfill run for ${sample_id}' \
        > '${sample_id}.no-gapfill.log'
    """
}
