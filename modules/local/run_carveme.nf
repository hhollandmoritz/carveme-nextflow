process RUN_CARVEME_DEFAULT {

    tag "${sample_id}:${medium}"

    input:
    tuple val(sample_id),
          path(faa),
          val(medium)

    output:
    tuple val(sample_id),
          val(medium),
          path("${sample_id}.${medium}.sbml", glob: false),
          emit: models

    tuple val(sample_id),
          val(medium),
          path("${sample_id}.${medium}.log", glob: false),
          emit: logs

    script:
    def model = "${sample_id}.${medium}.sbml"
    def carveme_args = params.carveme_args?.toString() ?: ""
    def log   = "${sample_id}.${medium}.log"


    def universe_arg = params.universe_file
        ? "--universe-file '${file(params.universe_file, checkIfExists: true)}'"
        : params.universe
            ? "--universe '${params.universe}'"
            : ""

    """
    set -o pipefail

    carve '${faa}' \
        --output '${model}' \
        --gapfill '${medium}' \
        --init '${medium}' \
        ${universe_arg} \
        --solver '${params.solver}' \
        ${carveme_args} \
        --verbose \
        2>&1 | tee '${log}'
    """
    stub:
    """
    touch '${sample_id}.${medium}.sbml'
    echo 'Stub built-in medium run: ${sample_id} ${medium}' \
        > '${sample_id}.${medium}.log'
    echo 'Universe: ${params.universe ?: ''}' \
        >> '${sample_id}.${medium}.log'
    echo 'Universe file: ${params.universe_file ?: ''}' \
        >> '${sample_id}.${medium}.log'
    echo 'Extra args: ${params.carveme_args ?: ''}' \
        >> '${sample_id}.${medium}.log'
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
          path("${sample_id}.${medium}.sbml", glob: false),
          emit: models

    tuple val(sample_id),
          val(medium),
          path("${sample_id}.${medium}.log", glob: false),
          emit: logs

    script:
    def model = "${sample_id}.${medium}.sbml"
    def carveme_args = params.carveme_args?.toString() ?: ""
    def log   = "${sample_id}.${medium}.log"

    def universe_arg = params.universe_file
    ? "--universe-file '${file(params.universe_file, checkIfExists: true)}'"
    : params.universe
        ? "--universe '${params.universe}'"
        : ""

    """
    set -o pipefail

    carve '${faa}' \
        --output '${model}' \
        --gapfill '${medium}' \
        --init '${medium}' \
        --mediadb '${mediadb}' \
        ${universe_arg} \
        --solver '${params.solver}' \
        ${carveme_args} \
        --verbose \
        2>&1 | tee '${log}'
    """

    stub:
    """
    touch '${sample_id}.${medium}.sbml'
    echo 'Stub custom medium run: ${sample_id} ${medium} ${mediadb}' \
        > '${sample_id}.${medium}.log'
    echo 'Universe: ${params.universe ?: ''}' \
        >> '${sample_id}.${medium}.log'
    echo 'Universe file: ${params.universe_file ?: ''}' \
        >> '${sample_id}.${medium}.log'
    echo 'Extra args: ${params.carveme_args ?: ''}' \
        >> '${sample_id}.${medium}.log'
    """
}

process RUN_CARVEME_NO_GAPFILL {

    tag "${sample_id}:no-gapfill"

    input:
        tuple val(sample_id), 
        path(faa)

    output:
    tuple val(sample_id),
          val('NO_GAPFILL'),
          path("${sample_id}.no-gapfill.sbml", glob: false),
          emit: models

    tuple val(sample_id),
          val('NO_GAPFILL'),
          path("${sample_id}.no-gapfill.log", glob: false),
          emit: logs

    script:
    def model = "${sample_id}.no-gapfill.sbml"
    def carveme_args = params.carveme_args?.toString() ?: ""
    def log   = "${sample_id}.no-gapfill.log"

    def universe_arg = params.universe_file
        ? "--universe-file '${file(params.universe_file, checkIfExists: true)}'"
        : params.universe
            ? "--universe '${params.universe}'"
            : ""

    """
    set -o pipefail

    carve '${faa}' \
        --output '${model}' \
        ${universe_arg} \
        --solver '${params.solver}' \
        ${carveme_args} \
        --verbose \
        2>&1 | tee '${log}'
    """

    stub:
    """
    touch '${sample_id}.no-gapfill.sbml'
    echo 'Stub non-gapfill run for ${sample_id}' \
        > '${sample_id}.no-gapfill.log'
    echo 'Universe: ${params.universe ?: ''}' \
        >> '${sample_id}.no-gapfill.log'
    echo 'Universe file: ${params.universe_file ?: ''}' \
        >> '${sample_id}.no-gapfill.log'
    echo 'Extra args: ${params.carveme_args ?: ''}' \
        >> '${sample_id}.no-gapfill.log'
    """
}
