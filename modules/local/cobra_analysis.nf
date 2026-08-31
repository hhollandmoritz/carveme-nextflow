/*
 * Run COBRApy optimization on one CarveMe model per task, then aggregate
 * the standardized TSV outputs in a separate task.
 *
 * Expected model input tuple:
 *   tuple(sample_id, medium, model)
 */

process RUN_COBRA_MODEL {
    tag "${sample_id}:${medium}"
    label 'process_cobra_model'

    input:
    tuple val(sample_id), val(medium), path(model)
    path cobra_medium
    val options

    output:
    tuple val(sample_id), val(medium), path('*.model.tsv'),   emit: model_table
    tuple val(sample_id), val(medium), path('*.summary.tsv'), emit: summary_table
    tuple val(sample_id), val(medium), path('*.fluxes.tsv'),  emit: flux_table

    script:
    def prefix = "${sample_id}.${medium}"
        .replaceAll(/[^A-Za-z0-9._-]+/, '_')
        .replaceAll(/^[._-]+|[._-]+$/, '')

    if (!prefix) {
        error "COBRA output prefix is empty for sample=${sample_id}, medium=${medium}"
    }

    def args = []
    args << "--method '${options.method ?: 'pfba'}'"
    args << "--fraction-of-optimum '${options.fraction_of_optimum == null ? 1.0 : options.fraction_of_optimum}'"
    args << "--flux-threshold '${options.flux_threshold == null ? 1e-9 : options.flux_threshold}'"

    if (options.objective) {
        args << "--objective '${options.objective}'"
    }
    if (options.solver) {
        args << "--solver '${options.solver}'"
    }
    if (cobra_medium) {
        args << "--medium '${cobra_medium}'"
        args << "--medium-mode '${options.medium_mode ?: 'replace'}'"
    }
    if (options.ignore_missing_medium_reactions) {
        args << '--ignore-missing-medium-reactions'
    }
    if (options.fail_on_nonoptimal) {
        args << '--fail-on-nonoptimal'
    }
    if (options.extra_args) {
        args << options.extra_args.toString()
    }

    def argument_string = args.join(' \\\n        ')

    """
    run_cobra_model.py \\
        '${model}' \\
        --output-dir . \\
        --name '${sample_id}.${medium}' \\
        --prefix '${prefix}' \\
        --label 'sample_id=${sample_id}' \\
        --label 'medium=${medium}' \\
        ${argument_string}
    """

    stub:
    def stub_prefix = "${sample_id}.${medium}"
        .replaceAll(/[^A-Za-z0-9._-]+/, '_')
        .replaceAll(/^[._-]+|[._-]+$/, '')

    """
    printf 'model_name\tresult_id\tsample_id\tmedium\tstatus\tmaximum_biomass\n%s\t%s\t%s\t%s\toptimal\t1.0\n' \\
        '${sample_id}.${medium}' '${stub_prefix}' '${sample_id}' '${medium}' \\
        > '${stub_prefix}.model.tsv'

    printf 'model_name\tresult_id\tsample_id\tmedium\treaction\tmetabolite\tboundary_type\tfactor\tflux\tabsolute_flux\tdirection\n%s\t%s\t%s\t%s\tEX_stub_e\tstub_e\texchange\t-1\t1.0\t1.0\tuptake\n' \\
        '${sample_id}.${medium}' '${stub_prefix}' '${sample_id}' '${medium}' \\
        > '${stub_prefix}.summary.tsv'

    printf 'model_name\tresult_id\tsample_id\tmedium\treaction\treaction_name\tflux\tabsolute_flux\tlower_bound\tupper_bound\tequation\tgene_reaction_rule\n%s\t%s\t%s\t%s\tBIOMASS\tBiomass\t1.0\t1.0\t0\t1000\t\t\n' \\
        '${sample_id}.${medium}' '${stub_prefix}' '${sample_id}' '${medium}' \\
        > '${stub_prefix}.fluxes.tsv'
    """
}


process SUMMARIZE_COBRA_MODELS {
    tag 'all_models'
    label 'process_cobra_summary'

    input:
    path model_tables
    path summary_tables
    path flux_tables
    path bigg_metabolites
    val options

    output:
    path 'comparison', emit: comparison

    script:
    def requested_outputs = options.outputs ?: [
        'growth-table',
        'summary-long',
        'summary-matrix',
        'summary-heatmap'
    ]
    def row_columns = options.row_columns ?: ['sample_id', 'medium']
    def boundary_types = options.summary_boundary_types ?: ['exchange']

    def args = []
    args << '--outputs ' + requested_outputs.collect { "'${it}'" }.join(' ')
    args << '--row-columns ' + row_columns.collect { "'${it}'" }.join(' ')
    args << "--row-separator '${options.row_separator ?: ' | '}'"
    args << "--summary-column '${options.summary_column ?: 'metabolite'}'"
    args << '--summary-boundary-types ' + boundary_types.collect { "'${it}'" }.join(' ')
    args << "--min-abs-flux '${options.min_abs_flux == null ? 1e-9 : options.min_abs_flux}'"
    args << "--top-n '${options.top_n == null ? 50 : options.top_n}'"
    args << "--heatmap-transform '${options.heatmap_transform ?: 'signed-log1p'}'"
    args << "--figure-format '${options.figure_format ?: 'png'}'"
    args << "--dpi '${options.dpi ?: 300}'"

    if (options.include_nonoptimal) {
        args << '--include-nonoptimal'
    }

    /*
     * This option is added only when a BiGG file is supplied. It therefore
     * remains compatible with the non-BiGG summarizer when the parameter is null.
     */
    if (bigg_metabolites) {
        args << "--bigg-metabolites '${bigg_metabolites}'"
    }
    if (options.extra_args) {
        args << options.extra_args.toString()
    }

    def argument_string = args.join(' \\\n        ')

    """
    mkdir -p comparison

    summarize_cobra_models.py \\
        --results . \\
        --output-dir comparison \\
        ${argument_string}
    """

    stub:
    """
    mkdir -p comparison

    printf 'model_name\tresult_id\tstatus\tmaximum_biomass\nstub\tstub\toptimal\t1.0\n' \\
        > comparison/growth_table.tsv
    printf 'model_name\tresult_id\tmetabolite\tflux\n' \\
        > comparison/summary_long.tsv
    printf 'model_name\tstub_e\nstub\t1.0\n' \\
        > comparison/summary_flux_matrix.tsv
    printf 'model_name\tBIOMASS\nstub\t1.0\n' \\
        > comparison/reaction_flux_matrix.tsv
    : > comparison/summary_flux_heatmap.png
    : > comparison/reaction_flux_heatmap.png
    """
}
