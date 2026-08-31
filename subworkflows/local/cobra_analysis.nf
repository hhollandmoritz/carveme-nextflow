include {
    RUN_COBRA_MODEL
    SUMMARIZE_COBRA_MODELS
} from '../../modules/local/cobra_analysis'

/*
 * Input `models` must contain:
 *   tuple(sample_id, medium, model)
 *
 * `cobra_medium` and `bigg_metabolites` are value channels containing either
 * one staged file or an empty list (`[]`).
 */
workflow COBRA_ANALYSIS {
    take:
    models
    cobra_medium
    bigg_metabolites
    run_options
    summary_options

    main:
    RUN_COBRA_MODEL(
        models,
        cobra_medium,
        run_options
    )

    model_tables = RUN_COBRA_MODEL.out.model_table
        .map { sample_id, medium, table -> table }
        .collect()

    summary_tables = RUN_COBRA_MODEL.out.summary_table
        .map { sample_id, medium, table -> table }
        .collect()

    flux_tables = RUN_COBRA_MODEL.out.flux_table
        .map { sample_id, medium, table -> table }
        .collect()

    SUMMARIZE_COBRA_MODELS(
        model_tables,
        summary_tables,
        flux_tables,
        bigg_metabolites,
        summary_options
    )

    emit:
    model_tables = RUN_COBRA_MODEL.out.model_table
    summary_tables = RUN_COBRA_MODEL.out.summary_table
    flux_tables = RUN_COBRA_MODEL.out.flux_table
    comparison = SUMMARIZE_COBRA_MODELS.out.comparison
}
