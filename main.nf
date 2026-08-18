nextflow.enable.dsl = 2

include { 
    RUN_CARVEME_DEFAULT
    RUN_CARVEME_CUSTOM
    RUN_CARVEME_NO_GAPFILL
 } from './modules/local/run_carveme'

workflow {

    log.info "Project directory: ${projectDir}"
    /* 
    * Checks on universe specifications. 
    * 1) only one of universe options specified at a time
    */
    if (params.universe && params.universe_file) {
        error "Specify either --universe or --universe_file, not both."
    }
    /*
     * Check if the universe file exists if specified.
     */ 
    if (params.universe_file) {
        file(params.universe_file, checkIfExists: true)
    }

    /*
     * Emit:
     * tuple(sample_id, faa)
     */
    genomes_ch = channel
        .fromPath(params.input, checkIfExists: true)
        .map { faa ->
            tuple(faa.baseName, faa)
        }

    /*
     * Optional/modifiable arguments 
     * (media, extra carveme arguments)
     */
    media_config = params.media_databases ?: [:]

    /*
     * Run without gap-filling when:
     *
     * 1. No media configuration was supplied, or
     * 2. NO_GAPFILL was explicitly requested.
     */
    run_no_gapfill = media_config.isEmpty() ||
        media_config.containsKey('NO_GAPFILL')

    /*
     * Remove the special non-media condition before processing
     * actual media.
     */
    gapfill_config = media_config.findAll { medium, mediadb ->
        medium.toString() != 'NO_GAPFILL'
    }

    run_gapfill = !gapfill_config.isEmpty()

    /*
     * Non-gap-filled models
     */
    if (run_no_gapfill) {
        RUN_CARVEME_NO_GAPFILL(genomes_ch)
    }

    /*
     * Gap-filled models
     */
    if (run_gapfill) {

        /*
         * Built-in media have no custom database.
         */
        default_media = gapfill_config
            .findAll { medium, mediadb -> mediadb == null }
            .keySet()
            .collect { medium ->
                medium.toString()
            }

        /*
         * Custom media have an associated database file.
         */
        custom_media = gapfill_config
            .findAll { medium, mediadb -> mediadb != null }
            .collect { medium, mediadb ->
                tuple(
                    medium.toString(),
                    file(mediadb, checkIfExists: true)
                )
            }

        /*
         * Emit:
         * tuple(sample_id, faa, medium)
         */
        default_jobs_ch = genomes_ch
            .combine(channel.fromList(default_media))
            .map { sample_id, faa, medium ->
                tuple(sample_id, faa, medium)
            }

        /*
         * Emit:
         * tuple(sample_id, faa, medium, mediadb)
         */
        custom_jobs_ch = genomes_ch
            .combine(channel.fromList(custom_media))
            .map { sample_id, faa, medium, mediadb ->
                tuple(sample_id, faa, medium, mediadb)
            }

        RUN_CARVEME_DEFAULT(default_jobs_ch)
        RUN_CARVEME_CUSTOM(custom_jobs_ch)

        gapfill_models_ch = RUN_CARVEME_DEFAULT.out.models
            .mix(RUN_CARVEME_CUSTOM.out.models)

        gapfill_logs_ch = RUN_CARVEME_DEFAULT.out.logs
            .mix(RUN_CARVEME_CUSTOM.out.logs)
    }

    /*
     * Create unified output channels.
     */
    if (run_no_gapfill && run_gapfill) {

        models_ch = RUN_CARVEME_NO_GAPFILL.out.models
            .mix(gapfill_models_ch)

        logs_ch = RUN_CARVEME_NO_GAPFILL.out.logs
            .mix(gapfill_logs_ch)

    } else if (run_no_gapfill) {

        models_ch = RUN_CARVEME_NO_GAPFILL.out.models
        logs_ch   = RUN_CARVEME_NO_GAPFILL.out.logs

    } else {

        models_ch = gapfill_models_ch
        logs_ch   = gapfill_logs_ch
    }
}
