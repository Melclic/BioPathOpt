nextflow.enable.dsl=2

/*
 * Help message
 */
def helpMessageSink() {
    log.info """
╭──────────────────────────────────────-------─────────────────────╮
│                        RP2 Sink Extractor                        │
╰─────────────────────────────────────---------------──────────────╯

Options
    --model <path>      Path to SBML model xml or json (required)
    --out   <path>      Path to output sink csv        (required)
    --help              Show this help message

Example
    nextflow run main.nf \\
        --model models/iML1515.xml \\
        --out   results/sink.csv
""".stripIndent()
}

params.model = null
params.out   = null

process EXTRACT_SINK {

    tag { model_file.baseName }

    container "melclic/biopathopt:latest"

    publishDir {
        def outFile   = file(params.out)
        def outFolder = outFile.parent ?: file('.')
        // Use the directory part of out_tar for publishing
        outFolder ? outFolder.toString() : "."
    }, mode: 'copy', overwrite: true

    input:
        path model_file

    output:
        path 'sinkfile.csv', emit: sink_csv

    script:
        """
        python /home/rp2/extract_sink.py \\
            --model ${model_file} \\
            --out   sinkfile.csv
        """
}

workflow {

    if( params.help ) {
        helpMessageSink()
        System.exit(0)
    }

    if( !params.model || !params.out ) {
        log.error "Missing required arguments model and out"
        helpMessageSink()
        System.exit(1)
    }

    model_ch = channel.fromPath(params.model, checkIfExists: true)

    EXTRACT_SINK(model_ch)
}
