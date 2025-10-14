#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Enrich a stoichiometric SBML model using BioPathOpt.
 * Container: melclic/biopathopt:latest
 *
 * Example:
 *   nextflow run enrich_model.nf \
 *     --model iML1515.xml \
 *     --taxonomy_id 83333 \
 */


/* ----------------------------
 * Help message
 * ---------------------------- */
def helpMessageEnrich() {
    log.info """
╭────────────────────────────────────────────────────────────────────────╮
│                         BioPathOpt Enrich Model                        │
╰────────────────────────────────────────────────────────────────────────╯

    Usage:
      nextflow run enrich_model.nf --model <SBML_FILE> [options]

    Options:
      --input_model <path>           Path to SBML model (required)
      --taxonomy_id <int>      NCBI Taxonomy ID (default: None)
      --use_progressbar <bool> Show progress bar (default: false)
      --low_memory_mode <bool> Use low memory mode (default: true)
      --output_folder <path>    Path to the folder output (default: biopathopt)
      --help                   Show this help message

    Examples:
      nextflow run enrich_model.nf --model iML1515.xml

    ==========================================
""".stripIndent()
}

/* ----------------------------
 * Parameters
 * ---------------------------- */
params.input_model = null
params.taxonomy_id = 'None'
params.help = false
params.output_folder = "results"

/* ----------------------------
 * Process
 * ---------------------------- */
process enrich_model {
    errorStrategy 'terminate'

    container "melclic/biopathopt:latest"
    //container "biopathopt:cache"

    publishDir (
        path: { "${params.output_folder}/" },
        mode: "copy", 
        pattern: "enriched_model.xml"
    )

    input:
        path model_file
        val tax_id

    output:
        path "enriched_model.xml", emit: enriched_model

    script:
    """
#!/usr/bin/env python3

import os
import pathlib
import biopathopt

# Convert taxonomy_id safely
tax_id_val = None if str(${tax_id}) in ("None", "", "null", "NaN", "nan") else int(${tax_id})

bio_model = biopathopt.ModelBuilder(
    path_to_model="${model_file}",
    taxonomy_id=tax_id_val,
)
bio_model.save_model(file_path="enriched_model.xml")
    """
}

process extract_sink {
    errorStrategy 'terminate'

    container "melclic/biopathopt:latest"
    //container "biopathopt:cache"

    publishDir (
        path: { "${params.output_folder}/" },
        mode: "copy", 
        pattern: "sink.csv"
    )

    input:
        path model_file
    output:
        path "sink.csv", emit: sink

    script:
    """
#!/usr/bin/env python3

import os
import pathlib
import biopathopt
import pandas as pd
import csv

bio_model = biopathopt.ModelBuilder(
    path_to_model="${model_file}",
)
res = []
for m in bio_model.model.metabolites:
    #try to find the mnx
    mnx = m.annotation.get('metanetx.chemical')
    name = None
    if mnx:
        if isinstance(mnx, list):
            if len(mnx)==1:
                name = mnx[0]
            elif len(mnx)>1:
                name = mnx[0]
        elif isinstance(mnx, str):
            name = mnx
    #if not use the name
    if not name:
        name = m.name
    #try to find the inchi
    res.append({'Name': name, 'InChI': m.annotation.get('inchi')})
res = pd.DataFrame.from_records(res)
res = res.drop_duplicates()
res = res.set_index('Name')
res = res.fillna('None')
res.to_csv('sink.csv', quotechar='"', quoting=csv.QUOTE_ALL)
    """
}


/* ----------------------------
 * Workflow
 * ---------------------------- */
workflow enrich {
    if (params.help || !params.input_model) {
        helpMessageEnrich()
        exit 0
    }

    Channel
      .fromPath(params.input_model)
      .ifEmpty { exit 1, "ERROR: Model file not found: ${params.input_model}" }
      .set { ch_model }

    ch_tax = Channel.value(params.taxonomy_id)

    enrich_model(
        ch_model,
        ch_tax,
    )
}

workflow generate_sink {
    if (params.help || !params.input_model) {
        helpMessageEnrich()
        exit 0
    }

    Channel
      .fromPath(params.input_model)
      .ifEmpty { exit 1, "ERROR: Model file not found: ${params.input_model}" }
      .set { ch_model }

    extract_sink(
        ch_model,
    )
}
