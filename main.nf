#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include processes from processes.nf
include { downloadSCOPe; foldseekClustering } from './processes.nf'


workflow {
    
    // Launch the download process
    downloadSCOPe()
    
    // Label the folders for downstream use
    globulars_labeled = globulars_ch.map { folder -> tuple('globulars', folder) }
    moltens_labeled   = moltens_ch.map   { folder -> tuple('moltens', folder) }
    
    // Run Foldseek clustering for both folders
    foldseek_tsv_glob = globulars_labeled  | foldseekClustering
    foldseek_tsv_molt = moltens_labeled    | foldseekClustering

    // Continue with the rest of your workflow (parsing TSV, extracting sequences, clustering, etc.)
    
    // For example:
    // parsed_txt_all = (foldseek_tsv_glob.mix(foldseek_tsv_molt)) | parseTSV
    // ... and so on
}
