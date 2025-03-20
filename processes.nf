#!/usr/bin/env nextflow

params.scop_db = "/store/EQUIPES/BIM/MEMBERS/simon.herman/MicroData/SCOPe_database"

// Process to generate SCOPe data
process generateSCOPeData {

    // Outputs: two FASTA files and two tagged pdb directories
    output:
    path "final_output/abcde_full_sequences.fasta", emit: abcde_fasta
    path "final_output/g_full_sequences.fasta", emit: g_fasta
    tuple val("abcde"), path("final_output/pdb_files/abcde"), emit: abcde_pdb
    tuple val("g"), path("final_output/pdb_files/g"), emit: g_pdb

    script:
    """
    # Run your SCOPe Python script.
    # It must create the outputs in the final_output folder as expected.
    python SCOPe.py
    """
}

process foldseekClustering {

    clusterOptions '-q bim -l ncpus=60 -l mem=250gb -l walltime=48:00:00'

    input:
    tuple val(label), path(pdb_dir)

    output:
    tuple val(label), path("${label}_vs_afdb_micro.tsv")

    script:
    """
    #!/bin/bash

    foldseek="/store/EQUIPES/BIM/MEMBERS/ambre.baumann/Outils/foldseek/bin/foldseek"

    mkdir -p globular moltens database database/tmp

    rm -rf globular/* moltens/*

    if [ "\$label" = "abcde" ]; then
        cp ${pdb_dir}/*.ent globular/
    elif [ "\$label" = "g" ]; then
        cp ${pdb_dir}/*.ent moltens/
    else
        echo "Unknown label: \$label"
        exit 1
    fi

    \$foldseek createdb pdbs_microproteins/ database/afdb_20_to_100_fs_db
    \$foldseek createindex database/afdb_20_to_100_fs_db database/tmp

    if [ "\$label" = "abcde" ]; then
        \$foldseek createdb globular/ database/globular_db
        \$foldseek createindex database/globular_db database/tmp
        \$foldseek easy-search database/globular_db database/afdb_20_to_100_fs_db \${label}_vs_afdb_micro.tsv database/tmp -e 0.1 --num-iterations 2 --exhaustive-search 1 --format-output query,target,alntmscore,qtmscore,ttmscore,evalue,pident,qlen,qstart,qend,tlen,tstart,tend,qcov,tcov,prob,lddt,lddtfull --format-mode 4 --threads 60
    elif [ "\$label" = "g" ]; then
        \$foldseek createdb moltens/ database/moltens_db
        \$foldseek createindex database/moltens_db database/tmp
        \$foldseek easy-search database/moltens_db database/afdb_20_to_100_fs_db \${label}_vs_afdb_micro.tsv database/tmp -e 0.1 --num-iterations 2 --exhaustive-search 1 --format-output query,target,alntmscore,qtmscore,ttmscore,evalue,pident,qlen,qstart,qend,tlen,tstart,tend,qcov,tcov,prob,lddt,lddtfull --format-mode 4 --threads 60
    fi
    """
}


workflow {

    result = generateSCOPeData()

    pdb_dirs = Channel.merge(result.abcde_pdb, result.g_pdb)

    foldseek_results = pdb_dirs | foldseekClustering

    foldseek_results.view { "Foldseek result for label: \${it[0]} -> file: \${it[1]}" }
}
