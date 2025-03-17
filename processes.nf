process downloadSCOPe {
    output:
        file("globulars") into globulars_ch
        file("moltens")   into moltens_ch
    script:
        """
        ./download_scop.sh
        """
}

process foldseekClustering {
    input:
        tuple val(label), file(folder)
    output:
        tuple val(label), file("${label}.tsv")
    script:
        """
        ./foldseek_cluster.sh ${folder} AFDB > ${label}.tsv
        """
}

// Define additional processes (parseTSV, extractSequences, clusterSequences) here...
