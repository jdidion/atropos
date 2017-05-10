
process A {
    container "jdidion/atropos_paper_analysis"
    
    input:
    val animal from { [ 'horse_4_rna_20', 'cow_4_rna_20' ] }

    output:
    set val(animal), file("${animal}.txt") into A_results

    """
    /usr/bin/time -v -o ${animal}.txt echo -n "${animal}!"
    """
}

process B {
    container "jdidion/atropos_paper_analysis"
    
    input:
    val flower from { [ 'weed_4_rna_20', 'clover_4_rna_20' ] }

    output:
    set val(flower), file("${flower}.txt") into B_results

    """
    /usr/bin/time -v -o ${flower}.txt echo -n "${flower}!"
    """
}

Channel.empty().
    concat(A_results, B_results).
    set { C_results }

process ParseCombined {
    container "jdidion/python_bash"
    
    input:
    set val(item), file(timing) from C_results
    
    output:
    stdout merged
    
    script:
    template "parse_gtime.py"
}

process echoMerged {
    echo true
    
    input:
    val mergedRows from merged.toList()
    
    script:
    mergedData = mergedRows.join("")
    
    """
    echo '$mergedData'
    """
}