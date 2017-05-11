
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
    stdout timingParsed
    
    script:
    "parse_gtime.py -i $timing -p $item"
}

process ShowPerformance {
    container "jdidion/python_bash"
    
    input:
    val parsedRows from timingParsed.toList()
    
    output:
    file "timing.txt"
    file "timing.tex"
    file "timing.svg"
    file "timing.pickle"
    
    script:
    data = parsedRows.join("")
    """
    echo '$data' | show_performance.py -n foo -c bar -o timing -f txt tex svg pickle
    """
}
