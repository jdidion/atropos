
process A {
    input:
    val animal from { [ 'horse', 'cow' ] }

    output:
    file "${animal}.txt" into A_results

    """
    echo -n "${animal}!" > ${animal}.txt
    """
}

process B {
    input:
    val flower from { [ 'weed', 'clover' ] }

    output:
    file "${flower}.txt" into B_results

    """
    echo -n "${flower}!" > ${flower}.txt
    """
}

Channel.empty().
    concat(A_results, B_results).
    set { C_results }

process Combine {
    echo true

    input:
    file c from C_results.toList()

    script:
    
    """
    cat $c
    """
}