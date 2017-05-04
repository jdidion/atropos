
process A {
    input:
    val animal from { [ 'horse', 'cow' ] }

    output:
    set val(animal), file "${animal}.txt" into A_results

    """
    echo -n "${animal}!" > ${animal}.txt
    """
}

process B {
    input:
    val flower from { [ 'weed', 'clover' ] }

    output:
    set val(flower), file "${flower}.txt" into B_results

    """
    echo -n "${flower}!" > ${flower}.txt
    """
}

Channel.empty().
    concat(A_results, B_results).
    set { C_results }

process Combine {
    input:
    set val(item), file(combined) from C_results
    
    output:
    stdout merged
    
    script:
    template "parse_gtime.py"
}

process echoMerged {
    echo true
    
    input:
    val merged
    
    script:
    """
    echo $merged
    """
}