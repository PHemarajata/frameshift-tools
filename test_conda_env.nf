#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TestProcess1 {
    conda 'jq python samtools'
    
    script:
    """
    echo "Process 1: Using conda env: \$CONDA_DEFAULT_ENV"
    which jq || echo "jq not found"
    """
}

process TestProcess2 {
    conda 'jq python samtools'
    
    script:
    """
    echo "Process 2: Using conda env: \$CONDA_DEFAULT_ENV"
    which python || echo "python not found"
    """
}

workflow {
    TestProcess1()
    TestProcess2()
}