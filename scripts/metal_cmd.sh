#!/bin/bash

# Define the output script file
output_script="scripts/meta_analysis_script.sh"

# Create a new script file
{
    # Write the common header content
    echo "SCHEME STDERR"
    echo "COLUMNCOUNTING LENIENT"
    echo "SEPARATOR COMMA"
    echo "PVALUELABEL bacon.pval"
    echo "EFFECTLABEL bacon.es"
    echo "STDERRLABEL bacon.se"
    echo "MARKER cpgid"
    echo "WEIGHT n"

    # Process each provided input file path as PROCESSFILE commands
    while [ "$#" -gt 1 ]; do
        input_file="$1"
        echo "PROCESSFILE $input_file"
        shift
    done

    # Write the OUTFILE command using the last argument as the output file path
    output_file="$1"
    echo "OUTFILE $output_file  .txt"
    
    # Write the remaining content
    echo "ANALYZE HETEROGENEITY"
    echo "CLEAR"
} > "$output_script"

# Make the output script executable
chmod +x "$output_script"

echo "Script '$output_script' created and made executable."