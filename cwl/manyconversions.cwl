#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

inputs:
    input_image: File

outputs:
    output_image:
        type: File
        outputSource: convert3/output

steps:
    convert1:
        run: mrconvert.cwl
        in:
            input: input_image
        out: [output]
    convert2:
        run: mrconvert.cwl
        in:
            input: convert1/output
        out: [output]
    convert3:
        run: mrconvert.cwl
        in:
            input: convert2/output
        out: [output]
