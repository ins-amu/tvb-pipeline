#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: Mrconvert tool from MRtrix3
doc: This tool invokes mrconvert from MRtrix3

baseCommand: mrconvert

inputs:
    input:
        type: File
        inputBinding:
            position: 1
    output_fname:
        type: string
        default: image.mif
        inputBinding:
            position: 2

outputs:
    output:
        type: File
        outputBinding:
            glob: $(inputs.output_fname)
