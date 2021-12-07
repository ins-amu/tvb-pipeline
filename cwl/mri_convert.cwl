#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: mri_convert tool from FreeSurfer
doc: This tool invokes mri_convert from MRtrix3

baseCommand: mri_convert

inputs:
    resample_type:
        type: string
        default: nearest
        inputBinding:
            position: 1
            prefix: -rt
    out_orientation:
        type: string
        default: RAS
        inputBinding:
            position: 2
            prefix: --out_orientation
    input:
        type: File
        inputBinding:
            position: 3
    output_fname:
        type: string
        default: image.nii.gz
        inputBinding:
            position: 4

outputs:
    output:
        type: File
        outputBinding:
            glob: $(inputs.output_fname)
