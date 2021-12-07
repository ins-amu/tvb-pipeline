#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: recon-all tool from FreeSurfer
doc: This tool invokes recon-all from FreeSurfer

baseCommand: mri_convert

inputs:
    subject:
        type: string
        inputBinding:
            position: 1
            prefix: -s
    input:
        type: File?
        inputBinding:
            position: 2
            prefix: -i
    all:
        type: boolean
        inputBinding:
            position: 3
            prefix: -all
    parallel:
        type: boolean
        default: true
        inputBinding:
            position: 4
            prefix: -parallel
    openmp:
        type: int
        default: 4
        inputBinding:
            position: 5
            prefix: -openmp

# how to handle directory of outputs..?
# and directory of inputs?
outputs: []
