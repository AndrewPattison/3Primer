#!/usr/bin/env bash

makeblastdb -in reference.fa -parse_seqids -dbtype nucl

#samtools faidx reference.fa 