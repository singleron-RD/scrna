#!/usr/bin/bash

SRC_PATTERN="nextflow_schema.json"
if git diff --cached --name-only | grep --quiet "$SRC_PATTERN"
then
    nf-core schema docs -f -o docs/parameters.md
fi