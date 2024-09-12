#!/usr/bin/env bash
sphinx-apidoc -o api/ -f ../src/turbomoleio
make html
