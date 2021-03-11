#!/usr/bin/env bash
sphinx-apidoc -o api/ -f ../turbomoleio
make html
