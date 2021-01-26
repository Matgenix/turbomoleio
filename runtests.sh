#!/bin/bash

POSITIONAL=()
#while [[ $# -gt 0 ]]
#do
key="$1"

case $key in
    -u|--unit)
    pytest -m 'not integration' --ignore=BASFflows -v --cov-report term-missing:skip-covered --cov-config .coveragerc_unit --cov=turbomoleio
    ;;
    -i|--integration)
    pytest -m 'integration' --ignore=BASFflows -v --cov-report term-missing:skip-covered --cov-config .coveragerc_integration --cov=turbomoleio
    ;;
    *)
    pytest --ignore=BASFflows -v --cov-report term-missing:skip-covered --cov-config .coveragerc --cov=turbomoleio
    ;;
esac
#done
