#!/bin/bash


if test "$1" = ""; then
    plotDir="output"
else
    if test -d "$1"; then
        plotDir="$1"
    else 
        plotDir=""
        if test -f "$1"; then
            okular "$1" &
        else
            echo "something strange happened"
        fi
    fi
fi

shopt -s nullglob;
if test "$plotDir" != ""; then
    cd "$plotDir"
    for i in *.pdf; do
        okular "$i" &
    done
fi

