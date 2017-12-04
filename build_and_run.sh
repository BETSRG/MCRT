#! /bin/bash
export PATH="`pwd`/scripts:$PATH"
makefile.sh -c gfortran -f release -n MCRT
