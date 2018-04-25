#!/bin/bash
# <execname> <fNameSec> <nx> <ny> <nz> <absorb> <dx> <dy> <dz> <dt> <tmax>

if [ $# -lt 1 ]
then
	echo "Usage : $0 <size>"
	exit 1
fi

rsync -a -v ../src/ModelagemFletcher.exe .
time ./ModelagemFletcher.exe "ISO" $1 $1 $1 32 12.5 12.5 12.5 0.0010 0.2

