#!/bin/bash
# <execname> <fNameSec> <nx> <ny> <nz> <absorb> <dx> <dy> <dz> <dt> <tmax>

if [ $# -lt 1 ]
then
	echo "Usage : $0 <size>"
	exit 1
fi
export ACC_DEVICE_NUM=1;
rsync -a -v ../src/ModelagemFletcher.no.io.exe .
time ./ModelagemFletcher.no.io.exe $1 $1 $1
