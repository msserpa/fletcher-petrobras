#!/bin/bash
# <execname> <fNameSec> <nx> <ny> <nz> <absorb> <dx> <dy> <dz> <dt> <tmax>

if [ $# -lt 2 ]
then
	echo "Usage : $0 <size> <I/O step>"
	exit 1
fi
export ACC_DEVICE_NUM=1;
rsync -a -v ../src/ModelagemFletcher.io.exe .
time ./ModelagemFletcher.io.exe $1 $1 $1 $2

