#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2004/2004_2

/opt/NONMEM/nm75/run/nmfe75 2004_2.ctl  2004_2.lst  -parafile=2004_2.pnm -maxlim=2
