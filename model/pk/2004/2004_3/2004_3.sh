#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2004/2004_3

/opt/NONMEM/nm75/run/nmfe75 2004_3.ctl  2004_3.lst  -parafile=2004_3.pnm -maxlim=2
