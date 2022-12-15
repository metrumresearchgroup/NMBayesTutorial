#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2004/2004_4

/opt/NONMEM/nm75/run/nmfe75 2004_4.ctl  2004_4.lst  -parafile=2004_4.pnm -maxlim=2
