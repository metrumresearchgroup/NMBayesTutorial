#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2005/2005_3

/opt/NONMEM/nm75/run/nmfe75 2005_3.ctl  2005_3.lst  -parafile=2005_3.pnm -maxlim=2
