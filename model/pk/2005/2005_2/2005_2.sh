#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2005/2005_2

/opt/NONMEM/nm75/run/nmfe75 2005_2.ctl  2005_2.lst  -parafile=2005_2.pnm -maxlim=2
