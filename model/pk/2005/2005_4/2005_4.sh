#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2005/2005_4

/opt/NONMEM/nm75/run/nmfe75 2005_4.ctl  2005_4.lst  -parafile=2005_4.pnm -maxlim=2
