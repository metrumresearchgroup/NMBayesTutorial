#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2001/2001_3

/opt/NONMEM/nm75/run/nmfe75 2001_3.ctl  2001_3.lst  -parafile=2001_3.pnm -maxlim=2
