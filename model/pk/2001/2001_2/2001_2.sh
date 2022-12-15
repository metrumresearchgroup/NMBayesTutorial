#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2001/2001_2

/opt/NONMEM/nm75/run/nmfe75 2001_2.ctl  2001_2.lst  -parafile=2001_2.pnm -maxlim=2
