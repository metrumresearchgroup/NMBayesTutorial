#!/bin/bash

#$ -wd /data/bbr-nonmem-poppk-bayes/model/pk/2001/2001_4

/opt/NONMEM/nm75/run/nmfe75 2001_4.ctl  2001_4.lst  -parafile=2001_4.pnm -maxlim=2
