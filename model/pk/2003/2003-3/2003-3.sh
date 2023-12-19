#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2003/2003-3

/opt/NONMEM/nm75/run/nmfe75 2003-3.ctl  2003-3.lst  -parafile=2003-3.pnm -maxlim=2
