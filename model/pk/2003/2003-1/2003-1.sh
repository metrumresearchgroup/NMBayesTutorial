#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2003/2003-1

/opt/NONMEM/nm75/run/nmfe75 2003-1.ctl  2003-1.lst  -parafile=2003-1.pnm -maxlim=2
