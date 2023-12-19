#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2003/2003-2

/opt/NONMEM/nm75/run/nmfe75 2003-2.ctl  2003-2.lst  -parafile=2003-2.pnm -maxlim=2
