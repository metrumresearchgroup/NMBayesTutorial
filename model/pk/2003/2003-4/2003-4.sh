#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2003/2003-4

/opt/NONMEM/nm75/run/nmfe75 2003-4.ctl  2003-4.lst  -parafile=2003-4.pnm -maxlim=2
