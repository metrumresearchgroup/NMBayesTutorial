#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2004/2004-1

/opt/NONMEM/nm75/run/nmfe75 2004-1.ctl  2004-1.lst  -parafile=2004-1.pnm -maxlim=2
