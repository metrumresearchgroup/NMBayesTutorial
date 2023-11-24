#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2002/2002-1

/opt/NONMEM/nm75/run/nmfe75 2002-1.ctl  2002-1.lst  -parafile=2002-1.pnm -maxlim=2
