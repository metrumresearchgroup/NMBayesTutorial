#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2002/2002-4

/opt/NONMEM/nm75/run/nmfe75 2002-4.ctl  2002-4.lst  -parafile=2002-4.pnm -maxlim=2
