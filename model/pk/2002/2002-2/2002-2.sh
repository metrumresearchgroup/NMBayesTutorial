#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2002/2002-2

/opt/NONMEM/nm75/run/nmfe75 2002-2.ctl  2002-2.lst  -parafile=2002-2.pnm -maxlim=2
