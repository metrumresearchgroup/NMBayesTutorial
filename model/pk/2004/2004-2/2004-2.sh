#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2004/2004-2

/opt/NONMEM/nm75/run/nmfe75 2004-2.ctl  2004-2.lst  -parafile=2004-2.pnm -maxlim=2
