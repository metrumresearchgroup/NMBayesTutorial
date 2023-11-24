#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2004/2004-4

/opt/NONMEM/nm75/run/nmfe75 2004-4.ctl  2004-4.lst  -parafile=2004-4.pnm -maxlim=2
