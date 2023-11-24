#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2004/2004-3

/opt/NONMEM/nm75/run/nmfe75 2004-3.ctl  2004-3.lst  -parafile=2004-3.pnm -maxlim=2
