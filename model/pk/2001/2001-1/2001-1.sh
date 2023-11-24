#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2001/2001-1

/opt/NONMEM/nm75/run/nmfe75 2001-1.ctl  2001-1.lst  -parafile=2001-1.pnm -maxlim=2
