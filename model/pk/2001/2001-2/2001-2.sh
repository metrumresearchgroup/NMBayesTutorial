#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2001/2001-2

/opt/NONMEM/nm75/run/nmfe75 2001-2.ctl  2001-2.lst  -parafile=2001-2.pnm -maxlim=2
