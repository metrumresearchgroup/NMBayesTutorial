#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2001/2001-4

/opt/NONMEM/nm75/run/nmfe75 2001-4.ctl  2001-4.lst  -parafile=2001-4.pnm -maxlim=2
