#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2001/2001-3

/opt/NONMEM/nm75/run/nmfe75 2001-3.ctl  2001-3.lst  -parafile=2001-3.pnm -maxlim=2
