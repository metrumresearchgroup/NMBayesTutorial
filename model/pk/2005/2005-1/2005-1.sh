#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2005/2005-1

/opt/NONMEM/nm75/run/nmfe75 2005-1.ctl  2005-1.lst  -parafile=2005-1.pnm -maxlim=2
