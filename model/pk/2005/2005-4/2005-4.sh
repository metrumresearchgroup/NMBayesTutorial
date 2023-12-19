#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2005/2005-4

/opt/NONMEM/nm75/run/nmfe75 2005-4.ctl  2005-4.lst  -parafile=2005-4.pnm -maxlim=2
