#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/2005/2005-3

/opt/NONMEM/nm75/run/nmfe75 2005-3.ctl  2005-3.lst  -parafile=2005-3.pnm -maxlim=2
