#!/bin/bash

#$ -wd /data/NMBayesTutorial/model/pk/1000-FOCE

/opt/NONMEM/nm75/run/nmfe75 1000-FOCE.ctl  1000-FOCE.lst  -maxlim=2
