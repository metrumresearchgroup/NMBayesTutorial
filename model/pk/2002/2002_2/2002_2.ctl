$PROBLEM From bbr: see 2002.yaml for details

$INPUT C ID TIME AMT DV CMT DOSE EVID MDV ADDL II SS FORM WT PER BLQ NUM

$DATA ../../../../data/derived/atorvWrkShop3.csv IGNORE=(C='C', BLQ.EQ.1)

$SUBROUTINES ADVAN4 TRANS4

$PK

NWT = LOG(WT/70) ; normalized to 70 kg adult
MU_1 = THETA(1) + NWT * 0.75 ; CL
MU_2 = THETA(2) + NWT        ; V2
MU_3 = THETA(3) + NWT * 0.75 ; Q
MU_4 = THETA(4) + NWT        ; V3
MU_5 = THETA(5)              ; ka

CL = EXP(MU_1 + ETA(1)) ; atv cl 
V2 = EXP(MU_2 + ETA(2))
Q  = EXP(MU_3 + ETA(3))
V3 = EXP(MU_4 + ETA(4))
KA = EXP(MU_5 + ETA(5))

S2 = V2/1000          ; DOSE IN mmol & CONC IN nM

$THETA
( 5.65)     ;  THETA(1) TOTAL CL ATV
( 5.6 )      ;  THETA(2) V2
( 4.29)     ;  THETA(3) Q
( 7.07)     ;  THETA(4) V3
( -1.33)    ;  THETA(5) KA

$OMEGA BLOCK(2) ;INITIAL values of OMEGAs
0.1  ; cl
0.003 0.2  ; v2

$OMEGA 0 FIX
$OMEGA 0 FIX
$OMEGA 0 FIX
$OMEGA 0 FIX
$OMEGA 0 FIX

$SIGMA  ;Initial value of SIGMA
0.1; prop

$ERROR

IPRED = F;
Y = IPRED*(1+ERR(1)) ; ATV

$PRIOR NWPRI

$THETAP          ; Prior information of THETAS
( 5.65  FIX)     ;  THETA(1) TOTAL CL ATV
( 5.6  FIX)      ;  THETA(2) V2
( 4.29  FIX)     ;  THETA(3) Q
( 7.07  FIX)     ;  THETA(4) V3
( -1.33  FIX)    ;  THETA(5) KA


$THETAPV BLOCK(5)     ;  variances for priors on THETAS (var-cov)
2 FIX ; CL weakly informative
0.00 2  ; V2 weakly informative
0.00 0.00 2  ; Q weakly informative
0.00 0.00 0.00 .0156   ; V3 SE**2 informative
0.00 0.00 0.00 0.00 0.0684    ; KA SE**2 informative


$OMEGAP BLOCK(2)   ; prior information for OMEGA
.1003   FIXED         ; cl
.06034    .8477       ; v2

$OMEGAPD (4.0, FIXED)     ; df for OMEGA prior

$SIGMAP
0.1 FIX           ; 1 prop error

$SIGMAPD (1 FIX)

$EST METHOD=CHAIN FILE=../2002.CHN NSAMPLE=0 ISAMPLE=2 SEED=1 CTYPE=0 IACCEPT=0.25 DF=10 DFS=0

$EST MAXEVAL=9999 METHOD=NUTS INTERACTION AUTO=2 CTYPE=0 OLKJDF=2 OVARF=1 NUTS_DELTA=0.95 NBURN = 1000 NITER = 2000 SEED=2 PRINT=1 NOABORT THIN = 1 RANMETHOD=P PARAFPRINT=10000 BAYES_PHI_STORE=1

$TABLE NUM CL V2 Q KA V3 ETAS(1:LAST) EPRED IPRED NPDE EWRES NOPRINT ONEHEADER FILE=./2002.tab

