$PROB 104; NM 7.4 way of calling for individual posteriors

$INPUT C NUM ID TIME SEQ CMT EVID AMT DV AGE WT HT EGFR ALB BMI SEX AAG
  SCR AST ALT CP TAFD TAD LDOS MDV BLQ PHASE STUDYN

$DATA ../../../../data/derived/analysis3.csv IGNORE=(C='C', BLQ=1)

$ABBR DECLARE INTEGER FIRST_WRITE_IPAR

$SUBROUTINES ADVAN4 TRANS4

$PK
; Request extra information for Bayesian analysis.
; An extra call will then be made for accepted samples 
include '/opt/NONMEM/nm74/util/nonmem_reserved_general'
BAYES_EXTRA_REQUEST=1

;log transformed PK parms
;note, WT is non-time-dependent
 
V2WT  = LOG(WT/70)
CLWT  = LOG(WT/70)*0.75
V3WT  = LOG(WT/70)
QWT   = LOG(WT/70)*0.75

MU_1 = THETA(1) 
MU_2 = THETA(2) + V2WT 
MU_3 = THETA(3) + CLWT
MU_4 = THETA(4) + V3WT
MU_5 = THETA(5) + QWT

" CALL EXTRASEND()

KA   = EXP(MU_1+ETA(1))
V2   = EXP(MU_2+ETA(2))
CL   = EXP(MU_3+ETA(3))
V3   = EXP(MU_4+ETA(4))
Q    = EXP(MU_5+ETA(5)) 

S2 = V2/1000 ; dose in mg, conc in ng/mL

;Table to write out individual level ETA values
IF(BAYES_EXTRA==1 .AND. NDREC==1 .AND. ITER_REPORT>0) THEN 
  IF(FIRST_WRITE_IPAR==0) THEN
    " OPEN(unit=50,FILE='./ipar'//TFI(PNM_NODE_NUMBER)//'.txt') 
    FIRST_WRITE_IPAR=1
  ENDIF
  " WRITE(50,'(I12,1X,F14.0,10(1X,1PG19.10E3))') ITER_REPORT, ID, &
  " ETA(1), ETA(2), ETA(3), ETA(4), ETA(5)
ENDIF

;Table to write out individual level objective values (ppd)
IF(BAYES_EXTRA==1 .AND. NIREC==1 .AND. NDREC==1 .AND. &
   ITER_REPORT>0) THEN 
  IF(FIRST_WRITE_PAR==0) THEN
    " OPEN(unit=52,FILE='./iobj.txt')  
    FIRST_WRITE_PAR=1
  ENDIF
  " WRITE(52,'(I12,1X,500(1X,1PG19.10E3))') ITER_REPORT, &  
  " OBJI
ENDIF

$ERROR
IPRED = F 
Y = IPRED*(1+EPS(1))

$THETA
; log values
  (0.5)         ; 1 KA (1/hr) - 1.5
  (3.5)         ; 2 V2 (L) - 60
  (1)           ; 3 CL (L/hr) - 3.5
  (4)           ; 4 V3 (L) - 70
  (2)           ; 5 Q  (L/hr) - 4

$OMEGA BLOCK(3)
  0.2           ; ETA(KA)
  0.01 0.2      ; ETA(V2)
  0.01 0.01 0.2 ; ETA(CL)

$OMEGA
  0.01 FIX    ; ETA(V3)
  0.01 FIX    ; ETA(Q)

$SIGMA
  0.05          ; 1 pro error
  
$PRIOR NWPRI

$THETAP
  (0.5) FIX        ; 1 KA (1/hr) - 1.5
  (3.5) FIX        ; 2 V2 (L) - 60
  (1)   FIX        ; 3 CL (L/hr) - 3.5
  (4)   FIX        ; 4 V3 (L) - 70
  (2)   FIX        ; 5 Q  (L/hr) - 4
$THETAPV BLOCK(5) VALUES(10, 0) FIX

$OMEGAP BLOCK(3) VALUES(0.2, 0.01) FIX

$OMEGAPD (3 FIX)

$SIGMAP
  0.05 FIX           ; 1 pro error

$SIGMAPD (1 FIX)

; CHAIN options:
;   CTYPE=0: initial estimates for THETA are sampled from a uniform
;     distribution between (1-IACCEPT)*THETA and (1+IACCEPT)*THETA)
;   CTYPE=2: initial estimates for THETA are from a normal distribution with
;     mean from the initial estimate in $THETA and variance from $THETAPV
;   DF=0: initial estimates for OMEGA come from Wishart distribution using
;     values in $OMEGA and degrees of freedom equal to dimensions of OMEGA
;   DFS=0: initial estimates for SIGMA come from Wishart distribution using
;     values in $SIGMA and degrees of freedom equal to dimensions of SIGMA
$EST METHOD=CHAIN FILE=1000.chn NSAMPLE=4 ISAMPLE=0 SEED=1 CTYPE=0 IACCEPT=0.3 DF=10 DFS=0
;$EST METHOD=BAYES SEED=1 NBURN=1000 NITER=10000 PRINT=10 MSFO=./1000.msf RANMETHOD=P PARAFPRINT=10000 BAYES_PHI_STORE=1

;$TABLE NUM CL V2 Q V3 KA ETAS(1:LAST) EPRED IPRED NPDE EWRES NOPRINT ONEHEADER FILE=1000.tab RANMETHOD=P
