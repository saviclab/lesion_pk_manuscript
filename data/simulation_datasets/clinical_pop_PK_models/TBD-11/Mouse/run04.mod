;; 1. Based on: 02
;; 2. Description: mCLB073 mouse plasma pk 1-cmt, proportional error
;; x1. Author: user
;; Description: mCLB073 mouse plasma pk 1-cmt, proportional error
;; Author: aebustion
$PROBLEM    mouse plasma PK 1cmt prop error
$INPUT      ID TIME DOSEMGKG DVNGML WTG WTKG DV AMT EVID MDV AMTUG
$DATA      data_03.25.24.csv IGNORE=@
$SUBROUTINE ADVAN6 TOL=5
$MODEL      NCOMP=1 COMP=(CENTRAL,DEFDOSE)
$PK       
CL = THETA(1)*EXP(ETA(1))
V  = THETA(2)
KE = CL/V

$DES      
DADT(1)=-A(1)*KE

$ERROR         
IPRED=A(1)/V    ; predicting concentration
WA=THETA(3)    ; still including a theta for additive error. we'll fix it to zero in the THETA block
WP=THETA(4)
W = SQRT(WA**2+(WP*IPRED)**2) 
IRES=DV-IPRED
IWRES=IRES/W
Y = IPRED + W*EPS(1)            

$THETA  (0,0.00904) ; 1 CL L/h (based on excel estimates)
 (0,0.0669) ; 2 V L (based on excel estimates)
 0 FIX ; 3 Additive Error. Do (0)FIX for a proportional error model.
 (0,0.168) ; 4 Proportional Error. Do (0)FIX for an additive error model.
$OMEGA  0  FIX  ;     CL IIV
$SIGMA  1  FIX  ; eps residual variance
$ESTIMATION METHOD=1 INTERACTION PRINT=5 SIG=2 MAXEVAL=9999 NOABORT
$COVARIANCE UNCONDITIONAL
$TABLE      ID TIME IPRED PRED DV AMT WRES IWRES CWRES EVID
            FILE=sdtab04 NOPRINT ONEHEADER

