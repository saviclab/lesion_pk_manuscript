;; Based on: 130
;; Description: two compartment, metab, prop error SAME for parent + M1, nonlinear CL, IIV CLINT CLm, diff F for 1800 mg dose, M3 method
;; Author: Anu

$PROBLEM  SUTEZOLID population PK

$INPUT    ID TIME DV AMT EVID MDV CMT SEX AGE TRTPN MFL 
          RACEN ETHNICN WEIGHTSC BMISC BLOQ_new STRT

$DATA     SUTPK_COV_20211101_VPC.csv ; new dataset for VPC stratification
          IGNORE=@
          IGNORE(CMT.GT.4) ; only looking at parent + M1

; NOTES about data:
; AMT in mg
; DV in mg/L
; TIME in h

$SUBROUTINE   ADVAN6 TOL=6

$MODEL        NCOMP=4
              COMP(DEPOT)           ; dosing
              COMP(CENTRAL,DEFOBS)  ; parent
              COMP(PERIPH)          ; since 2 cmpt model
              COMP(METAB)           ; major metabolite

$PK
; Nonlinear CL, L/h
CLINT = THETA(1)*EXP(ETA(1))
KM = THETA(10)*EXP(ETA(6))
VMAX = KM*CLINT ; at low concs, assume CL = CLint

; Central Volume, L
TVV = THETA(2)
V = TVV*EXP(ETA(2))
S2 = V

; Peripheral Volume, L
TVV2 = THETA(6)
V2 = TVV2

; Bioavailability
IF(TRTPN.LT.4) F1 = 1
IF(TRTPN.EQ.4) F1 = THETA(11)

; Absorption rate, 1/h
TVKA = THETA(3)
KA = TVKA*EXP(ETA(3))

; Fraction of drug metabolized
Fm = 1 ; assume this to be 1

; Intercompartment CL, L/h
TVQ = THETA(7)
Q = TVQ

; Metabolite CL, L/h
TVCLm = THETA(8)
CLm = TVCLm*EXP(ETA(4))

; Metabolite Volume, L/h
TVVm = THETA(9)
Vm = TVVm*EXP(ETA(5))

$DES
CP = A(2)/S2
DADT(1) = -KA*A(1)
DADT(2) = KA*A(1) - Fm*(VMAX*CP)/(KM+CP) + Q/V2*A(3) - Q/V*A(2)
DADT(3) = Q/V*A(2) - Q/V2*A(3)
DADT(4) = Fm*(VMAX*CP)/(KM+CP) - (CLm/Vm)*A(4)

$ERROR
IPRED = F
IF(F.GT.0.AND.CMT.EQ.2) IPRED = A(2)/V
IF(F.GT.0.AND.CMT.EQ.4) IPRED = A(4)/Vm

PROP = THETA(4) ; Proportional error
ADD = THETA(5)  ; Additive error

W = SQRT(PROP**2*IPRED**2+ADD**2) ; combined parent and metabolite error model

IRES = DV - IPRED
IWRES = IRES/W

; IF(CMT.EQ.2) Y = IPRED + W*EPS(1)
; IF(CMT.EQ.4) Y = IPRED + W*EPS(2) ; take metab from different sigma distribution

; M3 method
LLOQ = 0.01 ; parent + M1
DUM = (LLOQ-IPRED)/W
CUMD = PHI(DUM)
;-------Above limit of quantification - parent-------
  IF(DV.GT.LLOQ.AND.CMT.EQ.2) THEN
  F_FLAG = 0
  Y = IPRED + W*EPS(1)
  ENDIF
;-------Above limit of quantification - M1-------
  IF(DV.GT.LLOQ.AND.CMT.EQ.4) THEN
  F_FLAG = 0
  Y = IPRED + W*EPS(2)
  ENDIF
;-------Below limit of quantification-------
  IF(DV.LE.LLOQ) THEN
  F_FLAG = 1
  Y = CUMD
  ENDIF

$THETA  
 (0,200)  ; CLINT
 (0,500)  ; Vc (V2)
 (0,1)    ; KA
 (0,0.1)  ; PROP
 0 FIX    ; ADD
 (0,500)  ; Vp (V3)
 (0,100)  ; Q
 (0,50)   ; CLm
 (0,100)  ; Vm
 (0,4)    ; KM
 (0,0.8)  ; F1 (1800 mg)
 ; KM is Cmax for the dose that begins to indicate saturation (1800 mg)

$OMEGA  
 0.2      ; IIV CLINT
 0 FIX    ; IIV Vc
 0 FIX    ; IIV KA
 0.2      ; IIV CLm
 0 FIX    ; IIV Vm
 0 FIX    ; IIV KM

$SIGMA  
 1 FIX    ; Residual variability
 1 FIX

$ESTIMATION METH=1 INTERACTION LAPLACIAN MAXEVAL=9999 NITER=100 NSIG=3 NOABORT PRINT=10
$COVARIANCE UNCONDITIONAL PRINT=E
$TABLE      ID TIME DV AMT EVID MDV CMT SEX AGE TRTPN MFL 
            RACEN ETHNICN WEIGHTSC BMISC BLOQ_new STRT
			      CLINT V V2 Q KA CLm Vm KM ETA1 ETA4 IPRED PRED
			      CWRES IWRES ONEHEADER NOPRINT FILE=tab135
