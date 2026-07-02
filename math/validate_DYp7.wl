(* ::Package:: *)

(* =================================================================== *)
(*  Validation of the Fortran FINITEDYp[7]  (fin_ns_v7_cc)             *)
(*  against the Mathematica  FINITEDYp[7]  of DY+_qq.wl.               *)
(*                                                                     *)
(*  WHAT IT CHECKS                                                     *)
(*   The Fortran ships the elastic soft [7] as                        *)
(*        respdf_4 = (4 zeta2 - 3 L) genGnloQCD + EBcode - dCcode      *)
(*   while term C (respdf_2) already ships the per-leg delta piece     *)
(*        dCcode  (the mySub_TCqq / _Lmu / _Lmu2 delta, x Q_IS^2).     *)
(*   Hence the TOTAL code contribution attributed to [7] is            *)
(*        codeV7 = respdf_4 + dCcode = (4 zeta2 - 3 L) genGnloQCD      *)
(*                 + EBcode                                            *)
(*   which must equal the notebook FINITEDYp[7] (up to CF*FLM).        *)
(*   => Check 1 isolates the explicit-bracket transcription EBcode.    *)
(*                                                                     *)
(*  Check 2 (optional) additionally tests that the Fortran soft        *)
(*  function calG_ONLOQCD_ns equals the notebook genGnloQCD.           *)
(*                                                                     *)
(*  HOW TO RUN                                                         *)
(*   1. Evaluate all cells of DY+_qq.wl first, so that genGnloQCD and  *)
(*      FINITEDYp[7] are defined in the current kernel session.        *)
(*   2. SetDirectory to this folder and evaluate:                     *)
(*         Get["validate_DYp7.wl"]                                     *)
(* =================================================================== *)

(* scale log used everywhere; matches Log[mu^2/(4 EC^2)] in the notebook *)
LL = Log[\[Mu]^2/(4 EC^2)];

(* ------------------------------------------------------------------- *)
(* Fortran transcription: explicit charge bracket EB (fin_ns_v7_cc).   *)
(*   charges:  Qq -> Q_IS(al),  Qqp -> Q_IS(be),  Qlp -> Qlept          *)
(*   Qd = Qlp - Qq - Qqp                                                *)
(* ------------------------------------------------------------------- *)
Qd = Qlp - Qq - Qqp;
EBcode = (
      -2/45 (Qq^2 + Qqp^2) Pi^4
    + (-6 Qd^2 - 16 (Qq^2 + Qqp^2) LL) Zeta[3]
    + 1/4 (-9 Qlp^2 + (-9 - 8 Pi^2/3) Qq^2 - 18 Qq Qqp
           - (27 + 8 Pi^2)/3 Qqp^2 + 18 Qlp (Qq + Qqp)) LL^2
    - Pi^2 Qd^2 (-21/8 - LL)
  );

(* ------------------------------------------------------------------- *)
(* Fortran transcription: term-C delta already shipped in respdf_2.    *)
(*   dCcode = (Q_IS(al)^2 + Q_IS(be)^2) * (per-leg xPij deltas)         *)
(*   with  mySub_TCqq      delta = -2 Pi^4/45                           *)
(*         mySub_TCqq_Lmu  delta = -2 (Pi^2 + 8 Zeta[3])   (coeff of L) *)
(*         mySub_TCqq_Lmu2 delta =  9/2 - 4 zeta2          (coeff of L^2)*)
(*   zeta2 = Pi^2/6                                                     *)
(* ------------------------------------------------------------------- *)
dCcode = (Qq^2 + Qqp^2) (
      -2 Pi^4/45
      - 2 (Pi^2 + 8 Zeta[3]) LL
      + (9/2 - 4 (Pi^2/6)) LL^2 );

(* ------------------------------------------------------------------- *)
(* Total code contribution to [7]:  respdf_4 + (term-C delta)          *)
(*   respdf_4       = (4 zeta2 - 3 L) genGnloQCD + EBcode - dCcode      *)
(*   term-C delta   = dCcode                                           *)
(*   => codeV7      = (4 zeta2 - 3 L) genGnloQCD + EBcode              *)
(* uses the notebook's OWN genGnloQCD, so the genG parts cancel in the *)
(* difference below and only the EB transcription is tested.           *)
(* ------------------------------------------------------------------- *)
codeV7 = (4 Zeta[2] - 3 LL) genGnloQCD + EBcode;

(* notebook reference, stripped of the common CF*FLM prefactor *)
nbV7 = FINITEDYp[7] /. {CF -> 1, FLM[___] -> 1};

Print["==================================================================="];
Print["Check 1:  codeV7 - FINITEDYp[7]   (expect 0)"];
diff1 = Simplify[codeV7 - nbV7];
Print["   -> ", diff1];
If[PossibleZeroQ[diff1] || diff1 === 0,
   Print["   PASS: explicit bracket EBcode matches the notebook."],
   Print["   MISMATCH. Residual (grouped by charge structure):"];
   Print["     coeff Qq^2   : ", Simplify[Coefficient[diff1 // Expand, Qq^2]]];
   Print["     coeff Qqp^2  : ", Simplify[Coefficient[diff1 // Expand, Qqp^2]]];
   Print["     coeff Qlp^2  : ", Simplify[Coefficient[diff1 // Expand, Qlp^2]]];
   Print["     coeff Qq Qqp : ", Simplify[Coefficient[diff1 // Expand, Qq Qqp]]];
   Print["     coeff Qlp Qq : ", Simplify[Coefficient[diff1 // Expand, Qlp Qq]]];
   Print["     coeff Qlp Qqp: ", Simplify[Coefficient[diff1 // Expand, Qlp Qqp]]];
];
Print["==================================================================="];

(* ------------------------------------------------------------------- *)
(* Check 2: the Fortran genGnloQCD (computed directly in fin_ns_v7_cc)  *)
(* against the notebook genGnloQCD, in the code variables               *)
(*    lE  = Log[EC/E4]  (charged lepton = leg 4)                        *)
(*    LL  = Log[mu^2/4EC^2]                                             *)
(*    le14=Log[eta14], le24=Log[eta24], li14=Li2(1-eta14), li24=...     *)
(* Notebook genGnloQCD is written in twoEConMu, twoE4onMu, eta[1,2];    *)
(* map with the notebook's own leg-replacement rules:                   *)
(*    Log[twoEConMu] -> -LL/2                                           *)
(*    Log[twoE4onMu] -> -lE - LL/2                                      *)
(*    eta[1,2] -> 1                                                     *)
(* (uses opaque symbols so nothing gets over-simplified)                *)
(* ------------------------------------------------------------------- *)
cQlp2   = 13/2 + 3 lE + 2 lE^2 - 2 Pi^2/3 + (3/2 + 2 lE) LL;
cQlpQq  = -3 lE - lE^2 + 3 le14 + 2 lE le14 + 2 li14 - Pi^2/3 - (3 + 2 lE) LL;
cQlpQqp = -3 lE - lE^2 + 3 le24 + 2 lE le24 + 2 li24 - Pi^2/3 - (3 + 2 lE) LL;
cQqQqp  = -2 Pi^2/3 + 3 LL;
codeGenG = Qlp^2 cQlp2 + Qlp Qq cQlpQq + Qlp Qqp cQlpQqp + Qq Qqp cQqQqp;

genGnb = genGnloQCD /. {
     Log[twoEConMu] -> -LL/2,
     Log[twoE4onMu] -> -lE - LL/2,
     eta[1, 2] -> 1,
     Log[eta[1, 4]] -> le14, Log[eta[2, 4]] -> le24,
     PolyLog[2, 1 - eta[1, 4]] -> li14, PolyLog[2, 1 - eta[2, 4]] -> li24};

Print["Check 2:  fin_ns_v7_cc genG - notebook genGnloQCD (mapped)  (expect 0)"];
diff2 = Simplify[codeGenG - genGnb];
Print["   -> ", diff2];
If[PossibleZeroQ[diff2] || diff2 === 0,
   Print["   PASS: Fortran genG == notebook genGnloQCD."],
   Print["   MISMATCH: residual = ", diff2]];
Print["==================================================================="];

(* ------------------------------------------------------------------- *)
(* Check 3/4: Fortran legfinalcalG_cc (used in [6]/[60] term A and the  *)
(* plus-subtraction) vs the notebook leg1finalcalG / leg2finalcalG.     *)
(*   legfinalcalG_cc(ileg=1) = genGnloQCD + 3 Q_IS(al)^2 * LL           *)
(*                           = genGnloQCD - 3 Qq^2 Log[4EC^2/mu^2].     *)
(* Map the notebook symbols (Log[fourEC2overmu2], Log[twoEConMu],        *)
(* Log[E4overEC], eta) to the code variables.                           *)
(* ------------------------------------------------------------------- *)
legmap = {
   Log[twoEConMu] -> -LL/2,
   Log[fourEC2overmu2] -> -LL,
   Log[E4overEC] -> -lE,
   eta[1, 2] -> 1,
   Log[eta[1, 4]] -> le14, Log[eta[2, 4]] -> le24,
   PolyLog[2, 1 - eta[1, 4]] -> li14, PolyLog[2, 1 - eta[2, 4]] -> li24};

codeLeg1 = codeGenG + 3 Qq^2 LL;
codeLeg2 = codeGenG + 3 Qqp^2 LL;
nbLeg1 = leg1finalcalG /. legmap;
nbLeg2 = leg2finalcalG /. legmap;

Print["Check 3:  legfinalcalG_cc(1) - notebook leg1finalcalG  (expect 0)"];
diff3 = Simplify[codeLeg1 - nbLeg1];
Print["   -> ", diff3];
If[PossibleZeroQ[diff3] || diff3 === 0,
   Print["   PASS."], Print["   MISMATCH: ", diff3]];

Print["Check 4:  legfinalcalG_cc(2) - notebook leg2finalcalG  (expect 0)"];
diff4 = Simplify[codeLeg2 - nbLeg2];
Print["   -> ", diff4];
If[PossibleZeroQ[diff4] || diff4 === 0,
   Print["   PASS."], Print["   MISMATCH: ", diff4]];
Print["==================================================================="];
