(* usage: NBDIR=<math/subtraction-notebooks> wolframscript -file verify_W_gq_fortran.wls
   Soft-photon term of gq s_12 for general charges (calG_gq_cc, mod_int_sub_nnlo.f90), written for two charged leptons
   (charges q3, q4, Q1 + Q2 = q3 + q4) and reduced to Z charges, against the single boost of Z.nb (FINITEgq).
   The lepton-lepton function X is not in the W code: it is solved for, and must depend on E3, E4, eta34, L only.
   Not a derivation: terms proportional to (Q1 + Q2) that do not depend on the lepton kinematics vanish for Z. *)
nb=Get[Environment["NBDIR"]<>"/Z.nb"];
cells=Cases[nb,Cell[BoxData[b_],"Input",___]:>b,Infinity];
Do[If[!FreeQ[c,"Ptildegq"|"TCgqBoost"|"Gq"]&&!FreeQ[c,":="]&&FreeQ[c,"FINITEgq"|"FBFortranForm"],ToExpression[c,StandardForm]],{c,cells}];
(* Z.nb, FINITEgqFB[2]: bracket multiplying Ptildegq *)
zb = Gq[Qq,Qe] - 2 Qe Qq Log[z] Log[(E3 eta[1,3])/(E4 eta[1,4])] - 4 Zeta[2] Qq^2;
(* genGnloQCD_cc + legfinalcalG_cc + calG_gq_cc, lepton l, ileg = 1; lmu = log(mu^2/4EC^2), s_1l/s12 = El eta_1l/EC *)
lep[l_,q_] := Module[{El=ToExpression["E"<>ToString[l]],lE,le1,le2,li1,li2},
  lE=Log[EC/El]; le1=Log[eta[1,l]]; le2=Log[eta[2,l]]; li1=PolyLog[2,1-eta[1,l]]; li2=PolyLog[2,1-eta[2,l]];
  q^2(13/2+3lE+2lE^2-2Pi^2/3+(3/2+2lE)L) + q Q1(-3lE-lE^2+3le1+2lE le1+2li1-Pi^2/3-(3+2lE)L) + q Q2(-3lE-lE^2+3le2+2lE le2+2li2-Pi^2/3-(3+2lE)L)
  - 2 Log[z] q Q1 Log[El eta[1,l]/EC]];
quark = Q1 Q2(-2Pi^2/3+3L) + (3L-4Zeta[2])Q1^2;
num = Join[{E3->37/10,E4->53/10,EC->71/10,L->13/10,z->3/10},Flatten@Table[eta[i,j]->Prime[3 i+j]/50,{i,1,3},{j,i+1,4}]];
Do[
  d = Expand[PowerExpand[zb - (lep[3,sg Qe]+lep[4,-sg Qe]+quark /. {Q1->Qq,Q2->-Qq})]];
  Print["q3 = ",sg," Qe:  Qq^2 terms ",N[Coefficient[d,Qq,2]/.num,20],"   Qq Qe terms ",Chop[N[Coefficient[d,Qq,1]/.Qe->1/.num,20],10^-15]];
  X = -Simplify[Coefficient[d,Qq,0]/Qe^2];
  Print["   lepton-lepton function X depends on ",Select[{z,L,E3,E4,EC,eta[1,3],eta[1,4],eta[2,3],eta[2,4],eta[3,4]},!FreeQ[X,#]&]];
  Print["   X = ",X],
 {sg,{1,-1}}]
