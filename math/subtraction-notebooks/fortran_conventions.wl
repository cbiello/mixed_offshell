(* ::Package:: *)

(* ::Title:: *)
(*Fortran conventions for the subtraction notebooks (PLAN M2)*)


(* ::Section:: *)
(*Conventions*)


(* Legs as in the Fortran: q(1) q'(2) -> l(3) l'(4), momenta p1..p4, leptons in lep_ids order.
   Charges: Qq = Q_IS(al), Qqp = Q_IS(be) (incoming, physical); Q3, Q4 outgoing; Qq + Qqp = Q3 + Q4.
   The charged lepton of a W is on leg ilept (W+: 4, W-: 3) with charge Qlept (Qlept_chg).
   Kinematics: E1 = E2 = EC, E3, E4, eta[i,j] (i<j), s_ij = 4 Ei Ej eta[i,j]; lmu = Log[mu^2/(4 EC^2)].
   Structures: FLM (Born), FLVqcdfin, FLVewfin, FVVqcdewfin; FLM[z p1,p2,p3,p4,z] = leg 1 boosted. *)
lepCharges["wp"] = {Q3 -> 0, Q4 -> Qlept};
lepCharges["wm"] = {Q3 -> Qlept, Q4 -> 0};
lepCharges["z"]  = {Q3 -> -Qlept, Q4 -> Qlept, Qqp -> -Qq};
ilept["wp"] = 4; ilept["wm"] = 3;
El[il_] := Symbol["E" <> ToString[il]];


(* ::Section:: *)
(*Maps from the other notebooks*)


(* W.nb: q qbp -> nu l, charges Qq, Qqbp outgoing (= -Q_IS), Qe of l; l is leg 4. For ilept = 3 the lepton moves to leg 3. *)
fromWnb[il_] := Join[
  {Qq -> -Qq, Qqbp -> -Qqp, Qe -> Qlept},
  {(h:FLM|FLVqcdfin|FLVewfin|FVVqcdewfin)[a_, b_, nu, l, r___] :> h[a /. q -> p1, b /. qbp -> p2, p3, p4, r]},
  If[il == 3, {E4 -> E3, eta[i_, 4] :> eta[i, 3]}, {}]];
(* notes/fabrizio/subtraction_f.wl: all-outgoing charges q1..q4, q1 + q2 + q3 + q4 = 0 *)
fromAllOut = {q1 -> -Qq, q2 -> -Qqp, q3 -> Q3, q4 -> Q4};


(* ::Section:: *)
(*Kernels (mod_hoppet_nnlo.f90, mod_splittings_bare.f90)*)


z2r = Pi^2/6;
Kns[z_, Lh_] := (-2(1+z)Log[1-z] - (1+z^2)/(1-z)Log[z] + (1-z)) + 4 D1[z] + 2 z2r delta[1-z] +
  Lh(-(1+z^2)/(1-z) + 2/(1-z) - 2 D0[z] - 3/2 delta[1-z]);                 (* mySub_ns_func, _Lmu *)
P[z_]   := (1-z) - 2(1+z)Log[1-z] + (1+z-2/(1-z))Log[z] + 4 D1[z];          (* PqqNLO_func *)
Pl[z_]  := (1+z) - 2 D0[z];                                                 (* PqqNLO_Lmu_func *)
Pb[z_]  := ((1-z) - 2(1+z)Log[1-z] + (1+z-2/(1-z))Log[z] + 4 Log[1-z]/(1-z))Log[z];   (* PqqNLO_b_func *)
Pbl[z_] := (1+z-2/(1-z))Log[z];                                             (* PqqNLO_b_Lmu_func *)
PA[z_]  := (1-z) - 2(1+z)Log[1-z] + 4 D1[z];                                (* Pqq_NLO *)
PlA[z_] := (1+z) - 2 D0[z];                                                 (* PqqAP_0_R *)
TC0[x_] := Module[{omx=1-x,opx=1+x,lx=Log[x],lomx=Log[1-x],li2x=PolyLog[2,x],li3x=PolyLog[3,x],li3xb=PolyLog[3,1-x],zeta2=z2r,zeta3=Zeta[3]},
  (-16*lomx^2*lx)/omx-(4*lx^3)/omx+12*li3x*opx - 6*li3xb*opx-8*lomx^3*opx-(13*lx^3*opx)/6 + li2x*(-6*omx-2*lomx*opx-4*lx*opx)+lomx^2*(-4*omx+lx*opx)+(4*lx*(-3*lx+4*lomx*lx + 2*(li2x-zeta2)))/omx+6*(3*omx+(5*lomx*opx)/3-(8*lx*opx)/3)*zeta2-28*opx*zeta3 + (12*lx^2*((-4*lx)/9+1))/omx+lx^2*(-x+2)+lx^2*(-4*lomx*opx+3*lx*opx-4*(x + 2))+lomx*(4*lx^2*opx-8*lx*(x+2)-3*(3*x-4))+lx*(-9*x+5) + lx*((-24*lomx)/omx-(24*lomx^2)/omx+12*lomx^2*opx-lx^2*opx+(8*zeta2)/omx-4*opx*zeta2 + 4*lx*(x+2)+lomx*(-8*lx*opx+8*(x+2))+2*(5*x-6))+2*(8*x-9) + (10*lx^3-4*li3xb+12*lomx*(2*lx-li2x+zeta2)+8*lx*(li2x+(3*lomx^2)/4 + 4*zeta2)-32*(li3x-zeta3))/omx
  + 16(D3[x]-zeta2 D1[x]+2 zeta3 D0[x]) - 2 Pi^4/45 delta[1-x]];          (* mySub_TCqq_func *)
TC1[x_] := Module[{omx=1-x,opx=1+x,lx=Log[x],lomx=Log[1-x],li2x=PolyLog[2,x],zeta2=z2r,zeta3=Zeta[3]},
  4*(-3*lx+4*lomx*lx+2*(li2x-zeta2))/omx - 8*lx^2/omx + 2*(-4*(2+x) + 12/omx + 8*lomx/omx + (3*lx-4*lomx)*opx)*lx - 4*zeta2*opx + 2*(-6 + 5*x) + 12*opx*lomx^2 + 4*(2 + x)*lx - opx*lx^2 + 8*lomx*(2+x-opx*lx)
  + 8(zeta2 D0[x]-3 D1[x]-3 D2[x]) - 2(Pi^2+8 zeta3) delta[1-x]];         (* mySub_TCqq_Lmu_func *)
TC2[x_] := Module[{omx=1-x,opx=1+x,lx=Log[x],lomx=Log[1-x]},
  -4*(lx/omx + 2 + x + opx*lomx) + 3*opx*lx + 12 D0[x] + 8 D1[x] + (9/2-4 z2r) delta[1-x]];   (* mySub_TCqq_Lmu2_func *)
TC[z_, Lh_] := TC0[z] + Lh TC1[z] + Lh^2 TC2[z];


(* ::Section:: *)
(*Elastic functions (mod_int_sub_nnlo.f90: fin_ns_vqcd_cc, genGnloQCD_cc, legfinalcalG_cc, fin_ns_v7_cc), charged lepton on leg il*)


lnks[il_] := {Log[eta[1,il]], Log[eta[2,il]], PolyLog[2,1-eta[1,il]], PolyLog[2,1-eta[2,il]]};
finNsVqcdCC[il_, lmu_] := Module[{l1,l2,li1,li2}, {l1,l2,li1,li2} = lnks[il];
  Qlept^2(13/2-2Pi^2/3+3/2 lmu) + Qlept Qq(-Pi^2/3+3l1+2li1-3lmu) + Qlept Qqp(-Pi^2/3+3l2+2li2-3lmu)
  + Qq Qqp(-2Pi^2/3+3lmu) + (Qq^2+Qqp^2)(-2 z2r+3/2 lmu)];
genGnloQCDcc[il_, L_] := Module[{lE=Log[EC/El[il]],l1,l2,li1,li2}, {l1,l2,li1,li2} = lnks[il];
  Qlept^2(13/2+3lE+2lE^2-2Pi^2/3+(3/2+2lE)L) + Qlept Qq(-3lE-lE^2+3l1+2lE l1+2li1-Pi^2/3-(3+2lE)L)
  + Qlept Qqp(-3lE-lE^2+3l2+2lE l2+2li2-Pi^2/3-(3+2lE)L) + Qq Qqp(-2Pi^2/3+3L)];
legfinalcalGcc[il_, L_, ileg_] := genGnloQCDcc[il, L] + 3 L If[ileg == 1, Qq, Qqp]^2;
finNsV7cc[il_, L_] := Module[{Qd=Qlept-Qq-Qqp,EB,dC},
  EB = -2/45(Qq^2+Qqp^2)Pi^4 + (-6Qd^2-16(Qq^2+Qqp^2)L)Zeta[3] + 1/4(-9Qlept^2+(-9-8Pi^2/3)Qq^2-18Qq Qqp-(27+8Pi^2)/3 Qqp^2+18Qlept(Qq+Qqp))L^2 - Pi^2 Qd^2(-21/8-L);
  dC = (Qq^2+Qqp^2)(-2Pi^4/45-2(Pi^2+8Zeta[3])L+(9/2-4 z2r)L^2);
  (2Pi^2/3-3L) genGnloQCDcc[il, L] + EB - dC];


(* ::Section:: *)
(*q qb' finite remainder per sector (mod_xsects_nnlo_s_ns.f90)*)


(* coefficient of alpha/2pi alpha_s/2pi; the delta(1-z) terms act on the boosted structures (applyd) *)
finW[il_] := Module[{lmu=Log[mu^2/(4 EC^2)], Lh1=Log[mu^2/(4 EC^2 z)], Lh12=Log[mu^2/(4 EC^2 z1 z2)],
    F0=FLM[p1,p2,p3,p4], F1=FLM[p1 z,p2,p3,p4,z], F2=FLM[p1,p2 z,p3,p4,z], F12=FLM[p1 z1,p2 z2,p3,p4,z1,z2],
    Vq0=FLVqcdfin[p1,p2,p3,p4], Vq1=FLVqcdfin[p1 z,p2,p3,p4,z], Vq2=FLVqcdfin[p1,p2 z,p3,p4,z],
    Ve1=FLVewfin[p1 z,p2,p3,p4,z], Ve2=FLVewfin[p1,p2 z,p3,p4,z], lb=Log[El[il] eta[1,il]/EC], lbp=Log[El[il] eta[2,il]/EC]},
  <|"s_vqcd" -> Vq0 finNsVqcdCC[il,lmu] + Qq^2 Kns[z,Lh1] Vq1 + Qqp^2 Kns[z,Lh1] Vq2,
    "s_vewk" -> CF Kns[z,Lh1](Ve1 + Ve2),
    "s_12"   -> CF((Qq^2+Qqp^2) F12 ((P[z1]+Lh12 Pl[z1])(P[z2]+Lh12 Pl[z2]) + (Pb[z1]+Lh12 Pbl[z1])Pl[z2] + Pl[z1](Pb[z2]+Lh12 Pbl[z2]) + Pbl[z1]Pbl[z2])
      - 2 Qlept (Pb[z]+Lh1 Pbl[z])(Qq lb F1 + Qqp lbp F2) + TC[z,Lh1](Qq^2 F1 + Qqp^2 F2)
      + (PA[z]+lmu PlA[z])(legfinalcalGcc[il,lmu,1] F1 + legfinalcalGcc[il,lmu,2] F2) + finNsV7cc[il,lmu] F0)|>];
stru = FLM|FLVqcdfin|FLVewfin;
applyd[e_] := Expand[e] /. delta[1-z] x_ :> (x /. {(h:stru)[a_ z, b_, c_, d_, z] :> h[a, b, c, d], (h:stru)[a_, b_ z, c_, d_, z] :> h[a, b, c, d]} /. z -> 1);
