(* usage: NBOUT=<dir with wnb.mx from run_W_nb.wls> wolframscript -file verify_Wp_qq_fortran.wls   -> every line must be 0 (elastic lines: ~1e-30) *)
(* W+ q qb' channel: Fortran (s_12, s_vqcd, s_vewk) transcribed, against the finite remainder of W.nb *)
Get[Environment["NBOUT"]<>"/wnb.mx"];
fin = finitesimpl /. {Qqbp->QqpOLD,Qq->QqOLD,Qe->QeOLD,l->lOLD} /. {QqpOLD->-Qqp,QqOLD->-Qq,QeOLD->Qlp,qbp->qp,lOLD->lp,nu->l};
pref = \[Alpha]ontwopi*\[Alpha]sontwopi;
com = {E1->EC,E2->EC,Emax->EC,eta[1,2]->1};
z2r = Pi^2/6;

(* ---------- Fortran transcription ---------- *)
(* mod_hoppet_nnlo.f90: kernels as distributions; Lh = log(mu^2/mll^2) *)
Kns[z_,Lh_] := (-2(1+z)Log[1-z]-(1+z^2)/(1-z)Log[z]+(1-z)) + 4 D1[z] + 2 z2r delta[1-z] + Lh(-(1+z^2)/(1-z) + 2/(1-z) - 2 D0[z] - 3/2 delta[1-z]);
P[z_]   := (1-z)-2(1+z)Log[1-z]+(1+z-2/(1-z))Log[z] + 4 D1[z];      (* PqqNLO_func *)
Pl[z_]  := (1+z) - 2 D0[z];                                          (* PqqNLO_Lmu_func *)
Pb[z_]  := ((1-z)-2(1+z)Log[1-z]+(1+z-2/(1-z))Log[z]+4Log[1-z]/(1-z))Log[z];  (* PqqNLO_b_func *)
Pbl[z_] := (1+z-2/(1-z))Log[z];                                      (* PqqNLO_b_Lmu_func *)
TC0[x_] := Module[{omx=1-x,opx=1+x,lx=Log[x],lomx=Log[1-x],li2x=PolyLog[2,x],li3x=PolyLog[3,x],li3xb=PolyLog[3,1-x],zeta2=z2r,zeta3=Zeta[3]},
  (-16*lomx^2*lx)/omx-(4*lx^3)/omx+12*li3x*opx - 6*li3xb*opx-8*lomx^3*opx-(13*lx^3*opx)/6 + li2x*(-6*omx-2*lomx*opx-4*lx*opx)+lomx^2*(-4*omx+lx*opx)+(4*lx*(-3*lx+4*lomx*lx + 2*(li2x-zeta2)))/omx+6*(3*omx+(5*lomx*opx)/3-(8*lx*opx)/3)*zeta2-28*opx*zeta3 + (12*lx^2*((-4*lx)/9+1))/omx+lx^2*(-x+2)+lx^2*(-4*lomx*opx+3*lx*opx-4*(x + 2))+lomx*(4*lx^2*opx-8*lx*(x+2)-3*(3*x-4))+lx*(-9*x+5) + lx*((-24*lomx)/omx-(24*lomx^2)/omx+12*lomx^2*opx-lx^2*opx+(8*zeta2)/omx-4*opx*zeta2 + 4*lx*(x+2)+lomx*(-8*lx*opx+8*(x+2))+2*(5*x-6))+2*(8*x-9) + (10*lx^3-4*li3xb+12*lomx*(2*lx-li2x+zeta2)+8*lx*(li2x+(3*lomx^2)/4 + 4*zeta2)-32*(li3x-zeta3))/omx
  + 16(D3[x]-zeta2 D1[x]+2 zeta3 D0[x]) - 2 Pi^4/45 delta[1-x]];
TC1[x_] := Module[{omx=1-x,opx=1+x,lx=Log[x],lomx=Log[1-x],li2x=PolyLog[2,x],zeta2=z2r,zeta3=Zeta[3],two=2},
  4*(-3*lx+4*lomx*lx+2*(li2x-zeta2))/omx -8*lx^2/omx + 2*(-4*(two+x) + 12/omx + 8*lomx/omx + (3*lx-4*lomx)*opx)*lx - 4*zeta2*opx + 2*(-6 + 5*x) + 12*opx*lomx^2 + 4*(two + x)*lx - opx*lx^2 + 8*lomx*(two+x-opx*lx)
  + 8(zeta2 D0[x]-3 D1[x]-3 D2[x]) - 2(Pi^2+8 zeta3) delta[1-x]];
TC2[x_] := Module[{omx=1-x,opx=1+x,lx=Log[x],lomx=Log[1-x],two=2},
  -4*(lx/omx + two + x + opx*lomx) + 3*opx*lx + 12 D0[x] + 8 D1[x] + (9/2-4 z2r) delta[1-x]];
TC[z_,Lh_] := TC0[z] + Lh TC1[z] + Lh^2 TC2[z];
(* mod_splittings_bare.f90: Pqq_NLO, PqqAP_0_R (explicit z, plus prescription by endpoint subtraction) *)
PA[z_]  := (1-z) - 2(1+z)Log[1-z] + 4 D1[z];
PlA[z_] := (1+z) - 2 D0[z];
(* mod_int_sub_nnlo.f90; lmu = log(mu^2/4EC^2), lE = log(EC/E4), Q1 = Q_IS(al), Q2 = Q_IS(be), Ql = Qlept *)
finvqcd[Q1_,Q2_,Ql_,lmu_] := Module[{l12=0,li12=0,l14=Log[eta[1,4]],li14=PolyLog[2,1-eta[1,4]],l24=Log[eta[2,4]],li24=PolyLog[2,1-eta[2,4]]},
  Ql^2(13/2-2Pi^2/3+3/2 lmu) + Ql Q1(-Pi^2/3+3l14+2li14-3lmu) + Ql Q2(-Pi^2/3+3l24+2li24-3lmu) + Q1 Q2(-2Pi^2/3-3l12-2li12+3lmu) + (Q1^2+Q2^2)(-2 z2r+3/2 lmu)];
genG[Q1_,Q2_,Ql_,L_] := Module[{lE=Log[EC/E4],le1=Log[eta[1,4]],le2=Log[eta[2,4]],li1=PolyLog[2,1-eta[1,4]],li2=PolyLog[2,1-eta[2,4]]},
  Ql^2(13/2+3lE+2lE^2-2Pi^2/3+(3/2+2lE)L) + Ql Q1(-3lE-lE^2+3le1+2lE le1+2li1-Pi^2/3-(3+2lE)L) + Ql Q2(-3lE-lE^2+3le2+2lE le2+2li2-Pi^2/3-(3+2lE)L) + Q1 Q2(-2Pi^2/3+3L)];
legG[Q1_,Q2_,Ql_,L_,i_] := genG[Q1,Q2,Ql,L] + 3 L If[i==1,Q1,Q2]^2;
finv7[Q1_,Q2_,Ql_,L_] := Module[{Qd=Ql-Q1-Q2,EB,dC},
  EB = -2/45(Q1^2+Q2^2)Pi^4 + (-6Qd^2-16(Q1^2+Q2^2)L)Zeta[3] + 1/4(-9Ql^2+(-9-8Pi^2/3)Q1^2-18Q1 Q2-(27+8Pi^2)/3 Q2^2+18Ql(Q1+Q2))L^2 - Pi^2 Qd^2(-21/8-L);
  dC = (Q1^2+Q2^2)(-2Pi^4/45-2(Pi^2+8Zeta[3])L+(9/2-4 z2r)L^2);
  (2Pi^2/3-3L)genG[Q1,Q2,Ql,L] + EB - dC];

(* mod_xsects_nnlo_s_ns.f90, Vcharge == 1 branches *)
lmu = Log[mu^2/(4 EC^2)];                       (* fill_sv_logs(muf^2, 4 EC^2) *)
Lh1 = Log[mu^2/(4 EC^2 z)]; Lh12 = Log[mu^2/(4 EC^2 z1 z2)];   (* hoppet log, mu2ref = mll^2 *)
F0=FLM[q,qp,l,lp]; F1=FLM[q z,qp,l,lp,z]; F2=FLM[q,qp z,l,lp,z]; F12=FLM[q z1,qp z2,l,lp,z1,z2];
Vq0=FLVqcdfin[q,qp,l,lp]; Vq1=FLVqcdfin[q z,qp,l,lp,z]; Vq2=FLVqcdfin[q,qp z,l,lp,z];
Ve0=FLVewfin[q,qp,l,lp]; Ve1=FLVewfin[q z,qp,l,lp,z]; Ve2=FLVewfin[q,qp z,l,lp,z];
logs = {s12->4 EC^2, s14->4 EC E4 eta[1,4], s24->4 EC E4 eta[2,4]};
model =
  (* s_vqcd *) Vq0 finvqcd[Qq,Qqp,Qlp,lmu] + Qq^2 Kns[z,Lh1] Vq1 + Qqp^2 Kns[z,Lh1] Vq2 +
  (* s_vewk *) CF (Kns[z,Lh1] Ve1 + Kns[z,Lh1] Ve2) +
  (* s_12 *) CF (
    (Qq^2+Qqp^2) F12 ((P[z1]+Lh12 Pl[z1])(P[z2]+Lh12 Pl[z2]) + (Pb[z1]+Lh12 Pbl[z1])Pl[z2] + Pl[z1](Pb[z2]+Lh12 Pbl[z2]) + Pbl[z1]Pbl[z2])
    + ((-2 Qlp (Qq Log[s14/s12] (Pb[z]+Lh1 Pbl[z]) F1 + Qqp Log[s24/s12] (Pb[z]+Lh1 Pbl[z]) F2)) /. logs)
    + Qq^2 TC[z,Lh1] F1 + Qqp^2 TC[z,Lh1] F2
    + (PA[z]+lmu PlA[z]) (legG[Qq,Qqp,Qlp,lmu,1] F1 + legG[Qq,Qqp,Qlp,lmu,2] F2)
    + finv7[Qq,Qqp,Qlp,lmu] F0 );
(* delta(1-z) F[z..] -> elastic *)
applyd[e_] := Expand[e] /. {delta[1-z] x_ :> (x /. {F1->F0,F2->F0,Vq1->Vq0,Vq2->Vq0,Ve1->Ve0,Ve2->Ve0} /. z->1)};
model = applyd[model] // PowerExpand;

(* ---------- comparison, structure by structure ---------- *)
num = Join[{EC->RandomReal[{20,90},WorkingPrecision->80],E4->RandomReal[{20,90},WorkingPrecision->80],E3->RandomReal[{20,90},WorkingPrecision->80],mu->RandomReal[{30,200},WorkingPrecision->80],
   Qq->RandomReal[{-1,1},WorkingPrecision->80],Qlp->RandomReal[{-1,1},WorkingPrecision->80],CF->4/3},
   Flatten@Table[eta[i,j]->RandomReal[{1/10,9/10},WorkingPrecision->80],{i,1,3},{j,i+1,4}]];
dnreg = {D0[x_] Log[x_]^n_.:>Log[x]^n/(1-x), D1[x_] Log[x_]^n_.:>Log[1-x]Log[x]^n/(1-x), D2[x_] Log[x_]^n_.:>Log[1-x]^2Log[x]^n/(1-x), D3[x_] Log[x_]^n_.:>Log[1-x]^3Log[x]^n/(1-x)};
(* normal form in one variable: c(1) Dk + regular *)
nf[e_,x_] := Module[{ex=Expand[e]//.dnreg,r,out=0},
  r=ex; Do[With[{c=Coefficient[ex,Dk[x]]}, If[c=!=0, With[{c1=Limit[c,x->1]}, out+=c1 Dk[x]; r=r-c Dk[x]+(c-c1)Log[1-x]^kk/(1-x)]]],{Dk,{D0,D1,D2,D3}},{kk,{Position[{D0,D1,D2,D3},Dk][[1,1]]-1}}];
  {out, r}];
cmp[name_,st_,vars_] := Module[{a,b,d},
  a = Coefficient[fin,st]/pref /. If[vars==={}, {E1->EC,E2->EC,Emax->EC,eta[1,2]->1-10^-30}, com] /. Qqp->Qlp-Qq;   (* eta12 -> 1 as a limit: 0*Log[0] terms *)
  b = Coefficient[model,st] /. Qqp->Qlp-Qq;
  If[vars==={}, {a,b} = {a,b} /. {E4->EC,E3->EC}];  (* Born kinematics in the partonic c.m. frame *)
  d = PowerExpand[a-b] /. num;
  d = If[vars==={}, N[d,20],
     If[Length[vars]==1,
       Module[{n=nf[d,vars[[1]]]}, {N[Table[Coefficient[n[[1]],Dk[vars[[1]]]],{Dk,{D0,D1,D2,D3}}],15], N[Table[n[[2]]/.vars[[1]]->zz,{zz,{1/7,1/2,9/10}}],15]}],
       Module[{e=Expand[d]//.dnreg, syms={D0[z1],D1[z1],D0[z2],D1[z2]}},
         N[Table[Max[Abs[Flatten[CoefficientList[e,syms]/.{z1->a1,z2->a2}]]],{a1,{1/5,3/5}},{a2,{2/7,5/7}}],15]]]];
  Print[name,": ",d]];
Print["remaining delta in fin: ",!FreeQ[fin,delta], "   remaining delta in model: ",!FreeQ[model,delta]];
cmp["FLVqcd elastic",Vq0,{}]; cmp["FLVqcd leg1",Vq1,{z}]; cmp["FLVqcd leg2",Vq2,{z}];
cmp["FLVew elastic",Ve0,{}]; cmp["FLVew leg1",Ve1,{z}]; cmp["FLVew leg2",Ve2,{z}];
cmp["FLM double boost",F12,{z1,z2}];
cmp["FLM leg1",F1,{z}]; cmp["FLM leg2",F2,{z}];
cmp["FLM elastic",F0,{}];
Print["VV coefficient (notebook): ",Coefficient[fin,FVVqcdewfin[q,qp,l,lp]]/pref];
