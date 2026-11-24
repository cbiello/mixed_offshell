(* usage: NBDIR=<math/subtraction-notebooks> NBOUT=<dir with wnb.mx from run_W_nb.wls> wolframscript -file bench_fortran_conventions_W.wls
   q qb' finite remainder of fortran_conventions.wl (Fortran charges and legs, charged lepton on leg ilept) against W.nb
   mapped with fromWnb[ilept], for W+ (ilept = 4) and W- (ilept = 3), general charges. Every number must be 0 (elastic: ~1e-30).
   The last block uses wrong conventions on purpose and must be nonzero. *)
Get[Environment["NBDIR"]<>"/fortran_conventions.wl"];
Get[Environment["NBOUT"]<>"/wnb.mx"];
pref = \[Alpha]ontwopi \[Alpha]sontwopi;
num[il_] := Join[{EC->RandomReal[{20,90},WorkingPrecision->80], E3->RandomReal[{20,90},WorkingPrecision->80], E4->RandomReal[{20,90},WorkingPrecision->80], mu->RandomReal[{30,200},WorkingPrecision->80],
   Qq->RandomReal[{-1,1},WorkingPrecision->80], Qlept->RandomReal[{-1,1},WorkingPrecision->80], CF->4/3},
   Flatten@Table[eta[i,j]->RandomReal[{1/10,9/10},WorkingPrecision->80],{i,1,3},{j,i+1,4}]];
dnreg = {D0[x_] Log[x_]^n_.:>Log[x]^n/(1-x), D1[x_] Log[x_]^n_.:>Log[1-x]Log[x]^n/(1-x), D2[x_] Log[x_]^n_.:>Log[1-x]^2Log[x]^n/(1-x), D3[x_] Log[x_]^n_.:>Log[1-x]^3Log[x]^n/(1-x)};
(* one variable: endpoint coefficients of D0..D3 and the regular rest at three points *)
nf[e_,x_] := Module[{ex=Expand[e]//.dnreg,r,out=0},
  r=ex; Do[With[{c=Coefficient[ex,Dk[x]]}, If[c=!=0, With[{c1=Limit[c,x->1]}, out+=c1 Dk[x]; r=r-c Dk[x]+(c-c1)Log[1-x]^kk/(1-x)]]],{Dk,{D0,D1,D2,D3}},{kk,{Position[{D0,D1,D2,D3},Dk][[1,1]]-1}}];
  {out, r}];
diff[fin_,model_,st_,vars_,il_] := Module[{a,b,d,nm=num[il]},
  a = Coefficient[fin,st]/pref /. {E1->EC,E2->EC,Emax->EC,eta[1,2]->If[vars==={},1-10^-30,1]} /. Qqp->Qlept-Qq;
  b = Coefficient[model,st] /. Qqp->Qlept-Qq;
  If[vars==={}, {a,b} = {a,b} /. El[il]->EC];           (* Born kinematics in the partonic c.m. frame *)
  d = PowerExpand[a-b] /. nm;
  Chop[Max[Abs[Flatten[{Which[vars==={}, N[d,20],
     Length[vars]==1, With[{n=nf[d,z]}, {N[Table[Coefficient[n[[1]],Dk[z]],{Dk,{D0,D1,D2,D3}}],15], N[Table[n[[2]]/.z->zz,{zz,{1/7,1/2,9/10}}],15]}],
     True, With[{e=Expand[d]//.dnreg}, N[Table[Max[Abs[Flatten[CoefficientList[e,{D0[z1],D1[z1],D0[z2],D1[z2]}]/.{z1->a1,z2->a2}]]],{a1,{1/5,3/5}},{a2,{2/7,5/7}}],15]]]}] /. (D0|D1|D2|D3)[_]->1]], 10^-25]];
sts = {{"FLVqcd elastic",FLVqcdfin[p1,p2,p3,p4],{}}, {"FLVqcd leg1",FLVqcdfin[p1 z,p2,p3,p4,z],{z}}, {"FLVqcd leg2",FLVqcdfin[p1,p2 z,p3,p4,z],{z}},
  {"FLVew elastic",FLVewfin[p1,p2,p3,p4],{}}, {"FLVew leg1",FLVewfin[p1 z,p2,p3,p4,z],{z}}, {"FLVew leg2",FLVewfin[p1,p2 z,p3,p4,z],{z}},
  {"FLM double boost",FLM[p1 z1,p2 z2,p3,p4,z1,z2],{z1,z2}}, {"FLM leg1",FLM[p1 z,p2,p3,p4,z],{z}}, {"FLM leg2",FLM[p1,p2 z,p3,p4,z],{z}}, {"FLM elastic",FLM[p1,p2,p3,p4],{}}};
bench[name_,map_,il_,ilmodel_] := Module[{fin=finitesimpl/.map, model=PowerExpand[applyd[Total[finW[ilmodel]]]]},
  Print[name]; Do[Print["   ",s[[1]],": ",diff[fin,model,s[[2]],s[[3]],il]],{s,sts}]];
SeedRandom[1];
bench["W+ (ilept = 4)", fromWnb[4], 4, 4];
bench["W- (ilept = 3)", fromWnb[3], 3, 3];
Print["sensitivity, must be nonzero:"];
bench["  W.nb charges taken as incoming (no sign flip)", Join[{Qqbp->Qqp,Qe->Qlept}, Rest[Rest[Rest[fromWnb[4]]]]], 4, 4];
bench["  W- legs not swapped in the Fortran remainder", fromWnb[3], 3, 4];
