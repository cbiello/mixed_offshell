(* usage: checks/extract_kernels_fortran.py <mod_hoppet_nnlo.f90> > k.m; NBDIR=<math/subtraction-notebooks> FKERN=k.m wolframscript -file verify_Z_gq_aq_ag_fortran.wls
   Z: integrated subtraction of the gq, photon-quark and photon-gluon channels (s_12, s_vewk, s_vqcd), Fortran against Z.nb
   (FINITEgq, FINITE\[Gamma]q, FINITE\[Gamma]g; terms without a resolved parton). Every line must be 0 (40-digit numbers, six values of z).
   Kernels are read from the Fortran; the assembly of mod_xsects_nnlo_s_gq/aq/ag.f90 (Vcharge == 0) and calG_QqQl are transcribed. *)
dir=Environment["NBDIR"];
nb=Get[dir<>"/Z.nb"];
cells=Cases[nb,Cell[BoxData[b_],"Input",___]:>b,Infinity];
Do[If[!FreeQ[c,"Ptildegq"|"TCgqBoost"|"Gq"]&&!FreeQ[c,":="]&&FreeQ[c,"FINITEgq"|"FBFortranForm"],ToExpression[c,StandardForm]],{c,cells}];
def[name_]:=ToExpression[First@Cases[nb,RowBox[{name," "...,"=",r__}]:>RowBox[{r}],Infinity],StandardForm];
Fgq=def["FINITEgq"]; Faq=def["FINITE\[Gamma]q"]; Fag=def["FINITE\[Gamma]g"];
Get[Environment["FKERN"]];
k[n_,z_]:=f["mySub_"<>n<>"_func"][z];

(* quark kernels with distributions, mod_hoppet_nnlo.f90 (checked against W.nb in verify_Wp_qq_fortran.wls) *)
P[z_]   := (1-z)-2(1+z)Log[1-z]+(1+z-2/(1-z))Log[z] + 4 D1[z];
Pl[z_]  := (1+z) - 2 D0[z];
Pb[z_]  := ((1-z)-2(1+z)Log[1-z]+(1+z-2/(1-z))Log[z]+4Log[1-z]/(1-z))Log[z];
Pbl[z_] := (1+z-2/(1-z))Log[z];
(* double boost, tables A-D; a = leg 1, b = leg 2, each {K, K_Lmu, K_b, K_b_Lmu}; Lh = log(mu^2/mll^2) *)
dboost[a_,b_,Lh_] := (a[[1]]+Lh a[[2]])(b[[1]]+Lh b[[2]]) + (a[[3]]+Lh a[[4]])b[[2]] + a[[2]](b[[3]]+Lh b[[4]]) + a[[4]]b[[4]];
gqk[z_] := {k["gq",z],k["gq_Lmu",z],k["gq_b",z],k["gq_b_Lmu",z]};
qqk[z_] := {P[z],Pl[z],Pb[z],Pbl[z]};
tc[n_,z_,Lh_] := k[n,z] + Lh k[n<>"_Lmu",z] + Lh^2 k[n<>"_Lmu2",z];
(* calG_QqQl(icoll = 1, chrg = +1) and the gq soft-photon term *)
calG[Qq_,Ql_,z_] := Module[{a=0,b,c},
  b = 2(PolyLog[2,1-eta[1,3]]-PolyLog[2,1-eta[1,4]]-PolyLog[2,1-eta[2,3]]+PolyLog[2,1-eta[2,4]]) + (3+2Log[EC/E3])Log[eta[1,3]/eta[2,3]] + (3+2Log[EC/E4])Log[eta[2,4]/eta[1,4]]
      - 2Log[z]Log[E3/E4 eta[1,3]/eta[1,4]];
  c = 13-4Zeta[2]+Log[E3/E4]^2+(3+2Log[EC/E3]+2Log[EC/E4])Log[eta[3,4]]+2PolyLog[2,1-eta[3,4]];
  a Qq^2 + b Qq Ql + c Ql^2];
Pqg0[z_] := (1-z)^2+z^2;
softgq[Qq_,Ql_,z_] := (2((1-z)z+Pqg0[z]Log[1-z]) - Pqg0[z]Log[mu2/(4EC^2)]) calG[Qq,Ql,z];

(* notebook: FLM[..]/z -> F, alpha = alphas = 1, scale logs as in the code *)
pp[i_]:=Subscript[ToExpression["p"<>ToString[i]],_]|ToExpression["p"<>ToString[i]];
prep[e_] := e /. {\[Alpha]->1,\[Alpha]s->1} /. ONLO[__]->0 /. {
   FLM[z1 . _,z2 . _,__] :> z1 z2 F12, FLM[z . a_,b_,p3,p4] :> z If[FreeQ[a,\[Gamma]],F1,F1a], FLM[a_,z . b_,p3,p4] :> z If[FreeQ[b,\[Gamma]],F2,F2a],
   FLVQEDfin[z . _,__] :> z V1, FLVQCDfin[z . _,__] :> z V1};
single = Log[(2 EC)/\[Mu]] -> -Log[mu2/s]/2 - Log[z]/2;            (* 4 EC^2 = s/z *)
double = EC -> Sqrt[s/z1/z2]/2;
num = (#[[1]]->N[#[[2]],40])&/@Join[{E3->37/10,E4->53/10,EC->71/10,s->9000,mu2->4321,\[Mu]->Sqrt[4321],Qq->2/3,Qe->-1,TR->1/2,CF->4/3,Nc->3},
  Flatten@Table[eta[i,j]->Prime[3 i+j]/50,{i,1,3},{j,i+1,4}]];
dnreg = {D0[x_] Log[x_]^n_.:>Log[x]^n/(1-x), D1[x_] Log[x_]^n_.:>Log[1-x]Log[x]^n/(1-x)};
pts=N[{1/9,3/10,1/2,7/10,9/10,99/100},40];
show[name_,d_,vars_] := Module[{e=Expand[PowerExpand[d]]//.dnreg},
  Print[name,": ", Chop[Max[Abs[N[Flatten[{0,If[Length[vars]==1,
      Table[CoefficientList[e,{D0[z],D1[z]}]/.num/.z->p,{p,pts}],
      Table[CoefficientList[e,{D0[z2],D1[z2]}]/.num/.{z1->p,z2->q},{p,pts},{q,N[{2/7,5/7},40]}]]}],30]]],10^-25]]];
Lh = Log[mu2/s];

(* ---------- gq ---------- *)
g = prep[Fgq];
show["gq  s_vewk   FLVew(z.p1,p2)", (Coefficient[g,V1]/.single) - TR (k["gq",z]+Lh k["gq_Lmu",z]), {z}];
show["gq  s_12     double boost", (Coefficient[g,F12]/.double/.\[Mu]^2->mu2) - TR Qq^2 dboost[gqk[z1],qqk[z2],Lh], {z1,z2}];
gs = Coefficient[g,F1];
show["gq  s_12     single boost, triple-collinear kernel (Qq^2)", ((gs/.Qe->0)/.single) - (TR Qq^2 tc["TCgq",z,Lh]+TR softgq[Qq,0,z]/.EC->Sqrt[s/z]/2), {z}];
show["gq  s_12     single boost, soft-photon term (Q_lep = Qe)", ((gs-(gs/.Qe->0))/.single/.EC->Sqrt[s/z]/2) - (TR (softgq[Qq,Qe,z]-softgq[Qq,0,z])/.EC->Sqrt[s/z]/2), {z}];

(* ---------- photon-quark ---------- *)
aq = prep[Faq];
show["aq  s_vqcd   FLVqcd(z.p1,p2)", (Coefficient[aq,V1]/.single) - Nc Qq^2 (k["gq",z]+Lh k["gq_Lmu",z]), {z}];
show["aq  s_12     double boost", (Coefficient[aq,F12]/.double/.\[Mu]^2->mu2) - CF Nc Qq^2 dboost[gqk[z1],qqk[z2],Lh], {z1,z2}];
show["aq  s_12     single boost, q qb", (Coefficient[aq,F1]/.single) - CF Nc Qq^2 tc["TCgq",z,Lh], {z}];
show["aq  s_12     single boost, gamma gamma", (Coefficient[aq,F2a]/.single) - CF Qq^2 tc["TCaa",z,Lh], {z}];

(* ---------- photon-gluon ---------- *)
h = prep[Fag];
show["ag  s_12     double boost (each of the two flavour orderings)", (Coefficient[h,F12]/2/.double/.\[Mu]^2->mu2) - TR Nc Qq^2 dboost[gqk[z1],gqk[z2],Lh], {z1,z2}];
show["ag  s_12     single boost, gamma gamma", (Coefficient[h,F2a]/.single) - TR Qq^2 tc["TCag",z,Lh], {z}];
Print["terms left unassigned: ", {g,aq,h} /. {F1->0,F2->0,F12->0,F1a->0,F2a->0,V1->0} // Simplify];
