(* usage: see run_verify_W_polecheck.sh. Photon-quark and photon-gluon channels of W_NNLO_polecheck_cleaned.nb against the Fortran
   (mod_xsects_nnlo_s_aq.f90, _ag.f90, W branch; calG_CF_hard; kernels of mod_hoppet_nnlo.f90). Sector weights of the notebook
   (omega tilde) are mapped to the Fortran damping factors as printed. eta12 = 1 - 1e-36 (terms log(1-eta12) log(eta12) vanish). *)
s=Get[Environment["NBOUT"]<>"/wpole_finite.m"];
Get[Environment["FKERN"]];
Get[DirectoryName[$InputFileName]<>"W_polecheck_fortran_defs.wl"];
nbv[nm_,ch_]:=Expand[s[{nm,ch}] /. {\[Alpha]ontwopi->1,\[Alpha]sontwopi->1,E1->EC,E2->EC,Emax->EC} /. Subscript[p_,_]:>p];
L=-2Log[2 EC]; Lh1=L-Log[z]; Lh12=L-Log[z1 z2]; Lel=-2Log[EC];    (* mu = 1; Lel = log(mu^2/EC^2) *)
rnd=Join[{TR->1/2,CF->4/3,Nc->3,EC->71/10,E5->23/10,E6->31/10,eta[1,2]->1-10^-36,eta[1,5]->3/10,eta[2,5]->7/10,eta[1,6]->2/5,eta[2,6]->3/5,
   DQ->37/100,DE->41/100,DA->43/100,W->47/100},Flatten@Table[eta[i,j]->Prime[3 i+j]/50,{i,3,5},{j,i+1,6}]];
symeta={eta[5,1]->eta[1,5],eta[5,2]->eta[2,5],eta[6,1]->eta[1,6],eta[6,2]->eta[2,6]};
dnreg={D0[x_]Log[x_]^n_.:>Log[x]^n/(1-x),D1[x_]Log[x_]^n_.:>Log[1-x]Log[x]^n/(1-x)};
pts=N[{1/9,3/10,1/2,7/10,9/10},40];
posQ[f_]:=TrueQ[N[f/.symeta/.rnd/.{z->3/10,z1->3/10,z2->17/20}]>0];
splitlog[x_]:=x//.{Log[a_Times]/;AllTrue[List@@a,posQ]:>Total[Log/@List@@a],Log[a_^n_]/;posQ[a]:>n Log[a]};
diff[d_]:=Module[{ee=Expand[splitlog[Expand[d]]]//.dnreg},
  Chop[Max[0,Abs[Flatten[Table[CoefficientList[ee,{D0[z],D1[z],D0[z2],D1[z2]}]/.{Qq->qa,Qqp->qb}/.{z->p,z1->p,z2->(1-p/2)}/.symeta/.rnd,
     {p,pts},{qa,{2/3,-1/3,1/7}},{qb,{2/3,-1/3,3/11}}]//N[#,40]&]]],10^-25]];
report[cmp_]:=Do[Print[c[[1]],":  notebook - Fortran = ",diff[c[[2]]-c[[3]]],"     (size of the notebook term ",N[diff[c[[2]]],3],")"],{c,cmp}];

(* ---------------- photon-quark: gamma(1) q(2) -> nu l + g(5) + q'(6); Qqp = charge of the quark from the photon ---------------- *)
e=nbv["aqchri",2]; w=\[Omega]tilde;
Print["== photon-quark (aq)"];
Print["coefficient of FLVQCDfin[gamma,q,nu,l,qp] = ",Coefficient[e,FLVQCDfin[\[Gamma],q,nu,l,qp]]," (= ONLO, sector rvqcd_is: verify_W_polecheck_rv_fortran.wls)"];
usedaq={FLM[qbp z,q,nu,l,z],FLM[qbp z1,q z2,nu,l,z1,z2],FLVQCDfin[qbp z,q,nu,l,z],ONLO[FLM[qbp z,q,nu,l,g,z]],ONLO[FLM[\[Gamma],q z,nu,l,qp,z]],ONLO[FLM[\[Gamma],q,nu,l,qp]],FLVQCDfin[\[Gamma],q,nu,l,qp]};
Print["coverage: notebook remainder minus all compared structures = ",Simplify[e/.Thread[usedaq->0]]];
Q2=Qqp^2;
elnb=Coefficient[e,ONLO[FLM[\[Gamma],q,nu,l,qp]]]/.w[{51,61},{5,6}]->W;
sb=Coefficient[e,FLM[qbp z,q,nu,l,z]]; Print["terms with Log[1-eta12]: ",Select[List@@Expand[sb],!FreeQ[#,Log[1-eta[1,2]]]&]]; Print["aq single boost depends on: ",Select[{eta[1,2],E5,E6,EC,eta[1,5],eta[2,5],D0[z],D1[z]},!FreeQ[sb,#]&]];
report[{
 {"s_12   single boost      FLM[z p1,p2]", sb, CF Nc Q2 tc["TCgq",z,Lh1]},
 {"s_12   double boost      FLM[z1 p1,z2 p2]", Coefficient[e,FLM[qbp z1,q z2,nu,l,z1,z2]], CF Nc Q2 dboost[gqk[z1],qqk[z2],Lh12]},
 {"s_vqcd                   FLVqcd[z p1,p2]", Coefficient[e,FLVQCDfin[qbp z,q,nu,l,z]], Nc Q2 (k["gq",z]+Lh1 k["gq_Lmu",z])},
 {"s_oqcd (resolved gluon)  ONLO FLM[z p1,p2,g]", Coefficient[e,ONLO[FLM[qbp z,q,nu,l,g,z]]]/.w[{51,61},{6,1}]->DQ, Nc Q2 (PqgNLO[z]-Pqg0[z]L+Pqg0[z]Log[eta[1,5]/2]DQ)},
 {"s_oewk boosted           ONLO FLM[p1,z p2,q']", Coefficient[e,ONLO[FLM[\[Gamma],q z,nu,l,qp,z]]]/.w[{52,62},{5,2}]->DE, CF (PA[z]+L PlA[z]-PlA[z]Log[eta[2,6]/2]DE)},
 {"s_oewk elastic, damp2 = omega[{51,61},{5,6}]",   elnb, CF (calGCFhard[EC,E6,eta[1,6],eta[2,6],1-W,W]-3/2 Lel)},
 {"s_oewk elastic, damp2 = 1 - omega[{51,61},{5,6}]", elnb, CF (calGCFhard[EC,E6,eta[1,6],eta[2,6],W,1-W]-3/2 Lel)}}];
Print["sensitivity, aq elastic with 13/2 -> 6: ",N[diff[elnb-CF (calGCFhard[EC,E6,eta[1,6],eta[2,6],W,1-W]-1/2-3/2 Lel)],6]];

(* ---------------- photon-gluon: gamma(1) g(2) -> nu l + q + qbar'; two flavour orderings ---------------- *)
Print["== photon-gluon (ag)"];
e=Expand[s[{"agchri",3}]/.{\[Alpha]ontwopi->1,\[Alpha]sontwopi->1,E1->EC,E2->EC,Emax->EC}];
Do[{q1,q2,Qa,fs,lab,om}=o;
 report[{
  {"s_12   double boost      (leg-1 quark charge "<>ToString[Qa]<>")", Coefficient[e,FLM[z1 Subscript[p1,q1],z2 Subscript[p2,q2],p3,p4,z1,z2]], TR Nc Qa^2 dboost[gqk[z1],gqk[z2],Lh12]},
  {"s_oqcd (resolved quark)  (leg-1 quark charge "<>ToString[Qa]<>")", Coefficient[e,ONLO[FLM[z Subscript[p1,q1],Subscript[p2,g],p3,p4,fs,z]]], Nc Qa^2 (PqgNLO[z]-Pqg0[z]L)},
  {"s_oewk                   (quark from the gluon "<>ToString[q1]<>")", Coefficient[e,ONLO[FLM[Subscript[p1,\[Gamma]],z Subscript[p2,q1],p3,p4,fs,z]]]/.om->1-DA, TR (PqgNLO[z]-Pqg0[z]L+Pqg0[z]Log[eta[2,lab]/2]DA)}}],
 {o,{{q,qp,Qq,Subscript[p6,OverBar[qp]],6,\[Omega]tilde[{52,61},{5,2}]},{qp,q,Qqp,Subscript[p5,OverBar[q]],5,\[Omega]tilde[{51,62},{6,2}]}}}];
Print["all ag structures used: ",Simplify[e/.{FLM[__]->0,ONLO[__]->0}]];
