(* usage: see run_verify_W_polecheck.sh. gq channel of W_NNLO_polecheck_cleaned.nb (g(1) q(2) -> nu(3) l(4) + q'(5) + photon(6),
   Qq = quark on leg 2, Qqp = flavour of the quark from the gluon, Ql = Qq - Qqp on leg 4) against the Fortran of mod_xsects_nnlo_s_gq.f90
   (W branch: calG_gq_cc, calG_ONLOQCD_gq, kernels of mod_hoppet_nnlo.f90), structure by structure; general charges, five z, 40 digits.
   W-: the notebook with the lepton moved to leg 3 against the Fortran with ilept = 3. *)
s=Get[Environment["NBOUT"]<>"/wpole_finite.m"];
Get[Environment["FKERN"]];
Get[DirectoryName[$InputFileName]<>"W_polecheck_fortran_defs.wl"];
e=Expand[s[{"gqchri",1}] /. {\[Alpha]ontwopi->1,\[Alpha]sontwopi->1, E1->EC,E2->EC,Emax->EC}];
(* notebook charges: Qq = quark on leg 2, Qqp = flavour of the quark from the gluon (qbp enters on leg 1), lepton on leg 4 *)
L = -2 Log[2 EC]; Lh1 = L - Log[z];  Lh12 = L - Log[z1 z2];      (* mu = 1 in the notebook *)
(* notebook -> Fortran: Q1 = -Qqp (the qbar' from the gluon), Q2 = Qq, Ql = Qq - Qqp; photon label 6 (notebook) = 5 (Fortran) *)
F1 = TR(Qqp^2 tc["TCgq",z,Lh1] + (PqgNLO[z]-Pqg0[z]L) calGgqcc[-Qqp,Qq,Qq-Qqp,4,L,z]);
F12 = TR Qq^2 dboost[gqk[z1],qqk[z2],Lh12];
V1  = TR (k["gq",z]+Lh1 k["gq_Lmu",z]);
OQe = calGONLOQCDgq[Qq,Qq-Qqp,4,L];
OQb = Qq^2 (PA[z]+L PlA[z]);
OEW = TR (PqgNLO[z]-Pqg0[z]L+Pqg0[z]Log[eta[1,6]/2] DAMP);
nbONLOel = Coefficient[e,ONLO[FLM[p1,p2,p3,p4,pq]]] + Coefficient[e,ONLO[FLM[p1,p2,p3,p4,pq]Log[1-eta[5,1]]]] LOG1m[eta[1,5]] + Coefficient[e,ONLO[FLM[p1,p2,p3,p4,pq]Log[eta[5,1]]]] LOG[eta[1,5]];
cmp = {
 {"s_12   single boost      FLM[z p1,p2]",   Coefficient[e,FLM[p1 z,p2,p3,p4,z]], F1},
 {"s_12   double boost      FLM[z1 p1,z2 p2]",Coefficient[e,FLM[p1 z1,p2 z2,p3,p4,z1,z2]], F12},
 {"s_vewk                   FLVew[z p1,p2]", Coefficient[e,FLVQEDfin[p1 z,p2,p3,p4,z]], V1},
 {"s_oqcd elastic           ONLO FLM[pq]",   nbONLOel, OQe},
 {"s_oqcd boosted           ONLO FLM[p1,z p2,pq]", Coefficient[e,ONLO[FLM[p1,p2 z,p3,p4,pq,z]]], OQb},
 {"s_oewk (is + fs_53/54)   ONLO FLM[z p1,p2,pa]", Coefficient[e,ONLO[FLM[p1 z,p2,p3,p4,pa,z]]]/.\[Omega]tilde[{\[Gamma],q},{q,1}]->DAMP, OEW}};
Print["ONLO[FLVQEDfin[pq]] coefficient = ",Coefficient[e,ONLO[FLVQEDfin[p1,p2,p3,p4,pq]]]," (sector rvewk_is: verify_W_polecheck_rv_fortran.wls)"];
rnd=Join[{TR->1/2,EC->71/10,E4->53/10,E5->23/10,DAMP->37/100,eta[1,2]->1},Flatten@Table[eta[i,j]->Prime[3 i+j]/50,{i,1,5},{j,i+1,6}]];
symeta={eta[5,1]->eta[1,5],eta[5,2]->eta[2,5],eta[5,4]->eta[4,5],eta[6,1]->eta[1,6],eta[6,2]->eta[2,6]};
dnreg={D0[x_]Log[x_]^n_.:>Log[x]^n/(1-x),D1[x_]Log[x_]^n_.:>Log[1-x]Log[x]^n/(1-x)};
pts=N[{1/9,3/10,1/2,7/10,9/10},40];
ev[x_,zz_]:=N[x/.{LOG1m[a_]:>Log[1-a],LOG[a_]:>Log[a]}/.symeta/.rnd/.zz,40];
posQ[f_]:=TrueQ[N[f/.symeta/.rnd/.{z->3/10,z1->3/10,z2->17/20}]>0];
splitlog[x_]:=x//.{Log[a_Times]/;AllTrue[List@@a,posQ]:>Total[Log/@List@@a],Log[a_^n_]/;posQ[a]:>n Log[a]};
diff[d_]:=Module[{ee=Expand[splitlog[Expand[d]]]//.dnreg,cl},
  Chop[Max[0,Abs[Flatten[Table[ CoefficientList[#,{D0[z],D1[z],D0[z2],D1[z2]}]&/@{ee} /.{Qq->qa,Qqp->qb}/.{z->p,z1->p,z2->(1-p/2)}/.symeta/.rnd,
     {p,pts},{qa,{2/3,-1/3,1/7}},{qb,{2/3,-1/3,3/11}}]/.{LOG1m[a_]:>Log[1-a],LOG[a_]:>Log[a]}//N[#,40]&]]],10^-25]];
Do[Print[c[[1]],":  notebook - Fortran = ",diff[c[[2]]-c[[3]]],"     (size of the notebook term ",N[Max[Abs[Flatten[{diff[c[[2]]]}]]],3],")"],{c,cmp}];
(* sensitivity *)
Print["sensitivity, -4 zeta2 term dropped:  ",N[diff[cmp[[1,2]]-(F1+TR(PqgNLO[z]-Pqg0[z]L)4Zeta[2]Qqp^2)],6]];
(* W-: the same notebook with the lepton moved to leg 3 against the Fortran with ilept = 3 *)
to3={E4->E3,eta[i_,4]:>eta[i,3],eta[4,i_]:>eta[3,i]};
rnd=Join[rnd,{E3->53/10}];
Print["W- (lepton on leg 3), single boost:  ",diff[(cmp[[1,2]]/.to3)-(TR(Qqp^2 tc["TCgq",z,Lh1]+(PqgNLO[z]-Pqg0[z]L)calGgqcc[-Qqp,Qq,Qq-Qqp,3,L,z]))]];
Print["W- (lepton on leg 3), s_oqcd elastic: ",diff[(nbONLOel/.to3)-calGONLOQCDgq[Qq,Qq-Qqp,3,L]]];
(* coverage: every structure of the notebook remainder is compared above or is the RV term *)
used={FLM[p1 z,p2,p3,p4,z],FLM[p1 z1,p2 z2,p3,p4,z1,z2],FLVQEDfin[p1 z,p2,p3,p4,z],ONLO[FLM[p1,p2,p3,p4,pq]],ONLO[FLM[p1,p2,p3,p4,pq]Log[1-eta[5,1]]],
  ONLO[FLM[p1,p2,p3,p4,pq]Log[eta[5,1]]],ONLO[FLM[p1,p2 z,p3,p4,pq,z]],ONLO[FLM[p1 z,p2,p3,p4,pa,z]],ONLO[FLVQEDfin[p1,p2,p3,p4,pq]]};
Print["coverage: notebook remainder minus all compared structures = ",Simplify[e/.Thread[used->0]]];
