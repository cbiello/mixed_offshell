(* Mathematica transcriptions of the Fortran pieces compared with W_NNLO_polecheck_cleaned.nb, shared by
   verify_W_polecheck_{gq,aq_ag,rv,compiled}_fortran.wls. Load the HOPPET kernels first (FKERN, from extract_kernels_fortran.py).
   verify_W_polecheck_compiled_fortran.wls checks every function below against the compiled code (fortran_polecheck_driver.f90),
   except finEWK/finQCD (private to mod_amplitudes_loop_ppll; their call arguments sit next to OpenLoops calls). *)

k[n_,z_]:=f["mySub_"<>n<>"_func"][z];
tc[n_,z_,Lh_]:=k[n,z]+Lh k[n<>"_Lmu",z]+Lh^2 k[n<>"_Lmu2",z];
gqk[z_]:={k["gq",z],k["gq_Lmu",z],k["gq_b",z],k["gq_b_Lmu",z]};

(* splittings used by the sectors: Pgq_spav, 2((1-z)z + Pqg0 log(1-z)) in s_gq/s_aq/s_ag; Pqq_NLO(0|1), PqqAP_0_R(0|1) of
   mod_splittings_bare with the plus parts written as D0, D1 *)
Pqg0[x_]:=(1-x)^2+x^2; PqgNLO[x_]:=2((1-x)x+Pqg0[x]Log[1-x]);
PA[z_]:=(1-z)-2(1+z)Log[1-z]+4 D1[z]; PlA[z_]:=(1+z)-2 D0[z];
(* the hoppet Pqq kernels (PqqNLO_func, _Lmu, _b, _b_Lmu; checked in fortran_conventions.wl) *)
P[z_]:=(1-z)-2(1+z)Log[1-z]+(1+z-2/(1-z))Log[z]+4 D1[z]; Pl[z_]:=(1+z)-2 D0[z];
Pb[z_]:=((1-z)-2(1+z)Log[1-z]+(1+z-2/(1-z))Log[z]+4Log[1-z]/(1-z))Log[z]; Pbl[z_]:=(1+z-2/(1-z))Log[z];
qqk[z_]:={P[z],Pl[z],Pb[z],Pbl[z]};
(* double boost assembled by get_respdf_hoppet from the four convolutions xPij_A..D (kernel, _Lmu, _b, _b_Lmu) *)
dboost[a_,b_,Lh_]:=(a[[1]]+Lh a[[2]])(b[[1]]+Lh b[[2]])+(a[[3]]+Lh a[[4]])b[[2]]+a[[2]](b[[3]]+Lh b[[4]])+a[[4]]b[[4]];

(* mod_int_sub_nnlo.f90. Energies E<i>, angles eta[i,j]; Q1 = Q_IS(leg 1), Q2 = Q_IS(leg 2), Ql = Qlept on leg il *)
El[il_]:=ToExpression["E"<>ToString[il]];
genG[Q1_,Q2_,Ql_,il_,L_]:=Module[{lE=Log[EC/El[il]],l1=Log[eta[1,il]],l2=Log[eta[2,il]],li1=PolyLog[2,1-eta[1,il]],li2=PolyLog[2,1-eta[2,il]]},
  Ql^2(13/2+3lE+2lE^2-2Pi^2/3+(3/2+2lE)L)+Ql Q1(-3lE-lE^2+3l1+2lE l1+2li1-Pi^2/3-(3+2lE)L)
  +Ql Q2(-3lE-lE^2+3l2+2lE l2+2li2-Pi^2/3-(3+2lE)L)+Q1 Q2(-2Pi^2/3+3L)];                           (* genGnloQCD_cc *)
calGgqcc[Q1_,Q2_,Ql_,il_,L_,z_]:=genG[Q1,Q2,Ql,il,L]+3L Q1^2-4Zeta[2]Q1^2-2Log[z]Ql Q1 Log[El[il] eta[1,il]/EC];   (* calG_gq_cc, ileg = 1 *)
(* calG_ONLOQCD_gq(ig=1,iq=2): resolved quark on leg 5 (energy E5 = Lim_Ei(2) of the hard configuration), Qqp = Qq - Ql *)
calGONLOQCDgq[Qq_,Ql_,il_,L_]:=Module[{Qqp=Qq-Ql,Eq=EC,Emax=EC,LE5oEmax,LE5oEl,LEloEq,LEmaxsqoEqEl,LEmaxsqoElE5,LEmaxsqoEqE5,etal5,etaql,etaq5,etag5,cq,cqp,cqqp},
  LE5oEmax=Log[E5/Emax]; LE5oEl=Log[E5/El[il]]; LEloEq=Log[El[il]/Eq]; LEmaxsqoEqEl=Log[Emax^2/Eq/El[il]]; LEmaxsqoElE5=Log[Emax^2/El[il]/E5]; LEmaxsqoEqE5=Log[Emax^2/Eq/E5];
  etal5=eta[5,il]; etaql=eta[2,il]; etaq5=eta[2,5]; etag5=eta[1,5];
  cq=13/2+LEloEq^2+2PolyLog[2,1-etaql]-Pi^2+(3+2LEmaxsqoEqEl)Log[etaql]-3L/2;
  cqp=13+LE5oEl^2+2PolyLog[2,1-etal5]-2Pi^2/3+3Log[2]-4LE5oEmax Log[2]+(3+2LEmaxsqoElE5)Log[etal5]+(3/2-2LE5oEmax)(LOG1m[etag5]-LOG[etag5]);
  cqqp=-13+2LE5oEl LEloEq-2PolyLog[2,1-etal5]+2PolyLog[2,1-etaq5]-2PolyLog[2,1-etaql]+2Pi^2/3-(3+2LEmaxsqoElE5)Log[etal5]-(3+2LEmaxsqoEqEl)Log[etaql]+(3+2LEmaxsqoEqE5)Log[etaq5];
  Qq^2 cq+Qqp Qq cqqp+Qqp^2 cqp];
calGCFhard[EC_,Ek_,etaik_,etajk_,damp1_,damp2_]:=Module[{lg=Log[EC/Ek]},
  13/2-Pi^2+lg^2+(3+2lg)Log[4etajk]-(3/2+2lg)(damp1 Log[etaik/(1-etaik)]+damp2 Log[etajk/(1-etajk)])+2PolyLog[2,1-etajk]];   (* calG_CF_hard *)

(* mod_limits_rv.f90, without the factor 2/s5i. zF = Lim_KinInv(1) = E5/E1 (= 1 - z of the notebook), eta = eta_15,
   lsi3/lsi4 = log(s_i3/s12), log(s_i4/s12) of the reduced momenta *)
colRVewk[zF_,eta_,lsi3_,lsi4_]:=Module[{omz=1-zF,lz=Log[zF],lomzOz=Log[(1-zF)/zF],leta=Log[eta],ptree},
  ptree=(zF^2+(1-zF)^2)/(1-zF);
  {-1+(5/3 Pi^2+2lomzOz leta-lz^2-lomzOz^2+3lomzOz-2PolyLog[2,zF])ptree, (-2Pi^2-2lsi3 lomzOz)ptree, (-2Pi^2-2lsi4 lomzOz)ptree, ptree}];  (* col_rvewk_is_qg_gen *)
colRVqcd[zF_,eta_]:=Module[{r=colRVewk[zF,eta,0,0]},{r[[1]],r[[4]]}];   (* col_rvqcd_is_qa *)

(* mod_amplitudes_loop_ppll.f90: finite part subtracted from the OpenLoops eps^0 coefficient. Lij = log(s_ij/mu^2), Lij2 = Lij^2 - pi^2
   when i,j are both incoming (<= 2) or both outgoing (> 2) *)
LL2[i_,j_,lij_]:=If[(i<=2&&j<=2)||(i>2&&j>2),lij^2-Pi^2,lij^2];
finEWK[{i1_,i2_,i3_,i4_},{Qq_,Qqbp_,Ql_,Qlb_},Lf_]:=Module[{pr={{i1,i2,-Qq Qqbp},{i1,i3,Ql Qq},{i1,i4,Qlb Qq},{i2,i3,Ql Qqbp},{i2,i4,Qlb Qqbp},{i3,i4,-Ql Qlb}}},
  Total[#[[3]](3 Lf[#[[1]],#[[2]]]-LL2[#[[1]],#[[2]],Lf[#[[1]],#[[2]]]])&/@pr]];      (* get_fin_ewk_qqbllb_gen *)
finQCD[i1_,i2_,Lf_]:=CF(3Lf[i1,i2]-LL2[i1,i2,Lf[i1,i2]]);                             (* get_fin_qcd *)
