(* usage: see run_verify_W_polecheck.sh. The transcriptions in W_polecheck_fortran_defs.wl (used by the other verify_W_polecheck
   scripts to compare the notebook with the Fortran) against the compiled code: fortran_polecheck_driver.f90 evaluates the HOPPET
   kernels (private copies), Pqq_NLO, PqqAP_0_R, Pgq_spav, calG_gq_cc, calG_ONLOQCD_gq, calG_CF, col_rvewk_is_qg_gen, col_rvqcd_is_qa
   at the points below (FVALS = its output). Relative differences, double precision. *)
Get[Environment["FKERN"]];
lx=Log[x];lxp=Log[1+x];lxb=Log[1-x];xb=1-x;xp=1+x;zeta2=Zeta[2];zeta3=Zeta[3];
li2x=PolyLog[2,x];li2mx=PolyLog[2,-x];li3x=PolyLog[3,x];li3xb=PolyLog[3,1-x];li3mx=PolyLog[3,-x];
li3xoxp=PolyLog[3,x/(1+x)];li3xbxp=PolyLog[3,(1-x)(1+x)];
Get[Environment["F4Q"]];   (* f["mySub_ns_qq(b)(_Lmu)_func"] as expressions in x *)
k4[n_,xx_]:=f["mySub_"<>n<>"_func"]/.x->xx;
Get[DirectoryName[$InputFileName]<>"W_polecheck_fortran_defs.wl"];
Get[Environment["FVALS"]];
fv=N/@fv;
rel[a_,b_]:=N[Abs[a-b]/Max[1,Abs[b]]];
out[nm_,l_]:=Print[StringPadRight[nm,62],"max rel. difference ",ToString[NumberForm[Max[l],2]]];
xs=fv["xs"]; lmu=fv["lmu"];
plusv[e_,xx_]:={e/.{D1[x]->Log[1-x]/(1-x),D0[x]->1/(1-x)}, Coefficient[e,D1[x]]Log[1-x]/(1-x)+Coefficient[e,D0[x]]/(1-x)}/.x->xx;
out["HOPPET kernels gq, gq_Lmu, gq_b, gq_b_Lmu, TCgq(_Lmu,_Lmu2)",Flatten@Table[rel[k[n,xs[[i]]],fv[{n,i}]],{i,5},{n,{"gq","gq_Lmu","gq_b","gq_b_Lmu","TCgq","TCgq_Lmu","TCgq_Lmu2"}}]];
out["HOPPET kernels ns_qq, ns_qqb (and _Lmu)",Flatten@Table[rel[k4[n,xs[[i]]],fv[{n,i}]],{i,5},{n,{"ns_qq","ns_qq_Lmu","ns_qqb","ns_qqb_Lmu"}}]];
out["Pqq_NLO (0),(1) = PA regular + plus, plus part",Flatten@Table[MapThread[rel,{plusv[PA[x],xs[[i]]],{fv[{"PA0",i}],fv[{"PA1",i}]}}],{i,5}]];
out["PqqAP_0_R (0),(1) = PlA",Flatten@Table[MapThread[rel,{plusv[PlA[x],xs[[i]]],{fv[{"PlA0",i}],fv[{"PlA1",i}]}}],{i,5}]];
out["Pgq_spav = Pqg0",Table[rel[Pqg0[xs[[i]]],fv[{"Pgq_spav",i}]],{i,5}]];

{EC0,E30,E40,E50,e13,e14,e35,e45}=fv["kin"]; {x2,zz}=fv["x2z"];
kin={EC->EC0,E3->E30,E4->E40,E5->E50,eta[1,3]->e13,eta[2,3]->1-e13,eta[1,4]->e14,eta[2,4]->1-e14,eta[5,3]->e35,eta[5,4]->e45,eta[1,5]->x2,eta[2,5]->1-x2};
lg={LOG1m[a_]:>Log[1-a],LOG[a_]:>Log[a]};
out["calG_gq_cc (W+ dbar u, u dbar; W- d ubar, ubar d; two lmu)",Flatten@Table[{q1,q2,il,ql}=fv[{"chg",c}]; il=Round[il];
   Table[rel[calGgqcc[q1,q2,ql,il,lmu[[j]],zz]/.kin,fv[{"calG_gq_cc",c}][[j]]],{j,2}],{c,4}]];
out["calG_ONLOQCD_gq (same four, quark on leg 2)",Flatten@Table[{q1,q2,il,ql}=fv[{"chg",c}]; il=Round[il];
   Table[rel[calGONLOQCDgq[q2,ql,il,lmu[[j]]]/.lg/.kin,fv[{"calG_ONLOQCD_gq",c}][[j]]],{j,2}],{c,4}]];
out["calG_CF (hard)",{rel[calGCFhard[EC0,E50,0.31,0.69,0.41,0.59],fv["calG_CF"]]}];
{zF,et,zb}=fv["colkin"];
out["col_rvewk_is_qg_gen (x s5i/2)",MapThread[rel,{colRVewk[zF,et,Log[E30 e13/EC0],Log[E40 e14/EC0]],fv["col_rvewk"]}]];
out["col_rvqcd_is_qa (x s5i/2)",MapThread[rel,{colRVqcd[zF,et],fv["col_rvqcd"]}]];
(* sensitivity: a wrong transcription shows up *)
Print["sensitivity, calG_gq_cc without the -4 zeta2 term: ",ToString@NumberForm[{q1,q2,il,ql}=fv[{"chg",1}]; il=Round[il];rel[(calGgqcc[q1,q2,ql,il,lmu[[1]],zz]+4Zeta[2]q1^2)/.kin,fv[{"calG_gq_cc",1}][[1]]],2]];
