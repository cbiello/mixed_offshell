(* usage: see run_verify_W_polecheck.sh. The coefficient-1 terms ONLO[FLVQEDfin[p1,p2,p3,p4,pq]] (gq) and ONLO[FLVQCDfin[gamma,q,nu,l,qp]]
   (photon-quark) of W_NNLO_polecheck_cleaned.nb against the real-virtual sectors rvewk_is (xsect_nnlo_rvewk_is_gq) and rvqcd_is
   (xsect_nnlo_rvqcd_is_aq). ONLO[F] = Int (F - C F), so two things must agree:
   (1) the hard term: the finite remainder F that the Fortran builds from OpenLoops, res1(0) - fin * res0, against the notebook definition
       V = (pole term) FLM + F (FVQEDgq, FLRVQCDγq). OpenLoops normalisation = e^(eps gamma)/Gamma(1-eps) up to O(eps^3), so
       fin must be the eps^0 coefficient of the pole term with alpha_B -> 1; same convention as the four-point Born.
   (2) the counterterm: the collinear limit of F implied by the notebook. RVgq/RVγq subtract the collinear limit of the full RV,
       integrated:  -(aSB/eps) TR G(1-eps)^2/G(1-2eps) PqgNLO[z] FVQEDgqz1/(2E1)^(2eps)
                    +(aSB aB/(2eps)) TR Qqp^2 G(1-2eps)G(1-eps)/G(1-3eps) PqgRV[z] FLM[z p1]/(2E1)^(4eps)
       With s15 = (2E1)^2 (1-z) eta15, Int deta eta^(-1-eps)(1-eta)^(-eps) = -G(1-eps)^2/(eps G(1-2eps)),
       Int deta eta^(-1-2eps)(1-eta)^(-eps) = -G(1-2eps)G(1-eps)/(2eps G(1-3eps)) and PqgNLO = (1-z)^(-2eps) Pqg(z,eps), the unintegrated
       limit per unit of the tree splitting is  Pqg(z,eps) FVQEDgqz1 - aB Qqp^2 s15^(-eps) (1-z)^(3eps) PqgRV[z] FLM[z p1].
       Subtracting the limit of the pole term, (pole term)|_{5||1} Pqg(z,eps) FLM[z p1], leaves C FLVfin; its poles must cancel and its
       finite part must equal the Fortran col_rvewk_is_qg_gen / col_rvqcd_is_qa (the notebook frame has E1 = E2 = EC, s12 = 4 EC^2). *)
s=Get[Environment["NBOUT"]<>"/wpole_finite.m"];
Get[Environment["FKERN"]];
Get[DirectoryName[$InputFileName]<>"W_polecheck_fortran_defs.wl"];
ep=\[Epsilon];
DownValues[nbPqgRV]=s[{"PqgRV",1}]/.PqgRV->nbPqgRV; DownValues[nbPqgNLO]=s[{"PqgNLO",1}]/.PqgNLO->nbPqgNLO;
aB=Exp[ep EulerGamma]/Gamma[1-ep];
num[x_]:=N[x,30];
chk[name_,d_,vars_]:=Print[name,":  max |notebook - Fortran| = ",ToString[NumberForm[N[Max[Abs[Flatten[{num[d/.vars]}]]]],3]]];

(* ================= gq: FVQEDgq = P5 FLM[pq] + FLVQEDfin[pq] ================= *)
FV5=s[{"FVQEDgq",1}]; FV4=s[{"FVQEDgqz1dividedbyz",1}];
Print["FVQEDgq = P5 FLM + c FLVQEDfin with c = ",Coefficient[FV5,FLVQEDfin[p1,p2,p3,p4,pq]],";  FVQEDgqz1dividedbyz: c = ",Coefficient[FV4,FLVQEDfin[z p1,p2,p3,p4,z]]];
P5=Coefficient[FV5,FLM[p1,p2,p3,p4,pq]]/.\[Alpha]B->aB; P4z=Coefficient[FV4,FLM[z p1,p2,p3,p4,z]]/.\[Alpha]B->aB;

(* (1) hard finite remainder: Fortran res_ewkloop_g_gq_gen subtracts get_fin_ewk_qqbllb_gen(5,2,3,4,...) with Lij of legs 2..5 *)
sym=<|{2,5}->s52,{3,5}->s53,{4,5}->s54,{2,3}->s23,{2,4}->s24,{3,4}->s34,{1,2}->s12,{1,3}->s13,{1,4}->s14|>;
Lf5[i_,j_]:=Log[sym[Sort[{i,j}]]]; Lf4=Lf5;
delta5=SeriesCoefficient[P5/aB/.{s52->Exp[l52],s54->Exp[l54],s24->Exp[l42]},{ep,0,0}]/.{l52->Log[s52],l54->Log[s54],l42->Log[s24]};
(* charges in the Fortran call, read off res_ewkloop_g_gq_gen: W+ (lepton leg 4) "21 2 -> 12 -11 1": (-Qdn,Qup,0,-Q_lep);
   "21 -1 -> 12 -11 -2": (Qup,-Qdn,0,-Q_lep);  W- (lepton leg 3) "21 -2 -> 11 -12 -1": (Qdn,-Qup,Q_lep,0); "21 1 -> 11 -12 2": (-Qup,Qdn,Q_lep,0).
   Q_lep = -1. In notebook charges (Qq = leg 2, Qqp = outgoing quark 5, Ql = Qq - Qqp): (-Qqp, Qq, Ql on the lepton leg). *)
Qdn=-1/3; Qup=2/3; Qlep=-1;
fcalls={{"W+ g u",{-Qdn,Qup,0,-Qlep},{Qup,Qdn},4},{"W+ g dbar",{Qup,-Qdn,0,-Qlep},{-Qdn,-Qup},4},
        {"W- g ubar",{Qdn,-Qup,Qlep,0},{-Qup,-Qdn},3},{"W- g d",{-Qup,Qdn,Qlep,0},{Qdn,Qup},3}};
to3={s54->s53,s24->s23};
Do[{nm,args,{qq,qqp},il}=c;
  dnb=If[il==4,delta5,delta5/.{s54->s53,s24->s23}]/.{Qq->qq,Qqp->qqp,Ql->qq-qqp};
  dF=finEWK[{5,2,3,4},args,Lf5];
  chk["(1) hard, finite remainder of the 5-point EW loop, "<>nm<>" (eps^0 of the pole term - fin)",dnb-dF,{s52->7/3,s54->11/5,s53->13/7,s24->5/3,s23->17/9,s34->19/7}],
 {c,fcalls}];
(* four-point Born, same convention (res_ewkloop_qqb_gen, W+ "-1 2 -> 12 -11": (1,2,3,4,-Qdn,Qup,0,-Q_lep)), sanity *)
delta4=SeriesCoefficient[P4z/aB/.{z->1}/.{s12->Exp[l12],s14->Exp[l14],s24->Exp[l42]},{ep,0,0}]/.{l12->Log[s12],l14->Log[s14],l42->Log[s24]};
chk["(1) four-point Born W+, same relation (sanity)",(delta4/.{Qq->Qup,Qqp->Qdn,Ql->Qup-Qdn})-finEWK[{1,2,3,4},{-Qdn,Qup,0,-Qlep},Lf4],{s12->7/3,s14->11/5,s24->5/3,s13->13/7,s23->17/9,s34->19/7}];

(* (2) collinear counterterm *)
Pqgep=Normal[Series[nbPqgNLO[z] (1-z)^(2ep),{ep,0,3}]]/.ErrorPqgNLO->0;
P5lim=P5/.{s52->(1-z)s12,s54->(1-z)s14};
s15=(2E1)^2(1-z)eta;
ctEW=Pqgep (P4z-P5lim) - aB Qqp^2 s15^(-ep) (1-z)^(3ep) nbPqgRV[z];
serEW=Series[ctEW/.Ql->Qq-Qqp,{ep,0,0}];
vars=Flatten[Table[{Qq->qa,Qqp->qb,z->zz,eta->ee,E1->33/10,s14->29/10,s13->37/10,s24->41/10,s23->43/10},{qa,{2/3,-1/3}},{qb,{2/3,-1/3,1/7}},{zz,{1/9,3/10,7/10}},{ee,{1/13,1/2}}],3];
dd[x_]:=Max[Abs[num[x/.#]&/@vars]];
Print["(2) gq: poles of the implied C51 FLVQEDfin: 1/eps^2 -> ",Simplify[SeriesCoefficient[serEW,-2]],",  1/eps (max over the test points) -> ",N[dd[SeriesCoefficient[serEW,-1]/.s12->4E1^2]]];
finEW=SeriesCoefficient[serEW,0]/.s12->4E1^2;
(* Fortran: res0 (limgen(1) Qi^2 + limgen(2) Q3 Qi + limgen(3) Q4 Qi) + res1 limgen(4); Qi = Q_IS(leg 1) = -Qqp; ratio to the res1 term *)
{r1,r2,r3,r4}=colRVewk[1-z,eta,Log[s13/(4E1^2)],Log[s14/(4E1^2)]];
Print["(2) gq W+ (lepton leg 4): C51 FLVQEDfin, notebook - Fortran (coefficient of FLM[z p1] over that of FLVQEDfin[z p1]) = ",ToString[NumberForm[N[dd[finEW/Pqg0[z]-(r1 Qqp^2+r3 (Qq-Qqp)(-Qqp))/r4]],3]]];
finEW3=SeriesCoefficient[Series[(ctEW/.{s54->s53,s14->s13,s24->s23})/.Ql->Qq-Qqp,{ep,0,0}],0]/.s12->4E1^2;
Print["(2) gq W- (lepton leg 3): C51 FLVQEDfin, notebook - Fortran = ",ToString[NumberForm[N[dd[finEW3/Pqg0[z]-(r1 Qqp^2+r2 (Qq-Qqp)(-Qqp))/r4]],3]]];
Print["(2) sensitivity, 5/3 pi^2 -> 3/2 pi^2 in limgen(1): ",N[dd[finEW/Pqg0[z]-((r1-Pi^2/6 (z^2+(1-z)^2)/z) Qqp^2+r3 (Qq-Qqp)(-Qqp))/r4],4]];

(* ================= photon-quark: FLRVQCDγq = P6 FLM[gamma,q,nu,l,qp] + FLVQCDfin ================= *)
FQ6=s[{"FLRVQCD\[Gamma]q",2}]; FQ4=s[{"FLRVQCD\[Gamma]qz1",2}];
Print["FLRVQCDγq: c = ",Coefficient[FQ6,FLVQCDfin[\[Gamma],q,nu,l,qp]],";  FLRVQCDγqz1: c = ",Coefficient[FQ4,FLVQCDfin[z qbp,q,nu,l,z]]];
P6=Coefficient[FQ6,FLM[\[Gamma],q,nu,l,qp]]/.\[Alpha]sB->aB; P4q=Coefficient[FQ4,FLM[z qbp,q,nu,l,z]]/.\[Alpha]sB->aB;
(* (1) ol_res_qcdloop_a_aq_gen: get_fin_qcd(2,5) with Lij of legs 2,5 (Fortran leg 5 = notebook leg 6) *)
delta6=SeriesCoefficient[P6/aB/.s62->Exp[l62],{ep,0,0}]/.l62->Log[s62];
chk["(1) hard, finite remainder of the 5-point QCD loop, photon-quark",delta6-finQCD[2,5,Log[s62]&]/.CF->4/3,{s62->7/3}];
(* (2) C61: Fortran res0 limcol(1) CF + res1 limcol(2) *)
s16=(2E1)^2(1-z)eta;
ctQ=Pqgep (P4q-(P6/.s62->(1-z)s12)) - aB CF s16^(-ep) (1-z)^(3ep) nbPqgRV[z];
serQ=Series[ctQ,{ep,0,0}];
Print["(2) photon-quark: poles of the implied C61 FLVQCDfin: 1/eps^2 -> ",Simplify[SeriesCoefficient[serQ,-2]],",  1/eps (max over the test points) -> ",N[dd[SeriesCoefficient[serQ,-1]/.{s12->4E1^2,CF->4/3}]]];
finQ=SeriesCoefficient[serQ,0]/.s12->4E1^2;
{q1,q2}=colRVqcd[1-z,eta];
Print["(2) photon-quark: C61 FLVQCDfin, notebook - Fortran = ",ToString[NumberForm[N[dd[(finQ/Pqg0[z]-CF q1/q2)/.CF->4/3]],3]]];
