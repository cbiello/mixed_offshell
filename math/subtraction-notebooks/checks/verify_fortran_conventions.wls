(* usage: NBDIR=<math/subtraction-notebooks> FVALS=<fortran_values.m> wolframscript -file verify_fortran_conventions.wls
   (run_verify_fortran_conventions.sh builds the driver and runs this)
   fortran_conventions.wl against the compiled Fortran: hoppet kernels (pieces REAL, VIRT, DELTA), Pqq_NLO, PqqAP_0_R and
   the elastic functions fin_ns_vqcd_cc, genGnloQCD_cc, legfinalcalG_cc, fin_ns_v7_cc for W+, W- and non-conserved charges.
   Every line prints the largest relative difference; all must be ~1e-15. The last block uses a wrong kernel or leg on purpose. *)
Get[Environment["NBDIR"]<>"/fortran_conventions.wl"];
Get[Environment["FVALS"]];
xs = {7/100, 31/100, 62/100, 93/100}; lmus = {37/100, -121/100};
toL = {D0[x_]:>1/(1-x), D1[x_]:>Log[1-x]/(1-x), D2[x_]:>Log[1-x]^2/(1-x), D3[x_]:>Log[1-x]^3/(1-x), delta[_]->0};
(* hoppet pieces of a distribution K = R + sum c_k D_k + d delta: REAL = x (R + sum c_k L_k), VIRT = -x sum c_k L_k, DELTA = d *)
pieces[K_] := {z (K /. toL), -z ((K - (K /. {D0[_]->0,D1[_]->0,D2[_]->0,D3[_]->0})) /. toL), Coefficient[Expand[K], delta[1-z]]};
rel[a_,b_] := Abs[a-b]/Max[1,Abs[b]];
sci[x_] := ToString@ScientificForm[N[x],2,NumberFormat->(If[#3=="",#1,#1<>"e"<>#3]&)];
kern[Kn_] := <|"ns"->Kn[z,0], "nsL"->Coefficient[Kn[z,Lh],Lh], "P"->P[z], "Pl"->Pl[z], "Pb"->Pb[z], "Pbl"->Pbl[z],
   "TC0"->TC0[z], "TC1"->TC1[z], "TC2"->TC2[z]|>;
Print["hoppet kernels, max rel. diff per kernel:"];
ks = kern[Kns]; Do[Print["   ",nm,": ",sci@Max@Table[rel[fv[{nm,i,ip}], N[pieces[ks[nm]][[ip]] /. z->xs[[i]], 30]], {i,4}, {ip,3}]], {nm, Keys[ks]}];
(* Pqq_NLO, PqqAP_0_R: (-1) delta, (0) full function, (1) plus part *)
spl[K_] := {Coefficient[Expand[K],delta[1-z]], K /. toL, (K - (K /. {D0[_]->0,D1[_]->0})) /. toL};
Do[Print["   ",s[[1]]," (Pqq_NLO / PqqAP_0_R): ",sci@Max@Table[rel[fv[{s[[1]],i,k}], N[spl[s[[2]]][[k+2]] /. z->xs[[i]], 30]], {i,4}, {k,-1,1}]], {s,{{"PA",PA[z]},{"PlA",PlA[z]}}}];
kin = {EC->453/10, E3->381/10, E4->517/10, eta[1,3]->27/100, eta[2,3]->73/100, eta[1,4]->64/100, eta[2,4]->36/100, eta[3,4]->81/100};
labels = {"W+ u db","W+ db u","W+ c sb","W- d ub","W- ub d","u ub, Qlept = +1 (not conserved)"};
el[c_, L_] := Module[{q1,q2,il,ql}, {q1,q2,il,ql} = fv[{"chg",c}]; il = Round[il];
  N[{finNsVqcdCC[il,L], genGnloQCDcc[il,L], legfinalcalGcc[il,L,1], legfinalcalGcc[il,L,2], finNsV7cc[il,L]} /. kin /. {Qq->q1, Qqp->q2, Qlept->ql}, 30]];
fvel[c_, k_] := fv[{#,c}][[k]] & /@ {"vqcd","genG","legG1","legG2","v7"};
Print["elastic functions (vqcd, genG, legG1, legG2, v7), max rel. diff over both lmu:"];
Do[Print["   ",labels[[c]],": ",sci/@Max/@Transpose@Table[MapThread[rel, {fvel[c,k], el[c,lmus[[k]]]}], {k,2}]], {c,6}];
Print["sensitivity, must be far above 1e-15:"];
pert = kern[Kns]; pert["TC1"] = TC1[z] + 1/1000 D1[z];
Print["   TC1 with an extra 1e-3 D1: ",sci@Max@Table[rel[fv[{"TC1",i,ip}], N[pieces[pert["TC1"]][[ip]] /. z->xs[[i]], 30]], {i,4}, {ip,3}]];
Print["   genGnloQCD_cc with the charged lepton on the other leg: ",
  sci@Max@Table[rel[fv[{"genG",c}][[1]], N[genGnloQCDcc[7-Round[fv[{"chg",c}][[3]]],lmus[[1]]] /. kin /. {Qq->fv[{"chg",c}][[1]],Qqp->fv[{"chg",c}][[2]],Qlept->fv[{"chg",c}][[4]]},30]], {c,5}]];
