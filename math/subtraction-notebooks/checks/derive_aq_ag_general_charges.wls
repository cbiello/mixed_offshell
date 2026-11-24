(* usage: NOTES=<notes/fabrizio> wolframscript -file derive_aq_ag_general_charges.wls
   Photon-quark and photon-gluon channels with general charges: the derivations of notes/fabrizio/subtraction_f.wl
   (charges q1..q4 of the all-outgoing process) are evaluated without specifying the charges. Poles must cancel.
   In the photon-gluon channel q1 and q2 both stand for the charge of the quark (gamma gamma term), so q2 -> q1 there.
   The terms with a resolved parton (ONLO) are compared with the Fortran s_oqcd, s_oewk of aq and ag. All differences 0. *)
src = StringSplit[Import[Environment["NOTES"]<>"/subtraction_f.wl","Text"],"\n"];
line[s_] := First@FirstPosition[src,_?(StringContainsQ[#,s]&),{0},{1},Heads->False];
upto[a_,s_] := a - 1 + First@FirstPosition[src[[a;;]],_?(StringContainsQ[#,s]&),{0},{1},Heads->False];
block[a_,b_] := Get[StringToStream[StringRiffle[src[[a;;b]],"\n"]]];
block[line["Quit[]"]+1, line["(*NLO, unsimplified, c.o.m.*)"]-2];                          (* structures and substitutions *)
Do[With[{a=line[s]}, block[a, upto[a,";"]]], {s,{"conv[NLO[aqb],PAP0[qq]] =","conv[PAP0[qg],NLOQCD[ns]] =",
  "conv[NLO[aqb],PAP0[qg]] =","conv[NLO[aq],PAP0[qg]] =","conv[PAP0[qg],NLO[qbg]] =","conv[PAP0[qg],NLO[qg]] ="}}];
section[s_] := With[{a=line[s]}, block[a, upto[a,"(*Finite part check*)"]-2]];
section["\\[Gamma][1] qb[2] -> l[3] lb[4] g[5] qb[6]"]; totaq = rr+rv+cv //. subs;
section["\\[Gamma][1] g[2] -> l[3] lb[4] q[5] qb[6]"];  totag = rr+cv //. subs;
same = {four->4,e[1]->EC,e[2]->EC,em->EC,w[1,1]->1-w22,w[2,2]->w22,w56[1,1]->1-d2,w56[2,2]->d2};
fin[t_] := Expand[Normal[Series[t,{ep,0,0}]]] //. fixDeltaCH /. fixLog2 /. fixLi2 /. fixLogz /. same;
str[e_] := Cases[Variables[e],_FLM[__]|_FLVfin[__]];
rnd = (#[[1]]->N[#[[2]],40])&/@{EC->71/10,e[5]->23/10,mu2->4321,z->3/10,zb->7/10,eta15->3/10,eta25->7/10,w22->9/20,d2->13/20,
  tr->1/2,Cf->4/3,xn->3,DD[0,zb]->17/10,DD[1,zb]->19/10};
mx[e_] := Chop[Max[0,Abs[Last/@CoefficientRules[Expand[e]/.rnd,{q1,q2}]]],10^-25];
chg[e_] := Union[Cases[e,q1|q2|q3|q4,Infinity]];
report[name_,t_,rl_:{}] := Module[{on,rest},
  on = fin[Coefficient[t/.ONNLO->0,ONLO]]; rest = fin[t/.{ONLO->0,ONNLO->0}] /. rl;
  Print[name,": poles of the resolved-parton terms ",Table[mx[Coefficient[on,ep,k]],{k,-4,-1}],", of the other terms ",Table[mx[Coefficient[rest,ep,k]],{k,-4,-1}]];
  Do[Print["   charges in ",s,": ",chg[Coefficient[on,s]]],{s,str[on]}];
  Do[Print["   charges in ",s,": ",chg[Coefficient[rest,s]]],{s,str[rest]}];
  Coefficient[on,ep,0]];
onaq = report["photon-quark",totaq]; onag = report["photon-gluon",totag,q2->q1];

(* Fortran, mod_xsects_nnlo_s_aq.f90 and _ag.f90; Q1 = charge of the quark from the photon splitting *)
L = Log[mu2/(4 EC^2)]; Pqg0[x_] := (1-x)^2+x^2; PqgNLO[x_] := 2((1-x)x+Pqg0[x] Log[1-x]);
PqqNLO[x_] := (1-x) - 2(1+x)Log[1-x] + 4 DD[1,x]; Pqq0R[x_] := 1+x - 2 DD[0,x];
calGCF = 13/2 - Pi^2 + Log[EC/e[5]]^2 + (3+2Log[EC/e[5]])Log[4 eta25] -
  (3/2+2Log[EC/e[5]])((1-d2)Log[eta15/(1-eta15)] + d2 Log[eta25/(1-eta25)]) + 2PolyLog[2,1-eta25];
c[e_,s_] := Coefficient[e,s] /. {qae->1,qas->1};
Print["aq s_oqcd (resolved gluon):         derivation - Fortran = ",
  mx[c[onaq,FLM[q,qb,l,lb,g][z,2,3,4,5]] - xn q1^2 (PqgNLO[z] + Pqg0[z] Log[eta15/2](1-w22) - Pqg0[z] L)]];
Print["aq s_oewk (resolved quark), boosted: derivation - Fortran = ",
  mx[c[onaq,FLM[a,qb,l,lb,qb][1,zb,3,4,5]] - Cf (PqqNLO[zb] - Pqq0R[zb] Log[eta25/2] w22 + Pqq0R[zb] L)]];
Print["aq s_oewk (resolved quark), elastic: derivation - Fortran = ",
  mx[c[onaq,FLM[a,qb,l,lb,qb][1,2,3,4,5]] - Cf (calGCF - 3/2 Log[mu2/EC^2])]];
Print["ag s_oqcd (resolved quark, q g):    derivation - Fortran = ",
  mx[c[onag,FLM[q,g,l,lb,q][z,2,3,4,5]] - xn q1^2 (PqgNLO[z] - Pqg0[z] L)]];
Print["ag s_oewk (resolved quark, a q):    derivation - Fortran = ",
  mx[c[onag,FLM[a,q,l,lb,q][1,zb,3,4,5]] - tr (PqgNLO[zb] + Pqg0[zb] Log[eta25/2] w22 - Pqg0[zb] L)]];
Print["sensitivity, aq elastic with 13/2 -> 6: ",N[mx[c[onaq,FLM[a,qb,l,lb,qb][1,2,3,4,5]] - Cf (calGCF - 1/2 - 3/2 Log[mu2/EC^2])],6]];
