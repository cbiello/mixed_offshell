(* usage: checks/extract_kernels_fortran.py <mod_hoppet_nnlo.f90> > k.m; NOTES=<notes/fabrizio> FKERN=k.m wolframscript -file derive_gq_general_charges.wls
   gq channel with general charges: the derivation of notes/fabrizio/subtraction_f.wl (RR + RV + PDF renormalisation, charges
   q1..q4 of the all-outgoing process, q1 + q2 + q3 + q4 = 0) is evaluated without specifying the charges.
   Poles must cancel; the coefficient of F(z.p1,p2) is compared with gq s_12 of the Fortran for Z and for W (calG_gq_cc). All lines 0. *)
src = StringSplit[Import[Environment["NOTES"]<>"/subtraction_f.wl","Text"],"\n"];
line[s_] := First@FirstPosition[src,_?(StringContainsQ[#,s]&),{0},{1},Heads->False];
upto[a_,s_] := a - 1 + First@FirstPosition[src[[a;;]],_?(StringContainsQ[#,s]&),{0},{1},Heads->False];
block[a_,b_] := Get[StringToStream[StringRiffle[src[[a;;b]],"\n"]]];
block[line["Quit[]"]+1, line["(*NLO, unsimplified, c.o.m.*)"]-2];                          (* structures and substitutions *)
Do[With[{a=line[s]}, block[a, upto[a,";"]]], {s,{"conv[PAP0[qg],NLOEWK[ns]] =","conv[NLO[gqb],PAP0[qq]] ="}}];
With[{a=line["g[1] qb[2] -> l[3] lb[4] qb[5]"]}, block[a, upto[a,"(*Finite part check*)"]-2]];   (* rr, rv, cv *)
tot = (Collect[rr+rv+cv,{_FLM,ONLO,_FLVfin}] //. subs) /. {ONLO->0,ONNLO->0};
ttt = Expand[Normal[Series[tot,{ep,0,0}]]] //. fixDeltaCH /. fixLog2 /. fixLi2 /. fixLogz;
ttt = ttt /. {four->4,e[1]->EC,e[2]->EC,em->EC,q4->-q1-q2-q3};
c1 = Coefficient[ttt,FLM[q,qb][z,2,3,4]]/(qae qas);
rnd[zz_] := (#[[1]]->N[#[[2]],40])&/@{EC->71/10,e[3]->37/10,e[4]->53/10,mu2->4321,z->zz,eta13->11/50,eta14->13/50,eta23->23/50,eta24->29/50,eta34->41/50,tr->1/2};
pts = {1/9,3/10,1/2,7/10,9/10};
mx[e_,vars_] := Chop[Max[0,Abs[Flatten[Table[Last/@CoefficientRules[Expand[e]/.rnd[p],vars],{p,pts}]]]],10^-25];
Print["poles of the F(z.p1,p2) coefficient, general charges: ",Table[mx[Coefficient[c1,ep,k],{q1,q2,q3}],{k,-4,-1}]];
Print["double boost: charges ",Union[Cases[Coefficient[ttt,FLM[q,qb][z,zb,3,4]],q1|q2|q3,Infinity]],",  virtual term: charges ",
  Union[Cases[Coefficient[ttt,FLVfin[q,qb][EWK,z,2,3,4]],q1|q2|q3,Infinity]]];
fin = Coefficient[c1,ep,0];

(* Fortran: gq s_12, single boost *)
Get[Environment["FKERN"]];
k[n_,z_] := f["mySub_"<>n<>"_func"][z];
L = Log[mu2/(4 EC^2)]; Lh = L - Log[z];
tc = k["TCgq",z] + Lh k["TCgq_Lmu",z] + Lh^2 k["TCgq_Lmu2",z];
Pqg0 = (1-z)^2+z^2;
pt = 2((1-z)z+Pqg0 Log[1-z]) - Pqg0 L;
eta[i_,j_] := ToExpression["eta"<>ToString[i]<>ToString[j]];
calGZ[Qq_,Ql_] := Module[{b,c,E3=e[3],E4=e[4]},      (* calG_QqQl, icoll = 1 *)
  b = 2(PolyLog[2,1-eta13]-PolyLog[2,1-eta14]-PolyLog[2,1-eta23]+PolyLog[2,1-eta24]) + (3+2Log[EC/E3])Log[eta13/eta23] + (3+2Log[EC/E4])Log[eta24/eta14] - 2Log[z]Log[E3/E4 eta13/eta14];
  c = 13-4Zeta[2]+Log[E3/E4]^2+(3+2Log[EC/E3]+2Log[EC/E4])Log[eta34]+2PolyLog[2,1-eta34];
  b Qq Ql + c Ql^2];
calGW[Q1_,Q2_,Ql_,l_] := Module[{El=e[l],lE,le1,le2,li1,li2},      (* calG_gq_cc, ileg = 1, charged lepton l *)
  lE=Log[EC/El]; le1=Log[eta[1,l]]; le2=Log[eta[2,l]]; li1=PolyLog[2,1-eta[1,l]]; li2=PolyLog[2,1-eta[2,l]];
  Ql^2(13/2+3lE+2lE^2-2Pi^2/3+(3/2+2lE)L) + Ql Q1(-3lE-lE^2+3le1+2lE le1+2li1-Pi^2/3-(3+2lE)L) + Ql Q2(-3lE-lE^2+3le2+2lE le2+2li2-Pi^2/3-(3+2lE)L)
  + Q1 Q2(-2Pi^2/3+3L) + (3L-4Zeta[2])Q1^2 - 2 Log[z] Ql Q1 Log[El eta[1,l]/EC]];
Print["Z  (q1,q2,q3,q4) = (-Qq,Qq,Ql,-Ql):          derivation - Fortran = ",
  mx[(fin/.{q1->-Qq,q2->Qq,q3->Ql}) - tr (pt calGZ[Qq,Ql] + Qq^2 tc),{Qq,Ql}]];
(* W: incoming charges of the Fortran Q1 = -q1, Q2 = -q2, charged lepton Ql = Q1 + Q2 *)
Do[Print["W  charged lepton on leg ",l,", any q1, q2:        derivation - Fortran = ",
  mx[(fin/.q3->If[l==4,0,-q1-q2]) - tr (pt calGW[-q1,-q2,-q1-q2,l] + q1^2 tc),{q1,q2}]],{l,{4,3}}];
(* sensitivity: the alternative without the zeta2 term of calG_gq_cc *)
Print["W  same with the -4 zeta2 Q1^2 term removed:     derivation - Fortran = ",
  N[mx[(fin/.q3->0) - tr (pt (calGW[-q1,-q2,-q1-q2,4]+4Zeta[2]q1^2) + q1^2 tc),{q1,q2}],6]];
