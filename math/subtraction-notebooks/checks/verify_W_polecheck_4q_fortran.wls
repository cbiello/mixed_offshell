(* usage: see run_verify_W_polecheck.sh. Four-quark chapters of W_NNLO_polecheck_cleaned.nb: coefficient of FLM[p1, z p2] divided by
   CF Q^2, written in Lh = log(mu^2/mll^2), against mySub_ns_qqb, mySub_ns_qq (and _Lmu) of mod_hoppet_nnlo.f90. Ratio 1 = same kernel. *)
s=Get[Environment["NBOUT"]<>"/wpole_finite.m"];
Print["keys: ",Select[Keys[s],StringMatchQ[#[[1]],"four*"]&]];
lx=Log[x];lxp=Log[1+x];lxb=Log[1-x];xb=1-x;xp=1+x;zeta2=Zeta[2];zeta3=Zeta[3];
li2x=PolyLog[2,x];li2mx=PolyLog[2,-x];li3x=PolyLog[3,x];li3xb=PolyLog[3,1-x];li3mx=PolyLog[3,-x];
li3xoxp=PolyLog[3,x/(1+x)];li3xbxp=PolyLog[3,(1-x)(1+x)];
Get[Environment["F4Q"]];
fk[n_,zz_]:=f["mySub_ns_"<>n<>"_func"]/.x->zz; fkl[n_,zz_]:=f["mySub_ns_"<>n<>"_Lmu_func"]/.x->zz;
pts={1/10,3/10,1/2,7/10,9/10,99/100};
Do[{nm,ch}=key;
  e=s[key]/.{\[Alpha]ontwopi->1,\[Alpha]sontwopi->1}/.Subscript[p_,_]:>p;
  st=Union[Cases[e,FLM[___],Infinity]]; c=Coefficient[e,st[[1]]];
  qch=First[Select[{Qqpp,Qqp},!FreeQ[c,#]&]];
  (* notebook mu, E2: write in Lh = log(mu^2/mll^2), mll^2 = 4 E2^2 z (single boost of leg 2) *)
  cc=Expand[c/(CF qch^2)/.\[Mu]->Sqrt[4 E2^2 z Exp[Lh]]]; cc=Expand[PowerExpand[cc]]//Simplify;
  Print["== ",nm,", chapter ",ch,": structure ",InputForm[st],", charge ",qch,"^2, remaining variables ",Select[{E2,z,Lh},!FreeQ[cc,#]&]];
  e=e/.{\[Alpha]B->1,\[Alpha]sB->1};
  Do[dev=Max[Abs[N[Table[(Coefficient[cc/.{\[Alpha]B->1,\[Alpha]sB->1},Lh,i-1]/.z->p)/({fk,fkl}[[i]][n,p])-1,{p,pts}],30]]];
     Print["-> ",{"Lh^0","Lh^1"}[[i]],": max |notebook/Fortran mySub_ns_",n," - 1| over z = ",pts,": ",ScientificForm[N[dev,3]]],{n,{"qqb","qq"}},{i,2}];
  Print["-> Lh^2 coefficient: ",Simplify[Coefficient[cc,Lh,2]]],
 {key,Select[Keys[s],StringMatchQ[#[[1]],"four*"]&]}];
