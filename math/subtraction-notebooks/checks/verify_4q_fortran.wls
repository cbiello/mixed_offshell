(* usage: checks/extract_4q_fortran.py <mod_hoppet_nnlo.f90> > k.m; NBDIR=<math/subtraction-notebooks> F4Q=k.m wolframscript -file verify_4q_fortran.wls
   Finite remainders of the four-quark channels (Z.nb: FINITEqqbqqb, FINITEqqqq) against the kernels mySub_ns_qqb, mySub_ns_qq
   (and _Lmu) of f90/src/IntSub/mod_hoppet_nnlo.f90. F4Q holds the four Fortran expressions as f["name"]=...; *)
dir=Environment["NBDIR"];
nb=Get[dir<>"/Z.nb"];
def[name_]:=First@Cases[nb,RowBox[{name,"=",rhs_}]:>rhs,Infinity];
nbk[name_]:=Module[{e=ToExpression[def[name],StandardForm]},
  e=e/(CF Qq^2 \[Alpha] \[Alpha]s)/.FLM[__]->z/2/.Log[(2 EC)/\[Mu]]->-L/2-Log[z]/2;
  {Coefficient[e,L,0],Coefficient[e,L,1]}];
lx=Log[x];lxp=Log[1+x];lxb=Log[1-x];xb=1-x;xp=1+x;zeta2=Zeta[2];zeta3=Zeta[3];
li2x=PolyLog[2,x];li2mx=PolyLog[2,-x];li3x=PolyLog[3,x];li3xb=PolyLog[3,1-x];li3mx=PolyLog[3,-x];
li3xoxp=PolyLog[3,x/(1+x)];li3xbxp=PolyLog[3,(1-x)(1+x)];
Get[Environment["F4Q"]];
pts={1/10,3/10,1/2,7/10,9/10,99/100};
Do[k=nbk[c[[1]]];
  Do[Print[c[[1]]," ",{"finite","log(mu2/s)"}[[i]]," Fortran/notebook at z = ",pts,": ",
     N[Table[(f[c[[1+i]]]/.x->p)/(k[[i]]/.z->p),{p,pts}],12]],{i,2}],
 {c,{{"FINITEqqbqqb","mySub_ns_qqb_func","mySub_ns_qqb_Lmu_func"},{"FINITEqqqq","mySub_ns_qq_func","mySub_ns_qq_Lmu_func"}}}]
