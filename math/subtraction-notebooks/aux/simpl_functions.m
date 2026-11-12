(* ::Package:: *)

(*Useful functions for subtraction*)
 
  
ColorReduction = {
          ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],T1,T4]]-> 1/2 (-2 CF ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]]]
         -2 ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],T1,T2]]-2 ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],T1,T3]]),
          ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],T2,T3]]->  1/2 ((CA-CF-CF-CA) ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]]]
         -2 ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],T1,T2]] -2 ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],T1,T3]]),
          ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],T2,T4]]-> 1/2 ((CF+CA-CF-CA)ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]]]
         +2 ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],T1,T3]]),
          ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],T3,T4]]-> 1/2 ((CF+CF-CA-CA)ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]]]
          +2 ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],T1,T2]])
          };
          
 SectorReductionElastic = {
          ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]] wtilde4151[5,1]]->ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]]]
          -  ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]] wtilde4251[5,1]]- ONLO[2 \[CapitalDelta]34 FLM[ Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]] wtilde4351[5,1]],
          ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]] wtilde4252[5,2]]-> ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]]]
          -ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]] wtilde4152[5,2]]-ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]] wtilde4352[5,2]],
           ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]] wtilde4151[4,5]]-> ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]]]
          -ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]] wtilde4252[4,5]]-ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]] wtilde4353[4,5]],
          ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]] wtilde4353[5,3]]->ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]]]
          -  ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]] wtilde4153[5,3]]- ONLO[2 \[CapitalDelta]34 FLM[ Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g]] wtilde4253[5,3]]
          };
          
 SectorReductionBoost = {
          ONLO[2 \[CapitalDelta]34 FLM[z Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],z] wtilde4151[5,1]]->ONLO[2 \[CapitalDelta]34 FLM[z Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],z]]
          - ONLO[2 \[CapitalDelta]34 FLM[z Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],z] wtilde4251[5,1]]- ONLO[2 \[CapitalDelta]34 FLM[z Subscript[p1, q],Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],z] wtilde4351[5,1]],
          ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],z Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],z] wtilde4252[5,2]]-> ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],z Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],z]]
          -ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],z Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],z] wtilde4152[5,2]]-ONLO[2 \[CapitalDelta]34 FLM[Subscript[p1, q],z Subscript[p2, \!\(\*OverscriptBox[\(q\), \(_\)]\)],Subscript[p3, g],Subscript[p4, g],z] wtilde4352[5,2]]
          };

DnProperties[expr_] := Module[
    {res},
    res = expr//Expand;
    res = res/.{D0[z]Log[z]->Log[z]/(1-z)};
    res = res/.{D1[z] Log[z]->Log[1-z]/(1-z) Log[z]};
    res = res/.{D1[z] Log[z]^2 -> Log[1-z]/(1-z) Log[z]^2};
    res = res/.{D2[z] Log[z] -> Log[1-z]^2/(1-z) Log[z]};
    res = res/.{D0[z] Log[z]^2-> Log[z]^2/(1-z)};
    res = res/.{D0[z] Log[z]^3-> Log[z]^3/(1-z)};
 res   
];
 
ApplyDelta[expr_] := Module[
    {res},
    (* res = Collect[expr,{_delta},GG]; *)
    res = Collect[expr, {delta[a__] delta[b__], _delta}, GG] /. GG[x_] delta[y_] :>  GG[Factor[x]] * delta[y];
    res = res //. delta[1 - z1_] delta[1 - z2_] GG[b_] :> GG[b /. {z1 -> 1, z2 -> 1}];
    res = res //. delta[1 - a_] GG[b___] :>  GG[b /. a -> 1];
    res = res //. {FLM[a_, b_, c_, 1, 1] -> FLM[a, b, c], FLM[a_, b_, c_, 1] -> FLM[a, b, c]};
    res = res //. {ONLO[h__ FLM[a__, 1, 1]] -> ONLO[h FLM[a]], ONLO[h__ FLM[a__, 1]] -> ONLO[h FLM[a]]};
    (* res = res //. {FLM[a, b, c, d, 1] -> FLM[a, b, c, d]};*)
    res = res //. {FLVfin[a__,1]-> FLVfin[a]};
    res = res //. {GG[FLM[a__, 1, zn_] x___] :> (FLM[a, z] x /. zn -> z)};
    res = res //. {GG[FLM[a__, zn_, 1] x___] :> (FLM[a, z] x /. zn -> z)};
    res = res //. {GG[FLM[a__, 1, zn_]] :> (FLM[a, z]), GG[FLM[a__, zn_, 1]] :> (FLM[a, z] )};
   (* I added the line below *)
    res = res //. {GG[FLM[a__, z1] x___] :>  (FLM[a, z] x /. z1 -> z), GG[FLM[a__, z2] x___] :> (FLM[a, z] x /. z2 -> z)};
    res = res //. {GG[FLM[a__, z1]] :> (FLM[a, z]/. z1 -> z), GG[FLM[a__, z2]] :> (FLM[a, z] /. z2 -> z)};
   (*******)
    res = res /. GG -> Identity;
    res
    ];
    
Applydelta[expr_]:=Module[{res},
res=Collect[expr,{delta[a__] delta[b__],_delta},GG]/.GG[x_] delta[y_]:>Distribute[GG[x]]*delta[y];
res=res//Expand;
res=res//.delta[1-z1_] delta[1-z2_] GG[FLM[b___,z1_,z2_]]:>GG[FLM[b]/.z1->1/.z2->1];
res=res//.delta[1-z1_] delta[1-z2_] GG[x___ FLM[b___,z1_,z2_]]:>GG[(x FLM[b])/.{z1->1,z2->1}];
res=res//.delta[1-z_] GG[FLM[b___,z_]]:>GG[FLM[b]/.z->1];
res=res//.delta[1-z_] GG[x___ FLM[b___,z_]]:>GG[(x FLM[b])/.z->1];
res=res//.delta[1-z_] GG[FLM[b___,z_,c_]]:>GG[FLM[b,c]/.z->1];
res=res//.delta[1-z_] GG[x___ FLM[b___,z_,c_]]:>GG[(x FLM[b,c])/.z->1];
res=res//.GG[x_ FLM[a_,b_,c_, d_,z1]]:>GG[(x FLM[a,b,c,d,z])/.{z1->z}];
res=res//.GG[FLM[a_,b_,c_,d_,z1]]:>GG[FLM[a,b,c,d,z]/.{z1->z}];
res=res//.GG[x_ FLM[a_,b_,c_,d_,z2]]:>GG[(x FLM[a,b,c,d,z])/.{z2->z}];
res=res//.GG[FLM[a_,b_,c_,d_,z2]]:>GG[FLM[a,b,c,d,z]/.{z2->z}];
(* ONLO *)
res=res//.delta[1-z_] GG[x___ ONLO[FLM[b___,z_]]]:>GG[x ONLO[FLM[b]/.z->1]];
       res=res//.delta[1-z_] GG[ONLO[x___ FLM[b___,z_]]]:>GG[ONLO[x FLM[b]/.z->1]];
       res=res//.delta[1-z_] GG[y___ ONLO[x___ FLM[b___,z_]]]:>GG[y ONLO[x FLM[b]/.z->1]];
      (**)
       res=res//.delta[1-z_] GG[ONLO[y___ FLM[b___,z_,c_]]]:>GG[ONLO[y FLM[b,c]/.z->1]];
       res=res//.delta[1-z_] GG[x___ ONLO[FLM[b___,z_,c_]]]:>GG[x ONLO[FLM[b,c]/.z->1]];
       res=res//.delta[1-z_] GG[x___ ONLO[y___ FLM[b___,z_,c_]]]:>GG[x ONLO[y FLM[b,c]/.z->1]];
      (**)
       res=res//.GG[x___ ONLO[FLM[a_,b_,c_,z1]]]:>GG[x ONLO[FLM[a,b,c,z]]/.z1->z];  
       res=res//.GG[ONLO[x___ FLM[a_,b_,c_,z1]]]:>GG[ONLO[x FLM[a,b,c,z]]/.z1->z];
       res=res//.GG[y___ ONLO[x___ FLM[a_,b_,c_,z1]]]:>GG[y ONLO[x FLM[a,b,c,z]]/.z1->z];
     (**)
     res=res//.GG[x___ ONLO[FLM[a_,b_,c_,z2]]]:>GG[x ONLO[FLM[a,b,c,z]]/.z2->z];
     res=res//.GG[ONLO[x___ FLM[a_,b_,c_,z2]]]:>GG[ONLO[x FLM[a,b,c,z]]/.z2->z];
      res=res//.GG[y___ ONLO[x___ FLM[a_,b_,c_,z2]]]:>GG[y ONLO[x FLM[a,b,c,z]]/.z2->z];  
(* FLVfin *)
res=res//.delta[1-z_] GG[x___ FLVQEDfin[b___,z_]]:>GG[x FLVQEDfin[b]/.z->1];  
res=res//.delta[1-z_] GG[x___ FLVewfin[b___,z_]]:>GG[x FLVewfin[b]/.z->1];  
res=res//.delta[1-z_] GG[x___ FLVQCDfin[b___,z_]]:>GG[x FLVQCDfin[b]/.z->1];  
res=res//.delta[1-z_] GG[x___ FLVqcdfin[b___,z_]]:>GG[x FLVqcdfin[b]/.z->1];  
       
res=res/.GG->Identity;
res];


CleanMess[expr_] := Module[
    {res, FactFF},
    FactFF = {FF[a_Integer] :> a, FF[a_Rational] :> a,
    FF[a_Integer x_] :> a FF[x], FF[a_Rational x_] :> a FF[x]};
    res = expr;
    res = res /. PolyGamma[i_, j_] :> FunctionExpand[PolyGamma[i, j]];
    check = Coefficient[res,EulerGamma]//Simplify ;
    If[check === 0, res = res /. EulerGamma->0, Print["EulerGamma doesnt cancel PD!"]];
    Print["applying delta"];
    res = res // Applydelta;
    res = res /. Dispatch[LogSubs];
    Print["collecting 1..."];
    res = Collect[
      res, {asontwopi, bra\[Alpha]s, tag[a___], FLM[a___],
       FLVfin[a___], ONLO[a___],Onlo4[a___], _Log,
       PolyLog[a_, b_], \[Pi], _Zeta, _D0, _D1, \[Beta]0}, FF];
    res = res /. {Log[a_Integer] :> PowerExpand[Log[a]]};
    Print["factoring..."];
    res = res /. FF[x_] :> FF[Factor[x]] /. FactFF /. FF -> GG;
    Print["collecting 2..."];
    res = Collect[
        res, {asontwopi, bra\[Alpha]s, FLM[a___], FLVfin[a___],
         ONLO[a___],Onlo4[a___], _GG, _D0, _D1, \[Beta]0 (*,_Log,PolyLog[a_,
         b_],\[Pi],_Zeta*)}, FF] /. FF[x_] :> FF[Factor[x]] /. FactFF;
    res = res /. Dispatch[SpecialLogSubs] /. {FF -> Identity, GG -> Identity};
    res
    ];
    
ArrangeLogs[expr_] := Module[
    {res},
    res = expr;
    res = res /. Log[a_/ eta[c_, d_]] :> Log[a] - Log[eta[c, d]];
    res = res /. Log[a_ eta[c_, d_]] :> Log[a] + Log[eta[c, d]];
    res = res /. Log[Emax/a_] :> 1/2 Log[Emax^2/a^2];
    res = res /. Log[a_/mu2] :> -Log[Emax^2/a] + Log[Emax^2/mu2];
    res = res /. Log[a_ s_Rational] :> Log[a]+PowerExpand[Log[s]];
    res = Collect[
      res, {asontwopi, bra\[Alpha]s, FLM[a___], _Log,
       PolyLog[a_, b_], \[Pi], _Zeta, _D0, _D1}, Factor];
    res = res /. x___ Log[Emax^2/a_^2] + x___ Log[Emax^2/b_^2] :> 2 x Log[Emax^2/(a b)];
    res = Collect[
      res, {asontwopi, bra\[Alpha]s, FLM[a___], _Log,
       PolyLog[a_, b_], \[Pi], _Zeta, _D0, _D1}, Factor];
    res
    ];
