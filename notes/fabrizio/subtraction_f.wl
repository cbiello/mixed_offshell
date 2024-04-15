(* ::Package:: *)

Quit[]


(* taken almost verbatim from w_subtraction.wl, with some strategic changes *)
(* FLM, FLV assumed to contain the (d-dim) averaging factor *)


(* ::Title:: *)
(*Formulas for mixed QCD-EW NC*)


(* qi: charge of the all-outgoing process, with antiparticles having opposite charge *)
(* e.g.: u ub -> l lb --> 0 -> ub u l lb so
   q1 = -Qup, q2 = Qup, q3 = Qlep, q4 = -Qlep *)
   
(* checkmark: agrees with Chiara *)


(* charges always refer to the number in the FLM that multiplies them *)


(* general color rules: PqqI -> {Cf,q^2}; PqgI -> {tr, xn q^2}; PgqI -> q^2 [\[Gamma] -> q transition, QED only] *)


(* ::Chapter::Closed:: *)
(*Structures and substitutions*)


(* ::Subsection:: *)
(*Generic structures and substitutions*)


asontwopi = qas * Gamma[1-ep]/Exp[ep EulerGamma];
aemontwopi = qae * Gamma[1-ep]/Exp[ep EulerGamma];


(* ::Subsubsection::Closed:: *)
(*Single-collinear, tree-level*)


(* factor out either qas or qae and charge, either Cf or qi^2 *)
(* R/D: real/delta; Di/Df: act on the phase space *)

prefC = Gamma[1-ep]^2/Gamma[1-2ep];

getC = {
(* -- initial state -- *)
CiR[qq_,i_,z_] :>  1/ep * (four e[i]^2/mu2)^(-ep) * prefC * (-calP[qq,R,2][z]),
CiD[ij_,i_] :>     1/ep * (four e[i]^2/mu2)^(-ep) * prefC * 
	(If[ij===qq || ij===gg,2 ddt[1,em^2/e[i]^2],0]),
(* -- *)
Ci[qq_,i_,z_] :>  CiR[qq,i,z] + CiD[qq,i]delta[1-z],
CiDi[qq_,i_,z_]:> CiR[qq,i,z]/prefC + CiD[qq,i]/prefC delta[1-z],
(* -- *)
(* -- final state -- *)
Cf[qq_,i_] :> 1/ep * (four e[i]^2/mu2)^(-ep)* prefC * (ga[qq,R,2,2]) + CiD[qq,i],
CfDf[qq_,i_]:> Cf[qq,i]/prefC
};


(* ::Subsubsection::Closed:: *)
(*Single-collinear, one loop bit*)


(* factor out qas*qae and charge, Cf * qi^2 *)
(* angular averaging factor added below, so here we present CiDi for RV *)

getCRV = {
CiRVDi[qq_,i_,z_]:> 1/ep (four e[i]^2/mu2)^(-2ep) calPRV[qq,R,z]};

calPRV[qq_,R,z_] := -Gamma[1-ep]^3 Gamma[1+ep]/Gamma[1-2ep] Pt1L[qq,z];


(* ::Subsubsection::Closed:: *)
(*RR collinear compositions*)


(* factor out qas * qae, charges
   and overall prefactor coming from the angular ordering *)

getTimes = {
	times[qq,qq,i_,z_] :> 
	1/ep^2 *(4 e[i]^2/mu2)^(-2ep) * (
	 conv[cP[qq,R,2],zm2ep cP[qq,R,2]][z]-2 ddt[z^2,em^2/e[i]^2] calP[qq,R,2][z]+
	 -2 ddt[1,em^2/e[i]^2] calP[qq,R,2][z]+(2ddt[1,em^2/e[i]^2])^2 delta[1-z]),
	(* -- *)
	times[qq,qg,i_,z_] :>
	1/ep^2 * (4 e[i]^2/mu2)^(-2ep) * (
	 conv[cP[qq,R,2],zm2ep cP[qg,R,2]][z]- 2 ddt[z^2,em^2/e[i]^2] calP[qg,R,2][z]),
	(* -- *)
	times[gq,qq,i_,z_] :> 
	1/ep^2 * (4 e[i]^2/mu2)^(-2ep) * (
	 conv[cP[gq,R,2],zm2ep cP[qq,R,2]][z] +
	 -2 ddt[1,em^2/e[i]^2] calP[gq,R,2][z]),
	(* -- *)
	times[gq,qg,i_,z_] :> 
	1/ep^2 * (4 e[i]^2/mu2)^(-2ep) conv[cP[gq,R,2],zm2ep cP[qg,R,2]][z] };
	 
(* w.r.t. old notebook, factor out gamma factors Nc and use ga[qq,R,2,2] 
   instead of gatqq = Nc ga[qq,R,2,2] *)
getTimesIF = {
	timesIF[qgI,qqF,i_,z_] :>
	1/ep^2 * (4 e[i]^2/mu2)^(-2ep) * (
	 ga[qq,R,2,2]*(-calP[qg,R,4][z]) - 2 * 
	 (calP[qg,R,4][z]-(em^2/e[i]^2)^(-ep) calP[qg,R,2][z])/2/ep),
    (* -- *)
	timesIF[gqI,qqF,i_,z_] :>
	1/ep^2 * (4 e[i]^2/mu2)^(-2ep) * (
	 ga[qq,R,2,2]*(-calP[gq,R,4][z]) - 2 * 
	 (calP[gq,R,4][z]-(em^2/e[i]^2)^(-ep) calP[gq,R,2][z])/2/ep)
};


(* ::Subsubsection::Closed:: *)
(*Convolutions between AP and single-collinear*)


(* factor out all the charges *)
getConv = {
conv[Ci[qq],PAP0[qq]][i_,z_] :> 1/ep Gamma[1-ep]^2/Gamma[1-2ep] * 
	(four e[i]^2/mu2)^(-ep)* (-conv[cP[qq,R,2],zm2ep cAP0[qq]][z] + 
	PAP0R[qq,z]*2*ddt[z^2,em^2/e[i]^2]+PAP0D[qq]*2ddt[1,em^2/e[i]^2]delta[1-z]),
(* -- *)
conv[Ci[qq],PAP0[qg]][i_,z_] :> 1/ep Gamma[1-ep]^2/Gamma[1-2ep] * 
	(four e[i]^2/mu2)^(-ep)* (-conv[cP[qq,R,2],zm2ep cAP0[qg]][z] + 
	PAP0[qg,z]*2*ddt[z^2,em^2/e[i]^2]),
(* -- *)
conv[Ci[gq],PAP0[qq]][i_,z_] :> 1/ep Gamma[1-ep]^2/Gamma[1-2ep] * 
	(four e[i]^2/mu2)^(-ep)* (-conv[cP[gq,R,2],zm2ep cAP0[qq]][z]),
(* -- *)
conv[Ci[gq],PAP0[qg]][i_,z_] :> 1/ep Gamma[1-ep]^2/Gamma[1-2ep] * 
	(four e[i]^2/mu2)^(-ep)* (-conv[cP[gq,R,2],zm2ep cAP0[qg]][z])
};


(* ::Subsubsection::Closed:: *)
(*ddt*)


getddt = {ddt[a_, b_] :> (a^(-ep) - b^(-ep))/2/ep};


(* ::Subsubsection::Closed:: *)
(*Expressions for splitting functions and anomalous dimensions*)


(* calP[a,R,i] ->  Pi[a,z](1-z)^(-i ep), without delta function *)
(* ga[a,R,i,j]   -> -\int_0^1 z^(-i ep) (1-z)^(-j ep) P[a,z], without delta function *)

(* calP[qg,R,i] -> z * Pqq(1-1/z)/(1-ep) * (1-z)^(-i ep), where all Pij are color-stripped *)
(* calP[gq,R,i] -> z * <Pgq(1/z)> * (1-ep) * (1-z)^(-i ep), where all Pij are color-stripped *)
(* i.e. calP[gq,R,i] = [(1+(1-z)^2)/z -ep z] (1-z)^(-i ep) *)

getP = {calP[qq,R,2][z_]:> -1 - z + 2*DD[0, z] + ep*(-1 + z - 4*DD[1, z] + 
   (2 + 2*z)*Log[1 - z]) + 
 ep^2*(4*DD[2, z] + (2 - 2*z)*Log[1 - z] + 
   (-2 - 2*z)*Log[1 - z]^2) + 
 ep^3*((-8*DD[3, z])/3 + (-2 + 2*z)*Log[1 - z]^2 + 
   (4/3 + (4*z)/3)*Log[1 - z]^3) + 
 ep^4*((4*DD[4, z])/3 + (4/3 - (4*z)/3)*Log[1 - z]^3 + 
   (-2/3 - (2*z)/3)*Log[1 - z]^4),
(* -- *)
calP[qg,R,2][z_] :> 1-2*z+2*z^2 + ep*(-1+(1-2*z+ 2*z^2)*(1 - 2*Log[1 - z])) + 
  ep^2*(2*(-1 + z)*z - 4*(-1 + z)*z*Log[1 - z] + (2 - 4*z + 4*z^2)*Log[1 - z]^2) + 
  ep^3*(2*(-1 + z)*z - 4*(-1 + z)*z*Log[1 - z] + 4*(-1 + z)*z*Log[1 - z]^2 - 
   (4*(1 - 2*z + 2*z^2)*Log[1 - z]^3)/3) + 
 (2*ep^4*(3*(-1 + z)*z - 6*(-1 + z)*z*Log[1 - z] + 6*(-1 + z)*z*Log[1 - z]^2 - 
    4*(-1 + z)*z*Log[1 - z]^3 + (1 - 2*z + 2*z^2)*Log[1 - z]^4))/3,
(* -- *)
calP[qg,R,4][z_] :> 1 - 2*z + 2*z^2 + ep*(-1 + (1 - 2*z + 2*z^2)*(1 - 4*Log[1 - z])) + 
  ep^2*(2*(-1 + z)*z - 8*(-1 + z)*z*Log[1 - z] + 8*(1 - 2*z + 2*z^2)*Log[1 - z]^2) + 
  ep^3*(2*(-1 + z)*z - 8*(-1 + z)*z*Log[1 - z] + 16*(-1 + z)*z*Log[1 - z]^2 - 
   (32*(1 - 2*z + 2*z^2)*Log[1 - z]^3)/3) + 
  ep^4*(2*(-1 + z)*z - 8*(-1 + z)*z*Log[1 - z] + 16*(-1 + z)*z*Log[1 - z]^2 - 
   (64*(-1 + z)*z*Log[1 - z]^3)/3 + (32*(1 - 2*z + 2*z^2)*Log[1 - z]^4)/3),
(* -- *)
calP[gq,R,2][z_] :> (2 - 2*z + z^2)/z + ep*(-z - (2*(2 - 2*z + z^2)*Log[1 - z])/z) + 
 ep^2*(2*z*Log[1 - z] + (2*(2 - 2*z + z^2)*Log[1 - z]^2)/z) + 
 ep^3*(-2*z*Log[1 - z]^2 - (4*(2 - 2*z + z^2)*Log[1 - z]^3)/(3*z)) + 
 ep^4*((4*z*Log[1 - z]^3)/3 + (2*(2 - 2*z + z^2)*Log[1 - z]^4)/(3*z)),
(* -- *)
calP[gq,R,4][z_] :> (2 - 2*z + z^2)/z + ep*(-z - (4*(2 - 2*z + z^2)*Log[1 - z])/z) + 
 (4*ep^2*(z^2*Log[1 - z] + 4*Log[1 - z]^2 - 4*z*Log[1 - z]^2 + 2*z^2*Log[1 - z]^2))/z - 
 (8*ep^3*(3*z^2*Log[1 - z]^2 + 8*Log[1 - z]^3 - 8*z*Log[1 - z]^3 + 4*z^2*Log[1 - z]^3))/
  (3*z) + (32*ep^4*(z^2*Log[1 - z]^3 + 2*Log[1 - z]^4 - 2*z*Log[1 - z]^4 + 
    z^2*Log[1 - z]^4))/(3*z)
};
   
getGa = {ga[qq,R,2,2] -> 3/2 + ep*(13/2 - (2*Pi^2)/3) + 
	ep^2*(26 - Pi^2 - 16*Zeta[3]) + 
	ep^3*(104 - (13*Pi^2)/3 - (2*Pi^4)/5 - 24*Zeta[3]) + 
    ep^4*(416 - (52*Pi^2)/3 - (3*Pi^4)/5 - 104*Zeta[3] + 
   (32*Pi^2*Zeta[3])/3 - 192*Zeta[5])};


(* ::Subsubsection::Closed:: *)
(*Expressions for 1L splitting functions*)


getPRV = {
Pt1L[qq,z_] :> (((1 + 2*ep)*(ep*(1 - z) + z))/(1 - z)^(3*ep) + 
 (2*(-(ep*(1 - z)) + (1 + z^2)/(1 - z))*
   (-(Log[z]/ep) + PolyLog[2, 1 - z] + ep*PolyLog[3, 1 - z]))/
  (1 - z)^(3*ep)),
(* -- *)
Pt1L[qg,z_] :> 
 -(((1 + 2*ep)*(-ep + z))/((1 - ep)*(1 - z)^(3*ep))) - 
 (2*(1 - (2*(1 - z)*z)/(1 - ep))*Gamma[1 - ep]^2*Gamma[1 + ep]^2)/
  (ep^2*(1 - z)^(4*ep)*Gamma[1 - 2*ep]*Gamma[1 + 2*ep]) +  (2*
  (1 - (2*(1 - z)*z)/(1 - ep))*(ep^(-2) - Log[z]/ep - PolyLog[2, 1 - z] + 
    ep*PolyLog[3, 1 - z]))/(1 - z)^(3*ep),
(* -- *)
Pt1L[gq,z_] :> -2 * (1/ep^2+3/2/ep+4+8ep)(1-z)^(-3ep) * ((1+(1-z)^2)/z-ep z)
};


(* ::Subsubsection::Closed:: *)
(*Expressions for AP splitting functions*)


(* extract couplings *)

getPAP = {
PAP0[qq,z_]:>PAP0R[qq,z]+PAP0D[qq]delta[1-z],
PAP0R[qq,z_]:> -1-z + 2DD[0,z], PAP0D[qq]:> 3/2,
PAP0[qg,z_]:>z^2+(1-z)^2, 
PAP0[gq,z_]:> (1+(1-z)^2)/z,
(* -- *)
(* NLO AP below, factor out Cf*q^2 with Cf^2 -> 2*Cf*q^2, 
   i.e. take normal PAP1 and replace Cf^2 -> 2 *)
PAP1[nsMIX,z_]:>2*(3 - 2*z - (2*(1 + z^2)*Log[1 - z]*Log[z])/(1 - z) + 
  (2*(1 + z^2)*PolyLog[2, 1 - z]*(1 - tag[int]))/(1 - z) + (-8 + 7*z)*tag[int] + 
  (Log[z]*(2 - 2*z - 5*tag[int] + 2*z^2*tag[int]))/(1 - z) + 
  (Log[z]^2*(1 + 3*z^2 - 2*(1 + z^2)*tag[int]))/(2*(1 - z)) + 
  delta[1 - z]*(3/8 - Pi^2/2 + 6*Zeta[3])),
(* -- *)
PAP1[qqbVMix,z_]:>2*(4*(1 - z) + 2*(1 + z)*Log[z] + 
  (2*(1 + z^2)*(-Pi^2/6 + Log[z]^2/2 - 2*Log[z]*Log[1 + z] - 2*PolyLog[2, -z]))/(1 + z)),
(* -- *)
PAP1[qgAB,z_] :> ((4 - 9*z + 4*Log[1 - z] + 
	((1 - z)^2 + z^2)*(10 - (2*Pi^2)/3 - 4*Log[(1 - z)/z] + 
      2*Log[(1 - z)/z]^2) - (1 - 4*z)*Log[z] - (1 - 2*z)*Log[z]^2))/2,
(* -- this one is taken from ESW, with Cf^2 -> 1 -- *)
PAP1[gqAB,x_] :> -5/2 -7/2 x + (2+7/2x)Log[x]-(1-1/2x)Log[x]^2 + 
	-2x Log[1-x] - (3Log[1-x]+Log[1-x]^2)(1+(1-x)^2)/x,
(* -- this one is taken from ESW, with Ca->0, Cf->1, tr->1 -- *)
PAP1[agAB,x_] :> -16+8 x + 20/3 x^2 + 4/3/x - (6+10x)Log[x] - (2+2x)Log[x]^2
};


(* ::Subsubsection::Closed:: *)
(*Expressions for convolutions*)


(* taken from the DIS notebook, but divided by Cf^2 *)
(* or obtaied from ~/projects/gen/cv.wl *)


getPoP = {conv[cP[qq,R,2],zm2ep cP[qq,R,2]][z_] :> 
2*(-1 + z) + 8*DD[1, z] - (2*Pi^2*delta[1 - z])/3 - 
 4*(1 + z)*Log[1 - z] + ((1 + 3*z^2)*Log[z])/(-1 + z) + 
 ep^2*(10*(-1 + z) - (4*Pi^2*(1 + z))/3 - (32*Pi^2*DD[1, z])/3 + 
   (112*DD[3, z])/3 - (2*Pi^4*delta[1 - z])/5 - 
   (56*(1 + z)*Log[1 - z]^3)/3 + ((-15 + 2*Pi^2)*(1 + z)*Log[z])/3 - 
   4*Log[z]^2 + (2*(1 + 3*z^2)*Log[z]^3)/(3*(-1 + z)) + 
   Log[1 - z]^2*(4*(-1 + z) - (4*(-3 + z^2)*Log[z])/(-1 + z)) + 
   Log[1 - z]*(8*Pi^2*(1 + z) + 4*(-1 + z)*Log[z] + 
     (4*(1 + z^2)*Log[z]^2)/(-1 + z)) + 
   (8*(1 + z) - 16*(1 + z)*Log[1 - z])*PolyLog[2, z] - 
   16*(1 + z)*PolyLog[3, 1 - z] - 8*(1 + z)*PolyLog[3, z] - 
   24*(1 + z)*Zeta[3] + 64*DD[0, z]*Zeta[3]) + 
 ep*(-2*Pi^2*(1 + z) + (8*Pi^2*DD[0, z])/3 - 24*DD[2, z] + 
   12*(1 + z)*Log[1 - z]^2 + 4*Log[z] + ((1 + 3*z^2)*Log[z]^2)/
    (1 - z) + Log[1 - z]*(4 - 4*z - (4*(1 + z^2)*Log[z])/(-1 + z)) + 
   4*(1 + z)*PolyLog[2, z] - 16*delta[1 - z]*Zeta[3]),
(* -- *)
conv[cP[qq,R,2],zm2ep cP[qg,R,2]][z_] :>
(-2 + 5*z - 3*z^2 + 2*(1 - 2*z + 2*z^2)*Log[1 - z] + 
  (-1 + 2*z - 4*z^2)*Log[z] + ep*((-9 + 57*z + 4*(-12 + Pi^2)*z^2)/3 - 
    6*(1 - 2*z + 2*z^2)*Log[1 - z]^2 + (3 + 6*z - 4*z^2)*Log[z] + 
    (1 - 2*z + 4*z^2)*Log[z]^2 + Log[1 - z]*(8*(1 - 3*z + 2*z^2) + 
      4*(1 - 2*z + 2*z^2)*Log[z]) + 4*(1 - 2*z)*PolyLog[2, z]) + 
  ep^2*((28*(1 - 2*z + 2*z^2)*Log[1 - z]^3)/3 - 
    (2*(3 - 3*z + 6*z^2 + Pi^2*(-1 + 2*z))*Log[z])/3 + 
    (-3 - 6*z + 4*z^2)*Log[z]^2 - (2*(1 - 2*z + 4*z^2)*Log[z]^3)/3 + 
    Log[1 - z]^2*(-4*(4 - 13*z + 9*z^2) - 4*(3 - 6*z + 2*z^2)*Log[z]) + 
    Log[1 - z]*((-4*(-9 + 60*z + (-51 + 4*Pi^2)*z^2))/3 + 
      8*(-1 + z)*z*Log[z] - 4*(1 - 2*z + 2*z^2)*Log[z]^2) + 
    (4*(3 + 2*z) + 16*(-1 + 2*z)*Log[1 - z])*PolyLog[2, z] + 
    16*(-1 + 2*z)*PolyLog[3, 1 - z] + 8*(-1 + 2*z)*PolyLog[3, z] + 
    (2*(Pi^2*(-1 - 9*z + 5*z^2) + 6*(-4 + 6*Zeta[3] - 
         3*z*(-7 + 4*Zeta[3]) + z^2*(-17 + 8*Zeta[3]))))/3)),
(* -- *)
conv[cP[gq,R,2],zm2ep cP[qq,R,2]][z_] :> 
((-3 + 5*z - 2*z^2)/z + (2*(2 - 2*z + z^2)*Log[1 - z])/z - (-2 + z)*Log[z] + 
  ep*(((39 + 4*Pi^2*(-1 + z) - 3*z)*(-1 + z))/(3*z) + (2*(6 - 10*z + 3*z^2)*Log[1 - z])/
     z - (6*(2 - 2*z + z^2)*Log[1 - z]^2)/z - 2*(2 + z)*Log[z] + (-2 + z)*Log[z]^2 - 
    4*(-2 + z)*PolyLog[2, z]) + 
  ep^2*((-4*(39 + 4*Pi^2*(-1 + z) - 3*z)*(-1 + z)*Log[1 - z])/(3*z) + 
    (28*(2 - 2*z + z^2)*Log[1 - z]^3)/(3*z) + ((-2*Pi^2*(-2 + z) + 3*(4 + 5*z))*Log[z])/
     3 + 2*(2 + z)*Log[z]^2 - (2*(-2 + z)*Log[z]^3)/3 + 
    Log[1 - z]^2*((-2*(12 - 20*z + 5*z^2))/z + 8*(-2 + z)*Log[z]) + 
    (-8*(2 + z) + 16*(-2 + z)*Log[1 - z])*PolyLog[2, z] + 
    16*(-2 + z)*PolyLog[3, 1 - z] + 8*(-2 + z)*PolyLog[3, z] + 
    (Pi^2*(6 - 2*z + 6*z^2) + 3*(-52 + z*(65 - 16*Zeta[3]) + 32*Zeta[3] + 
        z^2*(-13 + 8*Zeta[3])))/(3*z))),
(* -- *)
conv[cP[gq,R,2],zm2ep cP[qg,R,2]][z_] :> 
1 + 4/(3*z) - z - (4*z^2)/3 + 2*(1 + z)*Log[z] + 
 ep*((-3 + 46/z + 3*z - 46*z^2 - 12*Pi^2*(1 + z))/9 + 
   (4*(-3 - 4/z + 3*z + 4*z^2)*Log[1 - z])/3 - 2*(1 + z)*Log[z] - 2*(1 + z)*Log[z]^2 + 
   8*(1 + z)*PolyLog[2, z]) + 
 ep^2*((4*(-46 + 3*(1 + 4*Pi^2)*z + 3*(-1 + 4*Pi^2)*z^2 + 46*z^3)*Log[1 - z])/(9*z) + 
   (2*(3 + 2*Pi^2)*(1 + z)*Log[z])/3 + 2*(1 + z)*Log[z]^2 + (4*(1 + z)*Log[z]^3)/3 + 
   Log[1 - z]^2*((8*(3 + 4/z - 3*z - 4*z^2))/3 - 16*(1 + z)*Log[z]) + 
   (-8*(1 + z) - 32*(1 + z)*Log[1 - z])*PolyLog[2, z] - 32*(1 + z)*PolyLog[3, 1 - z] - 
   16*(1 + z)*PolyLog[3, z] + (562 - 562*z^3 + 6*Pi^2*(-4 + 3*z + 9*z^2 + 4*z^3) + 
     6*z^2*(-7 + 72*Zeta[3]) + 6*z*(7 + 72*Zeta[3]))/(27*z))};


getPoPAP = {conv[cP[qq,R,2],zm2ep cAP0[qq]][z_] :> 
2*(-1 + z) + 8*DD[1, z] - (2*Pi^2*delta[1 - z])/3 - 4*(1 + z)*Log[1 - z] + 
 ((1 + 3*z^2)*Log[z])/(-1 + z) + ep*(-(Pi^2*(1 + z)) + (4*Pi^2*DD[0, z])/3 - 
   12*DD[2, z] - 2*(-1 + z)*Log[1 - z] + 6*(1 + z)*Log[1 - z]^2 - 
   (-3 + z)*Log[z] - ((1 + 3*z^2)*Log[z]^2)/(-1 + z) + 
   2*(1 + z)*PolyLog[2, z] - 8*delta[1 - z]*Zeta[3]) + 
 ep^2*((-8*Pi^2*DD[1, z])/3 + (32*DD[3, z])/3 - (8*Pi^4*delta[1 - z])/45 + 
   2*Pi^2*(1 + z)*Log[1 - z] - (16*(1 + z)*Log[1 - z]^3)/3 + 
   (2*(-3 + Pi^2*(1 + z))*Log[z])/3 - 2*(1 + z)*Log[1 - z]^2*Log[z] + 
   (-3 + z)*Log[z]^2 + (2*(1 + 3*z^2)*Log[z]^3)/(3*(-1 + z)) + 
   (-2*(-3 + z) - 4*(1 + z)*Log[1 - z])*PolyLog[2, z] - 
   4*(1 + z)*PolyLog[3, 1 - z] - 4*(1 + z)*PolyLog[3, z] + 
   16*DD[0, z]*Zeta[3] + (Pi^2*(-5 + 3*z) - 
     6*(1 + 2*Zeta[3] + z*(-1 + 2*Zeta[3])))/3) + (3*calP[qq, R, 2][z])/2,
(* -- *)
conv[cP[qq,R,2],zm2ep cAP0[qg]][z_] :> 
-2 + 5*z - 3*z^2 + 2*(1 - 2*z + 2*z^2)*Log[1 - z] + 
 (-1 + 2*z - 4*z^2)*Log[z] + 
 ep*((Pi^2*(-1 + 2*z) - 3*(1 - 9*z + 8*z^2))/3 + 
   2*(2 - 5*z + 3*z^2)*Log[1 - z] - 2*(1 - 2*z + 2*z^2)*Log[1 - z]^2 + 
   (3 + 2*z)*Log[z] + (1 - 2*z + 4*z^2)*Log[z]^2 + 
   2*(1 - 2*z)*PolyLog[2, z]) + 
 ep^2*((2*(3 + Pi^2 - 27*z - 2*Pi^2*z + 24*z^2)*Log[1 - z])/3 + 
   (4*(1 - 2*z + 2*z^2)*Log[1 - z]^3)/3 + 
   (2*(-3 + Pi^2*(1 - 2*z))*Log[z])/3 - (3 + 2*z)*Log[z]^2 - 
   (2*(1 - 2*z + 4*z^2)*Log[z]^3)/3 + Log[1 - z]^2*
    (-2*(2 - 5*z + 3*z^2) + 2*(-1 + 2*z)*Log[z]) + 
   (2*(3 + 2*z) + 4*(-1 + 2*z)*Log[1 - z])*PolyLog[2, z] + 
   4*(-1 + 2*z)*PolyLog[3, 1 - z] + 4*(-1 + 2*z)*PolyLog[3, z] + 
   (-(Pi^2*(3 + 2*z)) - 12*(1 + 4*z^2 - Zeta[3] + z*(-5 + 2*Zeta[3])))/3),
(* -- *)
conv[cP[gq,R,2],zm2ep cAP0[qq]][z_] :> 
(-3 + 5*z - 2*z^2)/z + (2*(2 - 2*z + z^2)*Log[1 - z])/z - 
 (-2 + z)*Log[z] + ep*((15*(-1 + z) + Pi^2*(4 - 6*z + 3*z^2))/
    (3*z) + (2*(3 - 5*z + z^2)*Log[1 - z])/z - 
   (4*(2 - 2*z + z^2)*Log[1 - z]^2)/z + (-2 + z)*Log[z] + 
   (-2 + z)*Log[z]^2 - 2*(-2 + z)*PolyLog[2, z]) + 
 ep^2*((-2*(15*(-1 + z) + Pi^2*(4 - 6*z + 3*z^2))*Log[1 - z])/
    (3*z) + (4*(2 - 2*z + z^2)*Log[1 - z]^3)/z - 
   (2*(-3 + Pi^2*(-2 + z))*Log[z])/3 - (-2 + z)*Log[z]^2 - 
   (2*(-2 + z)*Log[z]^3)/3 + Log[1 - z]^2*((2*(-3 + 5*z))/z + 
     2*(-2 + z)*Log[z]) + (2*(-2 + z) + 4*(-2 + z)*Log[1 - z])*
    PolyLog[2, z] + 4*(-2 + z)*PolyLog[3, 1 - z] + 
   4*(-2 + z)*PolyLog[3, z] + (-27 + 2*z*(15 + Pi^2 - 12*Zeta[3]) - 
     3*z^2*(1 + Pi^2 - 4*Zeta[3]) + 48*Zeta[3])/(3*z)) + 3/2 calP[gq,R,2][z],
(* -- *)
conv[cP[gq,R,2],zm2ep cAP0[qg]][z_] :> 
1 + 4/(3*z) - z - (4*z^2)/3 + 2*(1 + z)*Log[z] + 
 ep*((-15 + 26/z - 3*z - 8*z^2 - 6*Pi^2*(1 + z))/9 + (-2 - 8/(3*z) + 2*z + (8*z^2)/3)*
    Log[1 - z] - 2*(1 + z)*Log[z] - 2*(1 + z)*Log[z]^2 + 4*(1 + z)*PolyLog[2, z]) + 
 ep^2*((2*(15 - 26/z + 3*z + 8*z^2 + 6*Pi^2*(1 + z))*Log[1 - z])/9 + 
   (2 + (4*Pi^2*(1 + z))/3)*Log[z] + 2*(1 + z)*Log[z]^2 + (4*(1 + z)*Log[z]^3)/3 + 
   Log[1 - z]^2*(2 + 8/(3*z) - 2*z - (8*z^2)/3 - 4*(1 + z)*Log[z]) + 
   (-4*(1 + z) - 8*(1 + z)*Log[1 - z])*PolyLog[2, z] - 8*(1 + z)*PolyLog[3, 1 - z] - 
   8*(1 + z)*PolyLog[3, z] + (2*(-51 + 89/z - 35*z^2 + 9*Pi^2*(1 + z) + 108*Zeta[3] + 
      3*z*(-1 + 36*Zeta[3])))/27)
};


getPAPoPAP = {conv[cAP0[qq],cAP0[qq]][z_] :> 
	-5 - z + 6*DD[0, z] + 8*DD[1, z] + (9/4 - (2*Pi^2)/3)*delta[1 - z] +
	- 4*(1 + z)*Log[1 - z] + ((1 + 3*z^2)*Log[z])/(-1 + z),
(* -- *)
conv[cAP0[qq],cAP0[qg]][z_] :> (-2 + 5*z - 3*z^2) + 
	2*(1 - 2*z + 2*z^2)*Log[1 - z] - (1 - 2*z + 4*z^2)*Log[z] + 3/2 PAP0[qg,z],
(* -- *)
conv[cAP0[gq],cAP0[qq]][z_] :> (-3 + 5*z - 2*z^2)/z + (2*(2 - 2*z + z^2)*Log[1 - z])/z - (-2 + z)*Log[z] +
	3/2 PAP0[gq,z],
(* -- *)
conv[cAP0[gq],cAP0[qg]][z_] :> 1 + 4/(3*z) - z - (4*z^2)/3 + 2*(1 + z)*Log[z]
};


(* ::Subsubsection::Closed:: *)
(*Expressions for Max's integrals*)


(* taken verbatim from Max's ancillary file, with Ca -> 0, Cf -> 1 and multiplied by 2 where needed *)
(* Max's results don't contain any (1-ep) factors from mismatched d-dim averages *)


getMax = 
{
(* PqqID *)
twoISRab[z_,5] :> -2*2*(-1/3*((Pi^2*(1 + z^2) + 3*(8 - 15*z + 7*z^2))*Log[2])/(-1 + z) - 
  ((Pi^2*(1 + z^2) + 3*(8 - 15*z + 7*z^2))*Log[1 - z])/(3*(-1 + z)) - 
  ((66 - 57*z - 39*z^2 + 4*Pi^2*(1 + z^2))*Log[z])/(12*(-1 + z)) + 
  ((-5 + 2*z^2)*Log[2]*Log[z])/(-1 + z) + (3*(-1 + z)*Log[1 - z]*Log[z])/2 + 
  (2*(1 + z^2)*Log[2]*Log[1 - z]*Log[z])/(-1 + z) + (2*(1 + z^2)*Log[1 - z]^2*Log[z])/
   (-1 + z) + ((-25 + 12*z + 4*z^2)*Log[z]^2)/(8*(-1 + z)) - 
  ((1 + z^2)*Log[2]*Log[z]^2)/(-1 + z) - ((1 + z^2)*Log[1 - z]*Log[z]^2)/(4*(-1 + z)) - 
  (7*(1 + z^2)*Log[z]^3)/(12*(-1 + z)) - 2*(1 + z)*Log[z]*Log[1 + z] - 
  2*(1 + z)*PolyLog[2, -z] - (2*(1 + z^2)*Log[z]*PolyLog[2, -z])/(-1 + z) - 
  ((-13 + 6*z + z^2)*PolyLog[2, z])/(2*(-1 + z)) + (2*(1 + z^2)*Log[2]*PolyLog[2, z])/
   (-1 + z) + (2*(1 + z^2)*Log[1 - z]*PolyLog[2, z])/(-1 + z) - 
  (3*(1 + z^2)*Log[z]*PolyLog[2, z])/(2*(-1 + z)) + 
  ((Pi^2*(1 + z^2) + 3*(8 - 15*z + 7*z^2))/(12*(-1 + z)) - 
    ((-5 + 2*z^2)*Log[z])/(4*(-1 + z)) - ((1 + z^2)*Log[1 - z]*Log[z])/(2*(-1 + z)) + 
    ((1 + z^2)*Log[z]^2)/(4*(-1 + z)) - ((1 + z^2)*PolyLog[2, z])/(2*(-1 + z)))/ep + 
  (3*(1 + z^2)*PolyLog[3, 1 - z])/(-1 + z) + (4*(1 + z^2)*PolyLog[3, -z])/(-1 + z) + 
  (9*(1 + z^2)*PolyLog[3, z])/(2*(-1 + z)) - 
  (Pi^2*(11 - 6*z + z^2) + 3*(-1 - 6*z + 6*Zeta[3] + z^2*(7 + 6*Zeta[3])))/(12*(-1 + z))),
(* PqqbID *)
twoISRab[z_,6] :> -2*2*(((12*(-1 + z^2) + Pi^2*(1 + z^2))*Log[2])/(3*(1 + z)) + 
  ((12*(-1 + z^2) + Pi^2*(1 + z^2))*Log[1 - z])/(3*(1 + z)) - 
  ((57 + 90*z + 33*z^2 + 4*Pi^2*(1 + z^2))*Log[z])/(12*(1 + z)) - 
  2*(1 + z)*Log[2]*Log[z] - (-1 + z)*Log[1 - z]*Log[z] - (2 + z)*Log[z]^2 - 
  ((1 + z^2)*Log[2]*Log[z]^2)/(1 + z) - (7*(1 + z^2)*Log[z]^3)/(12*(1 + z)) + 
  (Pi^2*(1 + z^2)*Log[1 + z])/(2*(1 + z)) + 3*(1 + z)*Log[z]*Log[1 + z] + 
  (4*(1 + z^2)*Log[2]*Log[z]*Log[1 + z])/(1 + z) + 
  (4*(1 + z^2)*Log[1 - z]*Log[z]*Log[1 + z])/(1 + z) + 
  ((1 + z^2)*Log[z]^2*Log[1 + z])/(2*(1 + z)) + (3*(1 + z^2)*Log[z]*Log[1 + z]^2)/
   (1 + z) - ((1 + z^2)*Log[1 + z]^3)/(1 + z) + 3*(1 + z)*PolyLog[2, -z] + 
  (4*(1 + z^2)*Log[2]*PolyLog[2, -z])/(1 + z) + (4*(1 + z^2)*Log[1 - z]*PolyLog[2, -z])/
   (1 + z) - ((1 + z^2)*Log[z]*PolyLog[2, -z])/(1 + z) + 
  ((-12*(-1 + z^2) - Pi^2*(1 + z^2))/(12*(1 + z)) + ((1 + z)*Log[z])/2 + 
    ((1 + z^2)*Log[z]^2)/(4*(1 + z)) - ((1 + z^2)*Log[z]*Log[1 + z])/(1 + z) - 
    ((1 + z^2)*PolyLog[2, -z])/(1 + z))/ep + (3 + z)*PolyLog[2, z] - 
  ((1 + z^2)*Log[z]*PolyLog[2, z])/(1 + z) - (4*(1 + z^2)*PolyLog[3, 1 - z])/(1 + z) + 
  (9*(1 + z^2)*PolyLog[3, -z])/(1 + z) + (4*(1 + z^2)*PolyLog[3, z])/(1 + z) + 
  (6*(1 + z^2)*PolyLog[3, z/(1 + z)])/(1 + z) + (2*(1 + z^2)*PolyLog[3, 1 - z^2])/
   (1 + z) + (Pi^2*(-3 - 2*z + z^2) - 15*(3 + 2*Zeta[3] + z^2*(-3 + 2*Zeta[3])))/
   (12*(1 + z))),
(* Pgg, qqg emission, Cf->1, tr->1 *)
ISab[z_,7] :> 2*(15 - (Pi^2*(1 + 3*z))/3 + 10*(-1 + z)*Log[1 - z] + (2*(3 + Pi^2)*(1 + z)*Log[z])/3 + 
  ((-1 - 7*z)*Log[z]^2)/4 + (5*(1 + z)*Log[z]^3)/6 + 
  Log[2]^2*((-7*(-1 + z)*(4 + 7*z + 4*z^2))/(12*z) + (7*(1 + z)*Log[z])/2) + 
  Log[2]*((-46 + 3*(-59 + 4*Pi^2)*z + 3*(59 + 4*Pi^2)*z^2 + 46*z^3)/(18*z) + 
    (2 + 8/(3*z) - 2*z - (8*z^2)/3)*Log[1 - z] + (-1 - 5*z)*Log[z] + 
    3*(1 + z)*Log[z]^2) + ((-5*(-1 + z))/2 + ((1 + 3*z)*Log[z])/2 + 
    ((-1 - z)*Log[z]^2)/2 + Log[2]*(((-1 + z)*(4 + 7*z + 4*z^2))/(6*z) + 
      (-1 - z)*Log[z]))/ep + (2 + 6*z - 4*(1 + z)*Log[2])*PolyLog[2, z] - 
  4*(1 + z)*PolyLog[3, z] + 4*Zeta[3] + z*(-15 + 4*Zeta[3])),
(* Pgq *)
ISab[z_,8,e_] :> (3*(2 - 2*z + z^2)*Log[2]*Log[em/e]^2)/z + 
 ((9*z*(1 + z) + 2*Pi^2*(2 - 2*z + z^2))*Log[1 - z])/(3*z) + 
 (Log[em/e]*(-4*Pi^2*(2 - 2*z + z^2) + 3*Log[2]*(-12 + 12*z + 38*Log[2] - 38*z*Log[2] + 
      19*z^2*Log[2]) + 36*(2 - 2*z + z^2)*Log[2]*Log[1 - z]))/(6*z) + 
 (7/4 + (Pi^2*(-2 + z))/3 + z)*Log[z] + (-2 - (3*z)/8)*Log[z]^2 + 
 (5*(-2 + z)*Log[z]^3)/12 + Log[2]^2*(-71/4 + 57/(4*z) + 8*z + 
   (19 - 19/z - (19*z)/2)*Log[1 - z] + (7*(-2 + z)*Log[z])/4) + 
 Log[2]*((4*Pi^2*(3 - 4*z + 2*z^2) - 3*(48 - 57*z + 10*z^2))/(6*z) + 
   (-28 + 24/z + 10*z)*Log[1 - z] + (18 - 18/z - 9*z)*Log[1 - z]^2 + (-5 - z/2)*Log[z] + 
   (3*(-2 + z)*Log[z]^2)/2) + ((-3*(1 + z))/4 - (3*(2 - 2*z + z^2)*Log[2]*Log[em/e])/z + 
   (1 - z/4)*Log[z] + ((2 - z)*Log[z]^2)/4 + 
   Log[2]*((11 - 9/z - 5*z)/2 + (-6 + 6/z + 3*z)*Log[1 - z] + (1 - z/2)*Log[z]))/ep + 
 (4 - z + (4 - 2*z)*Log[2])*PolyLog[2, z] + (4 - 2*z)*PolyLog[3, z] + 
 (-4*Pi^2*(3 - z + z^2) + 3*z*(-19 - 16*Zeta[3] + z*(-19 + 8*Zeta[3])))/(12*z),
(* Pqg *)
ISab[z_,9,e1_] :> Pi^2*(-1/3 + z - z^2) + 3*(1 - 2*z + 2*z^2)*Log[2]*Log[em/e1]^2 + 
 ((3 + 15*z + Pi^2*(-2 + 4*z - 4*z^2))*Log[1 - z])/3 + 
 (Log[em/e1]*(Pi^2*(-4 + 8*z - 8*z^2) + 3*Log[2]*(6 + 19*Log[2] - 38*z*Log[2] + 
      38*z^2*Log[2]) + 36*(1 - 2*z + 2*z^2)*Log[2]*Log[1 - z]))/6 + 
 ((9 - 21*z - 4*Pi^2*(5 - 10*z + 8*z^2))/12 + Log[1 - z] + 
   (4 - 8*z + 8*z^2)*Log[1 - z]^2)*Log[z] + (9/8 + z + (1/2 - z + z^2)*Log[1 - z])*
  Log[z]^2 + ((-19 + 38*z - 28*z^2)*Log[z]^3)/12 + 
 Log[2]^2*(8 - (71*z)/4 + (57*z^2)/4 - (19*(1 - 2*z + 2*z^2)*Log[1 - z])/2 + 
   (7*(-1 + 2*z)*Log[z])/4) + Log[2]*((-4*Pi^2*(-1 + z)^2 - 3*(9 - 50*z + 39*z^2))/6 + 
   (7 - 22*z + 18*z^2)*Log[1 - z] - 9*(1 - 2*z + 2*z^2)*Log[1 - z]^2 + 
   (2*(2 + z) + (4 - 8*z + 8*z^2)*Log[1 - z])*Log[z] + (-7/2 + 7*z - 4*z^2)*Log[z]^2) + 
 (-2 + (6 - 12*z + 8*z^2)*Log[2] + (4 - 8*z + 8*z^2)*Log[1 - z] + 
   (-1 + 2*z - 2*z^2)*Log[z])*PolyLog[2, z] + 
 ((-3*(1 + 5*z) + Pi^2*(4 - 8*z + 8*z^2))/12 - 3*(1 - 2*z + 2*z^2)*Log[2]*Log[em/e1] + 
   (-3/4 + (-1 + 2*z - 2*z^2)*Log[1 - z])*Log[z] + (3/4 - (3*z)/2 + z^2)*Log[z]^2 + 
   Log[2]*((-5 + 11*z - 9*z^2)/2 + (3 - 6*z + 6*z^2)*Log[1 - z] + (1/2 - z)*Log[z]) + 
   (-1 + 2*z - 2*z^2)*PolyLog[2, z])/ep + (4 - 8*z + 8*z^2)*PolyLog[3, 1 - z] + 
 (9 - 18*z + 14*z^2)*PolyLog[3, z] + 4*z*(-2 + Zeta[3]) - 2*Zeta[3]
};


(* ::Subsubsection::Closed:: *)
(*Soft integrals*)


sij[eta_] := Gamma[1-ep]^2/Gamma[1-2ep] eta^(-ep) * kbar[eta];
sij[1] = Gamma[1-ep]^2/Gamma[1-2ep];

kbar[x_] := 1 + ep^2*PolyLog[2, 1 - x] + 
 ep^4*(-(Log[x]^2*(4*Pi^2 - 12*Log[1 - x]*Log[x] + Log[x]^2))/24 + 
   (Log[x]^2*PolyLog[2, 1 - x])/2 - PolyLog[4, (-1 + x)/x] + 
    Log[x]*(PolyLog[3, 1 - x] + PolyLog[3, x] - Zeta[3])) + 
 ep^3*(-(Log[x]*(Pi^2 - 3*Log[1 - x]*Log[x]))/6 + Log[x]*PolyLog[2, 1 - x] + 
   PolyLog[3, 1 - x] + PolyLog[3, x] - Zeta[3]) + ep^4 eporder[kbar,4];


(* ::Subsubsection::Closed:: *)
(*Fix delta functions, ln(2), ln(z) and Li2(z)*)


fixDelta = {Log[z_]^a_ delta[1-z_]:>0, 
	Log[z_] delta[1-z_]:>0,
	DD[a_,z_]Log[z_]^k_ :> Log[z]^k Log[1-z]^a/(1-z),
	DD[a_,z_]Log[z_] :> Log[z] Log[1-z]^a/(1-z),
	delta[1-z_]FLM[z_,a___]:>FLM[1,a],
	delta[1-z_] FLM[a_,z_,b___]:>FLM[a,2,b]};
	
fixDeltaCH = {Log[z_]^a_ delta[1-z_]:>0, 
	Log[z_] delta[1-z_]:>0,
	DD[a_,z_]Log[z_]^k_ :> Log[z]^k Log[1-z]^a/(1-z),
	DD[a_,z_]Log[z_] :> Log[z] Log[1-z]^a/(1-z),
	delta[1-z_]FLM[ch___][z_,a___]:>FLM[ch][1,a],
	delta[1-z_] FLM[ch___][a_,z_,b___]:>FLM[ch][a,2,b]};
	
fixLog2 = {Log[4]->2Log[2],Log[8]->3Log[2],Log[16]->4Log[2]};

fixLi2 = {PolyLog[2,z]->-PolyLog[2,1-z]-Log[1-z]Log[z]+Zeta[2],
	PolyLog[2,zb]->-PolyLog[2,1-zb]-Log[1-zb]Log[zb]+Zeta[2]};

fixLogz = {Log[z^a_]:>a Log[z],Log[-1+1/z]->Log[1-z]-Log[z]};


(* ::Subsubsection::Closed:: *)
(*Prefactors*)


getPref = {Nb -> Gamma[1-ep]Gamma[1-2ep]/Gamma[1-3ep],
Nc -> Gamma[1-ep]Gamma[1+2ep]/Gamma[1+ep]};


(* ::Subsubsection::Closed:: *)
(*Substitutions from Chiara's to my notation*)


fromChiara = {\[Alpha]s[\[Mu]] -> qas 2 Pi, \[Alpha][\[Mu]] -> qae 2Pi,  \[Alpha] -> qae, \[Alpha]s -> qas, CF -> Cf, TR -> tr, 
	Nc -> xn, Emax -> em,
	(* -- *)
	\[Mu]->Sqrt[mu2],
	D0[z_]:>DD[0,z],D1[z_]:>DD[1,z],D2[z_]:>DD[2,z],D3[z_]:>DD[3,z],
	eta12->eta[1,2],eta13->eta[1,3],eta14->eta[1,4],eta15->eta[1,5],
	eta23->eta[2,3],eta24->eta[2,4],eta25->eta[2,5],eta34->eta[3,4],eta35->eta[3,5],
	eta45->eta[4,5],
	z1->z,z2->zb,
	E1->e[1],E2->e[2],E3->e[3],E4->e[4],E5->e[5],EC->Sqrt[s]/2,
	s13->4e[1]e[3]eta[1,3],s14->4e[1]e[4]eta[1,4],s23->4e[2]e[3]eta[2,3],
	s24->4e[2]e[4]eta[2,4],	s35 -> 4 e[3]e[5]eta[3,5],s45->4 e[4] e[5] eta[4,5],
	Qe->Ql};


forChiara = {e[1] -> Sqrt[s]/2, e[2] -> Sqrt[s]/2, eta[1,2]->1,em->Sqrt[s]/2};


(* ::Subsection:: *)
(*Channel-dependent substitutions*)


(* ::Subsubsection::Closed:: *)
(*Max's non-singlet results, taken from Kirill*)


getMaxNS={maxresNS -> (-2)* mu2^(2ep)*
	((Cf*qae*qas*q1^2*FLM[z,2,3,4]*(-18 - 66*ep + 36*z + 132*ep*z + 
      8*ep*Pi^2*z - 18*z^2 - 66*ep*z^2 - 8*ep*Pi^2*z^2 + 72*ep*Log[2] - 
      144*ep*z*Log[2] + 72*ep*z^2*Log[2] + 72*ep*Log[1 - z] - 
      144*ep*z*Log[1 - z] + 72*ep*z^2*Log[1 - z] + 20*ep*Pi^2*Log[z] - 
      12*z*Log[z] + 12*ep*z*Log[z] + 12*z^2*Log[z] - 12*ep*z^2*Log[z] + 
      12*ep*Pi^2*z^2*Log[z] + 48*ep*z*Log[2]*Log[z] - 
      48*ep*z^2*Log[2]*Log[z] - 9*Log[z]^2 + 6*ep*z*Log[z]^2 - 
      3*z^2*Log[z]^2 - 6*ep*z^2*Log[z]^2 + 36*ep*Log[2]*Log[z]^2 + 
      12*ep*z^2*Log[2]*Log[z]^2 + 19*ep*Log[z]^3 + 9*ep*z^2*Log[z]^3 + 
      24*ep*(2*(-1 + z)*z + (1 + z^2)*Log[z])*PolyLog[2, z] - 
      24*ep*(5 + 3*z^2)*PolyLog[3, z] + 120*ep*Zeta[3] + 
      72*ep*z^2*Zeta[3]))/(12*e[1]^(4*ep)*ep*(-1 + z)) + 
   (Cf*qae*qas*q2^2*FLM[1,zb,3,4]*(-18 - 66*ep + 36*zb + 132*ep*zb + 
      8*ep*Pi^2*zb - 18*zb^2 - 66*ep*zb^2 - 8*ep*Pi^2*zb^2 + 
      72*ep*Log[2] - 144*ep*zb*Log[2] + 72*ep*zb^2*Log[2] + 
      72*ep*Log[1 - zb] - 144*ep*zb*Log[1 - zb] + 
      72*ep*zb^2*Log[1 - zb] + 20*ep*Pi^2*Log[zb] - 12*zb*Log[zb] + 
      12*ep*zb*Log[zb] + 12*zb^2*Log[zb] - 12*ep*zb^2*Log[zb] + 
      12*ep*Pi^2*zb^2*Log[zb] + 48*ep*zb*Log[2]*Log[zb] - 
      48*ep*zb^2*Log[2]*Log[zb] - 9*Log[zb]^2 + 6*ep*zb*Log[zb]^2 - 
      3*zb^2*Log[zb]^2 - 6*ep*zb^2*Log[zb]^2 + 36*ep*Log[2]*Log[zb]^2 + 
      12*ep*zb^2*Log[2]*Log[zb]^2 + 19*ep*Log[zb]^3 + 
      9*ep*zb^2*Log[zb]^3 + 24*ep*(2*(-1 + zb)*zb + (1 + zb^2)*Log[zb])*
       PolyLog[2, zb] - 24*ep*(5 + 3*zb^2)*PolyLog[3, zb] + 
      120*ep*Zeta[3] + 72*ep*zb^2*Zeta[3]))/(12*e[2]^(4*ep)*ep*(-1 + zb)))+
      (* --- *)
      tag[int] * Cf(2q1^2)*qas qae *FLM[z,2,3,4]*((-24 - Pi^2 + 45*z - 21*z^2 - Pi^2*z^2 + 
    3*(-5 + 2*z^2 + 2*(1 + z^2)*Log[1 - z])*Log[z] - 3*(1 + z^2)*Log[z]^2 + 
    6*(1 + z^2)*PolyLog[2, z])/(6*ep*(e[1]^2/musq)^(2*ep)*(-1 + z)) + 
  (-6 + 22*Pi^2 - 36*z - 12*Pi^2*z + 42*z^2 + 2*Pi^2*z^2 + 192*Log[2] + 8*Pi^2*Log[2] - 
    360*z*Log[2] + 168*z^2*Log[2] + 8*Pi^2*z^2*Log[2] + 192*Log[1 - z] + 
    8*Pi^2*Log[1 - z] - 360*z*Log[1 - z] + 168*z^2*Log[1 - z] + 8*Pi^2*z^2*Log[1 - z] + 
    132*Log[z] + 8*Pi^2*Log[z] - 114*z*Log[z] - 78*z^2*Log[z] + 8*Pi^2*z^2*Log[z] + 
    120*Log[2]*Log[z] - 48*z^2*Log[2]*Log[z] - 36*Log[1 - z]*Log[z] + 
    72*z*Log[1 - z]*Log[z] - 36*z^2*Log[1 - z]*Log[z] - 48*Log[2]*Log[1 - z]*Log[z] - 
    48*z^2*Log[2]*Log[1 - z]*Log[z] - 48*Log[1 - z]^2*Log[z] - 
    48*z^2*Log[1 - z]^2*Log[z] + 75*Log[z]^2 - 36*z*Log[z]^2 - 12*z^2*Log[z]^2 + 
    24*Log[2]*Log[z]^2 + 24*z^2*Log[2]*Log[z]^2 + 6*Log[1 - z]*Log[z]^2 + 
    6*z^2*Log[1 - z]*Log[z]^2 + 14*Log[z]^3 + 14*z^2*Log[z]^3 - 48*Log[z]*Log[1 + z] + 
    48*z^2*Log[z]*Log[1 + z] + 48*(-1 + z^2 + (1 + z^2)*Log[z])*PolyLog[2, -z] - 
    12*(13 - 6*z - z^2 + Log[16] + z^2*Log[16] + 4*(1 + z^2)*Log[1 - z] - 
      3*(1 + z^2)*Log[z])*PolyLog[2, z] - 72*PolyLog[3, 1 - z] - 
    72*z^2*PolyLog[3, 1 - z] - 96*PolyLog[3, -z] - 96*z^2*PolyLog[3, -z] - 
    108*PolyLog[3, z] - 108*z^2*PolyLog[3, z] + 36*Zeta[3] + 36*z^2*Zeta[3])/
   (12*(e[1]^2/musq)^(2*ep)*(-1 + z)))+
   tag[int] * Cf*(2q2^2)*qas qae*FLM[1,zb,3,4]*((-24 - Pi^2 + 45*zb - 21*zb^2 - Pi^2*zb^2 + 
    3*(-5 + 2*zb^2 + 2*(1 + zb^2)*Log[1 - zb])*Log[zb] - 3*(1 + zb^2)*Log[zb]^2 + 
    6*(1 + zb^2)*PolyLog[2, zb])/(6*ep*(e[2]^2/musq)^(2*ep)*(-1 + zb)) + 
  (-6 + 22*Pi^2 - 36*zb - 12*Pi^2*zb + 42*zb^2 + 2*Pi^2*zb^2 + 192*Log[2] + 
    8*Pi^2*Log[2] - 360*zb*Log[2] + 168*zb^2*Log[2] + 8*Pi^2*zb^2*Log[2] + 
    192*Log[1 - zb] + 8*Pi^2*Log[1 - zb] - 360*zb*Log[1 - zb] + 168*zb^2*Log[1 - zb] + 
    8*Pi^2*zb^2*Log[1 - zb] + 132*Log[zb] + 8*Pi^2*Log[zb] - 114*zb*Log[zb] - 
    78*zb^2*Log[zb] + 8*Pi^2*zb^2*Log[zb] + 120*Log[2]*Log[zb] - 
    48*zb^2*Log[2]*Log[zb] - 36*Log[1 - zb]*Log[zb] + 72*zb*Log[1 - zb]*Log[zb] - 
    36*zb^2*Log[1 - zb]*Log[zb] - 48*Log[2]*Log[1 - zb]*Log[zb] - 
    48*zb^2*Log[2]*Log[1 - zb]*Log[zb] - 48*Log[1 - zb]^2*Log[zb] - 
    48*zb^2*Log[1 - zb]^2*Log[zb] + 75*Log[zb]^2 - 36*zb*Log[zb]^2 - 12*zb^2*Log[zb]^2 + 
    24*Log[2]*Log[zb]^2 + 24*zb^2*Log[2]*Log[zb]^2 + 6*Log[1 - zb]*Log[zb]^2 + 
    6*zb^2*Log[1 - zb]*Log[zb]^2 + 14*Log[zb]^3 + 14*zb^2*Log[zb]^3 - 
    48*Log[zb]*Log[1 + zb] + 48*zb^2*Log[zb]*Log[1 + zb] + 
    48*(-1 + zb^2 + (1 + zb^2)*Log[zb])*PolyLog[2, -zb] - 
    12*(13 - 6*zb - zb^2 + Log[16] + zb^2*Log[16] + 4*(1 + zb^2)*Log[1 - zb] - 
      3*(1 + zb^2)*Log[zb])*PolyLog[2, zb] - 72*PolyLog[3, 1 - zb] - 
    72*zb^2*PolyLog[3, 1 - zb] - 96*PolyLog[3, -zb] - 96*zb^2*PolyLog[3, -zb] - 
    108*PolyLog[3, zb] - 108*zb^2*PolyLog[3, zb] + 36*Zeta[3] + 36*zb^2*Zeta[3])/
   (12*(e[2]^2/musq)^(2*ep)*(-1 + zb)))};


(* ::Subsubsection::Closed:: *)
(*Catani operators*)


getIcat = {
I1QCD[ns][e1_, e2_] :> 
  (qas*Cos[ep*Pi]*(-2*Cf)*(1/ep^2 + 3/2/ep))/(4*e1*e2/mu2)^ep,
I1EWK[ns][e1_,e2_] :>
	qae * (1/ep^2+3/2/ep) * 2 * (q1 q2 Cos[ep Pi] (mu2/4/e1/e2)^ep+
	q1 q3 (mu2/4/e1/e[3]/eta13)^ep + q1 q4 (mu2/4/e1/e[4]/eta14)^ep + 
	q2 q3 (mu2/4/e2/e[3]/eta23)^ep + q2 q4 (mu2/4/e2/e[4]/eta24)^ep + 
	Cos[ep Pi] q3 q4 (mu2/4/e[3]/e[4]/eta34)^ep),
(* -- *)
I1EWK[gqb] :>
	qae * (1/ep^2+3/2/ep) * 2 * (q2 q3 (mu2/4/e[2]/e[3]/eta23)^ep+
	q2 q4 (mu2/4/e[2]/e[4]/eta24)^ep + q2 q5 (mu2/4/e[2]/e[5]/eta25)^ep + 
	Cos[ep Pi] q3 q4 (mu2/4/e[3]/e[4]/eta34)^ep + 
	Cos[ep Pi] q3 q5 (mu2/4/e[3]/e[5]/eta35)^ep + 
	Cos[ep Pi] q4 q5 (mu2/4/e[4]/e[5]/eta45)^ep),
(* -- *)
I1QCD[aqb] :>
	qas * (1/ep^2 + 3/2/ep) * (-2*Cf)/(4*e[2]*e[5] eta25/mu2)^ep,
(* -- *)
I1EWK[aa] :> qae * (1/ep^2+3/2/ep) * 2 (q3 q4 Cos[ep Pi] *
	 (mu2/4/e[3]/e[4]/eta34)^ep)};


(* ::Subsubsection::Closed:: *)
(*Soft structures*)


getSoft = {
	aSgFull[ns] :> qas/ep^2 (four em^2/mu2)^(-ep)sij[1] 2 Cf,
	aSaFull[ns] :> qae/ep^2 (four em^2/mu2)^(-ep) * (-2) *(
	q1 q2 sij[1] + q1 q3 sij[eta13] + q1 q4 sij[eta14] + 
	q2 q3 sij[eta23] + q2 q4 sij[eta24] + q3 q4 sij[eta34]),
(* -- *)
	aSaFull[gqb] :> qae/ep^2 (four em^2/mu2)^(-ep) * (-2) *(
	q2 q3 sij[eta23] + q2 q4 sij[eta24] + q2 q5 sij[eta25] + 
	q3 q4 sij[eta34] + q3 q5 sij[eta35] + q4 q5 sij[eta45]),
(* -- *)
	aSgFull[aqb] :> qas/ep^2 (four em^2/mu2)^(-ep) * 2 Cf * sij[eta25],
(* -- *)
	aSaFull[aa] :> qae/ep^2 (four em^2/mu2)^(-ep) * (-2) *(
	q3 q4 sij[eta34]),
(* -- *)
	aSgFull[aqb] :> qas/ep^2 (four em^2/mu2)^(-ep) 2 Cf sij[eta25]};


(* ::Subsection:: *)
(*Join all the substitutions*)


subs=Join[getSoft,getIcat,getC,getddt,getP,getGa,getPAP,
	getTimes,getTimesIF,getConv,getPoP,getPoPAP,getPAPoPAP,getPref,
	getMaxNS,getMax,getCRV,getPRV];


(* ::Chapter::Closed:: *)
(*NLO, unsimplified, c.o.m.*)


(* ::Section::Closed:: *)
(*NS, QCD \[Checkmark] *)


r = ONLO FLM[1,2,3,4,g5] + 
	aSgFull[ns] FLM[1,2,3,4] + 
	qas * Cf * (Ci[qq,1,z] FLM[z,2,3,4] + Ci[qq,2,z] FLM[1,z,3,4]);


v = I1QCD[ns][e[1],e[2]] FLM[1,2,3,4] + FLVfin[QCD,1,2,3,4];


cv = asontwopi Cf * PAP0[qq,z]/ep * 
	(FLM[z,2,3,4] + FLM[1,z,3,4]);


nloqcd[ns] = r+v+cv//.subs;
Series[%,{ep,0,0}]//Normal//Expand;
%/.fixDelta//Simplify[#/.four->4,{em>0,e[1]>0,e[2]>0,mu2>0}]&;
PowerExpand[%];
nloqcd[ns] = 
	Collect[%,{ONLO,qas,FLM[a___],FLVfin[a___]},Simplify[#,{e[1]>0,e[2]>0,mu2>0}]&];
nloqcd[qqb] = nloqcd[ns]/.{
	FLM[a_,b_,c_,d_]:>FLM[q,qb,l,lb][a,b,c,d],
	FLM[1,2,3,4,g5]:>FLM[q,qb,l,lb,g][1,2,3,4,5],
	FLVfin[cor_,1,2,3,4]:>FLVfin[q,qb,l,lb][cor,1,2,3,4]};


(* ::Section::Closed:: *)
(*NS, EWK \[Checkmark]*)


r = ONLO FLM[1,2,3,4,a5] + 
	aSaFull[ns] FLM[1,2,3,4] + 
	qae*(q1^2 Ci[qq,1,z]FLM[z,2,3,4] + q2^2 Ci[qq,2,z] FLM[1,z,3,4]) + 
	qae*(q3^2 Cf[qq,3]+q4^2 Cf[qq,4])FLM[1,2,3,4];


v = I1EWK[ns][e[1],e[2]] FLM[1,2,3,4] + FLVfin[EWK,1,2,3,4];


cv = aemontwopi PAP0[qq,z]/ep * 
	(q1^2 FLM[z,2,3,4] + q2^2 FLM[1,z,3,4]);


nloewk[ns] = r+v+cv//.subs;
Series[%,{ep,0,0}]//Normal//Expand;
%/.fixDelta//Simplify[#/.four->4,{em>0,e[1]>0,e[2]>0,mu2>0}]&;
PowerExpand[%];
nloewk[ns] = 
	Collect[%,{ep,ONLO,qas,FLM[a___],FLVfin[a___]},Simplify[#,{e[1]>0,e[2]>0,mu2>0}]&];
nloewk[ns] = nloewk[ns]/.{q1+q2+q3+q4->0};
nloewk[qqb] = nloewk[ns]/.{
	q1 -> -Qq, q2 -> Qq, q3 -> Ql, q4 -> -Ql,
	FLM[a_,b_,c_,d_]:>FLM[q,qb,l,lb][a,b,c,d],
	FLM[1,2,3,4,a5]:>FLM[q,qb,l,lb,a][1,2,3,4,5],
	FLVfin[cor_,1,2,3,4]:>FLVfin[q,qb,l,lb][cor,1,2,3,4]};


(* ::Section::Closed:: *)
(*g qb \[Checkmark]*)


r = ONLO FLM[g,qb,l,lb,qb][1,2,3,4,5] + 
	qas * tr * (Ci[qg,1,z] FLM[q,qb,l,lb][z,2,3,4]);


cv = asontwopi tr * PAP0[qg,z]/ep * 
	(FLM[q,qb,l,lb][z,2,3,4]);


nlo[gqb] = r+cv//.subs;
Series[%,{ep,0,0}]//Normal//Expand;
%/.fixDelta//Simplify[#/.four->4,{em>0,e[1]>0,e[2]>0,mu2>0}]&;
PowerExpand[%];
nlo[gqb] = 
	Collect[%,{ep,ONLO,qas,FLM[a___],FLVfin[a___]},Simplify[#,{e[1]>0,e[2]>0,mu2>0}]&];


(* ::Section::Closed:: *)
(*\[Gamma] qb \[Checkmark]*)


r = ONLO FLM[a,qb,l,lb,qb][1,2,3,4,5] + 
	qae xn q1^2 Ci[qg,1,z]FLM[q,qb,l,lb][z,2,3,4] + 
	qae q1^2 Ci[gq,2,z] FLM[a,a,l,lb][1,zb,3,4];


cv = aemontwopi * xn q1^2 PAP0[qg,z]/ep FLM[q,qb,l,lb][z,2,3,4] + 
	 aemontwopi * q1^2 PAP0[gq,z]/ep FLM[a,a,l,lb][1,zb,3,4];


nlo[aqb] = r+cv//.subs;
Series[%,{ep,0,0}]//Normal;
nlo[aqb] = 
	Collect[%,{ep,ONLO,qas,FLM[a___],FLVfin[a___]},Simplify[#,{e[1]>0,e[2]>0,mu2>0}]&]/.{q1^2->Qq^2};


(* ::Section::Closed:: *)
(*\[Gamma] \[Gamma] \[Checkmark]*)


r = ONLO FLM[a,a,l,lb,a][1,2,3,4,5] + 
	aSaFull[aa] FLM[a,a,l,lb][1,2,3,4] + 
	qae * (q3^2 Cf[qq,3] + q4^2 Cf[qq,4]) FLM[a,a,l,lb][1,2,3,4];


v = I1EWK[aa] FLM[a,a,l,lb][1,2,3,4] + FLVfin[a,a,l,lb][EWK,1,2,3,4];


nlo[aa] = r+v//.subs/.q3->-q4;
Series[%,{ep,0,0}]/.four->4//Normal//Expand//FullSimplify[#,{em>0,mu2>0,e[3]>0,e[4]>0,eta34>0}]&;
%/.{e[3]->Sqrt[s]/2,e[4]->Sqrt[s]/2,eta34->1};
Collect[%,{ONLO,FLM[a___][b___],qae,q1,q2,q3,q4},Expand];
nlo[aa] = %/.{q4->Ql};


(* ::Chapter::Closed:: *)
(*NLO\[CircleTimes]PAP*)


(* ::Section::Closed:: *)
(*P0qq\[CircleTimes]NS[QCD] + NS[QCD]\[CircleTimes]P0qq*)


conv[PAP0[qq],NLOQCD[ns]] = q1^2 * PAP0[qq,z](ONLO FLM[q,qb][z,2,3,4,g5] + 
	aSgFull[ns] FLM[q,qb][z,2,3,4] + qas Cf Ci[qq,2,zb] FLM[q,qb][z,zb,3,4]) + 
	q1^2 qas Cf conv[Ci[qq],PAP0[qq]][1,z] FLM[q,qb][z,2,3,4] +
	(* -- *)
	q1^2 PAP0R[qq,z] (I1QCD[ns][z e[1],e[2]] FLM[q,qb][z,2,3,4] + 
		FLVfin[q,qb][QCD,z,2,3,4]) +
	q1^2 PAP0D[qq](I1QCD[ns][e[1],e[2]]FLM[q,qb][1,2,3,4] + 
		FLVfin[q,qb][QCD,1,2,3,4]) +
	(* -- *)
	q1^2 asontwopi Cf conv[cAP0[qq],cAP0[qq]][z]/ep FLM[q,qb][z,2,3,4] +
	q1^2 asontwopi Cf PAP0[qq,z] FLM[q,qb][z,zb,3,4] PAP0[qq,zb]/ep;


conv[NLOQCD[ns],PAP0[qq]] = q2^2 * PAP0[qq,zb](ONLO FLM[q,qb][1,zb,3,4,g5] + 
	aSgFull[ns] FLM[q,qb][1,zb,3,4] + qas Cf Ci[qq,1,z] FLM[q,qb][z,zb,3,4]) + 
	q2^2 qas Cf conv[Ci[qq],PAP0[qq]][2,zb] FLM[q,qb][1,zb,3,4] +
	(* -- *)
	q2^2 PAP0R[qq,zb] (I1QCD[ns][e[1], zb e[2]] FLM[q,qb][1,zb,3,4] + 
		FLVfin[q,qb][QCD,1,zb,3,4]) +
	q2^2 PAP0D[qq](I1QCD[ns][e[1],e[2]]FLM[q,qb][1,2,3,4] + 
		FLVfin[q,qb][QCD,1,2,3,4]) +
	(* -- *)
	q2^2 asontwopi Cf conv[cAP0[qq],cAP0[qq]][zb]/ep FLM[q,qb][1,zb,3,4] +
	q2^2 asontwopi Cf PAP0[qq,z]/ep FLM[q,qb][z,zb,3,4] PAP0[qq,zb];


(* ::Section::Closed:: *)
(*P0qq\[CircleTimes]NS[EWK] + NS[EWK]\[CircleTimes]P0qq*)


conv[PAP0[qq],NLOEWK[ns]] = Cf PAP0[qq,z](ONLO FLM[q,qb][z,2,3,4,a5] +
	aSaFull[ns] FLM[q,qb][z,2,3,4] + qae (q2^2 Ci[qq,2,zb] FLM[q,qb][z,zb,3,4] + 
		(q3^2 Cf[qq,3] + q4^2 Cf[qq,4]) FLM[q,qb][z,2,3,4])) +
	Cf qae q1^2 conv[Ci[qq],PAP0[qq]][1,z] FLM[q,qb][z,2,3,4] + 
	(* -- *)
	Cf PAP0R[qq,z] (I1EWK[ns][z e[1],e[2]] FLM[q,qb][z,2,3,4] + FLVfin[q,qb][EWK,z,2,3,4]) +
	Cf PAP0D[qq](I1EWK[ns][e[1],e[2]]FLM[q,qb][1,2,3,4] + FLVfin[q,qb][EWK,1,2,3,4]) +
	(* -- *)
	Cf aemontwopi q1^2 conv[cAP0[qq],cAP0[qq]][z]/ep FLM[q,qb][z,2,3,4] + 
	Cf aemontwopi q2^2 PAP0[qq,z] FLM[q,qb][z,zb,3,4] PAP0[qq,zb]/ep;


conv[NLOEWK[ns],PAP0[qq]] = Cf PAP0[qq,zb](ONLO FLM[q,qb][1,zb,3,4,a5] +
	aSaFull[ns] FLM[q,qb][1,zb,3,4] + qae (q1^2 Ci[qq,1,z] FLM[q,qb][z,zb,3,4] + 
		(q3^2 Cf[qq,3] + q4^2 Cf[qq,4]) FLM[q,qb][1,zb,3,4])) +
	Cf qae q2^2 conv[Ci[qq],PAP0[qq]][2,zb] FLM[q,qb][1,zb,3,4] + 
	(* -- *)
	Cf PAP0R[qq,zb] (I1EWK[ns][e[1],zb e[2]] FLM[q,qb][1,zb,3,4] + FLVfin[q,qb][EWK,1,zb,3,4]) +
	Cf PAP0D[qq](I1EWK[ns][e[1],e[2]]FLM[q,qb][1,2,3,4] + FLVfin[q,qb][EWK,1,2,3,4]) +
	(* -- *)
	Cf aemontwopi q2^2 conv[cAP0[qq],cAP0[qq]][zb]/ep FLM[q,qb][1,zb,3,4] + 
	Cf aemontwopi q1^2 PAP0[qq,z]/ep FLM[q,qb][z,zb,3,4] PAP0[qq,zb];


(* ::Section::Closed:: *)
(*Pqg\[CircleTimes]NS[EWK]*)


conv[PAP0[qg],NLOEWK[ns]] = tr PAP0[qg,z](ONLO FLM[q,qb][z,2,3,4,a5] +
	aSaFull[ns] FLM[q,qb][z,2,3,4] + qae (q2^2 Ci[qq,2,zb] FLM[q,qb][z,zb,3,4] + 
		(q3^2 Cf[qq,3] + q4^2 Cf[qq,4]) FLM[q,qb][z,2,3,4])) +
	tr qae q1^2 conv[Ci[qq],PAP0[qg]][1,z] FLM[q,qb][z,2,3,4] + 
	(* -- *)
	tr PAP0[qg,z] (I1EWK[ns][z e[1],e[2]] FLM[q,qb][z,2,3,4] + FLVfin[q,qb][EWK,z,2,3,4]) +
	(* -- *)
	tr aemontwopi q1^2 conv[cAP0[qq],cAP0[qg]][z]/ep FLM[q,qb][z,2,3,4] + 
	tr aemontwopi q2^2 PAP0[qg,z] FLM[q,qb][z,zb,3,4] PAP0[qq,zb]/ep;


(* ::Section::Closed:: *)
(*Pqg\[CircleTimes]NS[QCD]*)


conv[PAP0[qg],NLOQCD[ns]] = xn q1^2 * PAP0[qg,z](ONLO FLM[q,qb,l,lb,g][z,2,3,4,5] + 
	aSgFull[ns] FLM[q,qb,l,lb][z,2,3,4] + qas Cf Ci[qq,2,zb] FLM[q,qb,l,lb][z,zb,3,4]) + 
	xn q1^2 qas Cf conv[Ci[qq],PAP0[qg]][1,z] FLM[q,qb,l,lb][z,2,3,4] +
	(* -- *)
	xn q1^2 PAP0[qg,z] (I1QCD[ns][z e[1],e[2]] FLM[q,qb,l,lb][z,2,3,4] + 
		FLVfin[q,qb,l,lb][QCD,z,2,3,4])+
	(* -- *)
	xn q1^2 asontwopi Cf conv[cAP0[qq],cAP0[qg]][z]/ep FLM[q,qb,l,lb][z,2,3,4] +
	xn q1^2 asontwopi Cf PAP0[qg,z] FLM[q,qb,l,lb][z,zb,3,4] PAP0[qq,zb]/ep;


(* ::Section::Closed:: *)
(*gqb\[CircleTimes]Pqq*)


conv[NLO[gqb],PAP0[qq]] = q2^2 PAP0[qq,zb] ONLO FLM[g,qb][1,zb,3,4,5] +
	qas tr q2^2 Ci[qg,1,z] PAP0[qq,zb] FLM[q,qb][z,zb,3,4] + 
	(* -- *)
	asontwopi tr q2^2 PAP0[qg,z]/ep FLM[q,qb][z,zb,3,4] PAP0[qq,zb];


(* ::Section::Closed:: *)
(*aqb\[CircleTimes]Pqq*)


conv[NLO[aqb],PAP0[qq]] = Cf PAP0[qq,zb] ONLO FLM[a,qb,l,lb,qb][1,zb,3,4,5] +
	qae Cf xn q1^2 Ci[qg,1,z] PAP0[qq,zb] FLM[q,qb,l,lb][z,zb,3,4] + 
	qae Cf q2^2 conv[Ci[gq],PAP0[qq]][2,zb]FLM[a,a,l,lb][1,zb,3,4] + 
	(* -- *)
	aemontwopi Cf xn q1^2 PAP0[qg,z]/ep FLM[q,qb,l,lb][z,zb,3,4] PAP0[qq,zb] + 
	aemontwopi Cf q2^2 conv[cAP0[gq],cAP0[qq]][zb]/ep FLM[a,a,l,lb][1,zb,3,4];


(* ::Section::Closed:: *)
(*Pqa\[CircleTimes][qb g + qg]*)


conv[PAP0[qg],NLO[qbg]] = xn q1^2 PAP0[qg,z] ONLO FLM[qb,g,l,lb,qb][z,2,3,4,5] + 
	xn q1^2 PAP0[qg,z] qas tr Ci[qg,2,zb] FLM[qb,q,l,lb][z,zb,3,4] + 
	xn q1^2 PAP0[qg,z] asontwopi tr PAP0[qg,zb]/ep FLM[qb,q,l,lb][z,zb,3,4];
	
conv[PAP0[qg],NLO[qg]] = xn q1^2 PAP0[qg,z] ONLO FLM[q,g,l,lb,q][z,2,3,4,5] + 
	xn q1^2 PAP0[qg,z] qas tr Ci[qg,2,zb] FLM[q,qb,l,lb][z,zb,3,4] + 
	xn q1^2 PAP0[qg,z] asontwopi tr PAP0[qg,zb]/ep FLM[q,qb,l,lb][z,zb,3,4];


(* ::Section::Closed:: *)
(*[a q + a qb]\[CircleTimes]Pqg*)


conv[NLO[aq],PAP0[qg]] = tr PAP0[qg,zb] ONLO FLM[a,q,l,lb,q][1,zb,3,4,5] + 
	tr PAP0[qg,zb] qae xn q1^2 Ci[qg,1,z] FLM[qb,q,l,lb][z,zb,3,4] + 
	tr qae q1^2 conv[Ci[gq],PAP0[qg]][2,zb] FLM[a,a,l,lb][1,zb,3,4] + 
	tr aemontwopi xn q1^2 PAP0[qg,z]/ep PAP0[qg,zb] FLM[qb,q,l,lb][z,zb,3,4] + 
	tr aemontwopi q1^2 conv[cAP0[gq],cAP0[qg]][zb]/ep FLM[a,a,l,lb][1,zb,3,4];
	
conv[NLO[aqb],PAP0[qg]] = tr PAP0[qg,zb] ONLO FLM[a,qb,l,lb,qb][1,zb,3,4,5] + 
	tr PAP0[qg,zb] qae xn q1^2 Ci[qg,1,z] FLM[q,qb,l,lb][z,zb,3,4] + 
	tr qae q1^2 conv[Ci[gq],PAP0[qg]][2,zb] FLM[a,a,l,lb][1,zb,3,4] + 
	tr aemontwopi xn q1^2 PAP0[qg,z]/ep PAP0[qg,zb] FLM[q,qb,l,lb][z,zb,3,4] + 
	tr aemontwopi q1^2 conv[cAP0[gq],cAP0[qg]][zb]/ep FLM[a,a,l,lb][1,zb,3,4];


(* ::Chapter:: *)
(*NNLO, unsimplified, c.o.m.*)


(* ::Section::Closed:: *)
(*NS, ga final state \[Checkmark]*)


(* ::Subsection::Closed:: *)
(*RR*)


ds = aSgFull[ns] aSaFull[ns] FLM[1,2,3,4];


ss = aSgFull[ns] ONLO FLM[1,2,3,4,a5] + aSaFull[ns] ONLO FLM[1,2,3,4,g5] +
aSaFull[ns] * qas * Cf * (Ci[qq,1,z] FLM[z,2,3,4] + Ci[qq,2,zb] FLM[1,zb,3,4]) +
aSgFull[ns] * qae * (q1^2 Ci[qq,1,z]FLM[z,2,3,4] + q2^2 Ci[qq,2,zb]FLM[1,zb,3,4]+
	q3^2 Cf[qq,3]FLM[1,2,3,4] + q4^2 Cf[qq,4]FLM[1,2,3,4]);


sc = qas Cf CiDi[qq,1,z] * ONLO * (1-ww[1,1](1-eta15^(-ep))) FLM[z,2,3,4,a5] +
	 (* C1D1g Ci\[Gamma] *)
	 qas Cf CiDi[qq,1,z] * qae * (q2^2 Ci[qq,2,zb] FLM[z,zb,3,4] + 
					              q3^2 Cf[qq,3] FLM[z,2,3,4] + 
					              q4^2 Cf[qq,4] FLM[z,2,3,4]) +
	 qas qae Cf q1^2 Nb/2 times[qq,qq,1,z] FLM[z,2,3,4] + 
	 (* -- *)
	 qae q1^2 CiDi[qq,1,z] * ONLO * (1-ww[1,1](1-eta15^(-ep))) FLM[z,2,3,4,g5] + 
	 (* C1D1\[Gamma] Cig *)
	 qae q1^2 CiDi[qq,1,z] * qas * Cf * Ci[qq,2,zb]FLM[z,zb,3,4] +
	 qae qas q1^2 Cf Nb/2 times[qq,qq,1,z] FLM[z,2,3,4] + 
	 (* -- *)
     qas Cf CiDi[qq,2,zb] * ONLO * (1-ww[2,2](1-eta25^(-ep))) FLM[1,zb,3,4,a5] +
     (* C2D2g Ci\[Gamma] *)
	 qas Cf CiDi[qq,2,zb] * qae * (q1^2 Ci[qq,1,z] FLM[z,zb,3,4] + 
					               q3^2 Cf[qq,3] FLM[1,zb,3,4] + 
					               q4^2 Cf[qq,4] FLM[1,zb,3,4]) +
	 qas qae Cf q2^2 Nb/2 times[qq,qq,2,zb] FLM[1,zb,3,4] + 
	 (* -- *)
	 qae q2^2 CiDi[qq,2,zb] * ONLO * (1-ww[2,2](1-eta25^(-ep))) FLM[1,zb,3,4,g5] + 
	 (* C2D2\[Gamma] Cig *)
	 qae q2^2 CiDi[qq,2,zb] * qas * Cf * Ci[qq,1,z]FLM[z,zb,3,4] +
	 qae qas q2^2 Cf Nb/2 times[qq,qq,2,zb] FLM[1,zb,3,4] +
	 (* C3D3\[Gamma] Cig and same for 4 *)
	 qae * (q3^2 CfDf[qq,3] + q4^2 CfDf[qq,4]) * qas * Cf * (
	  Ci[qq,1,z] FLM[z,2,3,4] + Ci[qq,2,zb]FLM[1,zb,3,4])+
	 (* the remaining ONLO *)
	 qae (q3^2 CfDf[qq,3] + q4^2 CfDf[qq,4]) * ONLO FLM[1,2,3,4,g5];


dc = -qas * Cf * CiDi[qq,1,z] * qae * (q2^2 CiDi[qq,2,zb] FLM[z,zb,3,4] + 
                                 q3^2 CfDf[qq,3] FLM[z,2,3,4] + 
                                 q4^2 CfDf[qq,4] FLM[z,2,3,4])+
     -qas * Cf *  CiDi[qq,2,zb] * qae * (q1^2 CiDi[qq,1,z] FLM[z,zb,3,4] + 
                                  q3^2 CfDf[qq,3] FLM[1,zb,3,4] + 
                                  q4^2 CfDf[qq,4] FLM[1,zb,3,4]);


tc = maxresNS;


onnlo = ONNLO FLM[1,2,3,4,g5,a6];


rr = ds + ss + sc + dc + tc + onnlo;


(* ::Subsection::Closed:: *)
(*RV*)


(* convenient structures *)
myqcd = I1QCD[ns][e[1],e[2]] FLM[1,2,3,4] + FLVfin[QCD,1,2,3,4];
myqcdz = I1QCD[ns][z e[1],e[2]] FLM[z,2,3,4] + FLVfin[QCD,z,2,3,4];
myqcdzb = I1QCD[ns][e[1],zb e[2]] FLM[1,zb,3,4] + FLVfin[QCD,1,zb,3,4];

myewk = I1EWK[ns][e[1],e[2]] FLM[1,2,3,4] + FLVfin[EWK,1,2,3,4];
myewkz = I1EWK[ns][z e[1],e[2]] FLM[z,2,3,4] + FLVfin[EWK,z,2,3,4];
myewkzb = I1EWK[ns][e[1],zb e[2]] FLM[1,zb,3,4] + FLVfin[EWK,1,zb,3,4];


rv = ONLO (I1QCD[ns][e[1],e[2]] FLM[1,2,3,4,a5] + FLVfin[QCD,1,2,3,4,a5]) + 
(* -- *)
aSaFull[ns] myqcd + qae * (
	q1^2 * (CiR[qq,1,z] myqcdz + CiD[qq,1] myqcd) + 
	q2^2 * (CiR[qq,2,zb] myqcdzb + CiD[qq,2] myqcd) +
	q3^2 * Cf[qq,3] myqcd + q4^2 Cf[qq,4] myqcd) +
Nb/2 qae qas * Cf * (
	q1^2 CiRVDi[qq,1,z] FLM[z,2,3,4] + 
	q2^2 CiRVDi[qq,2,zb] FLM[1,zb,3,4]) +
(* -- *)
(* -- *)
ONLO (I1EWK[ns][e[1],e[2]] FLM[1,2,3,4,g5] + FLVfin[EWK,1,2,3,4,g5]) + 
(* -- *)
aSgFull[ns] myewk + qas * Cf * (
	CiR[qq,1,z] myewkz + CiD[qq,1] myewk + 
	CiR[qq,2,zb] myewkzb + CiD[qq,2] myewk ) +
Nb/2 qas qae * Cf * (
	q1^2 CiRVDi[qq,1,z] FLM[z,2,3,4] + 
	q2^2 CiRVDi[qq,2,zb] FLM[1,zb,3,4]);


(* ::Subsection::Closed:: *)
(*VV*)


HwQCDEWK = (Pi^2/2 - 6 Zeta[3] - 3/8) Cf(q1^2+q2^2);


vv = I1QCD[ns][e[1],e[2]] I1EWK[ns][e[1],e[2]] FLM[1,2,3,4] + 
	qas qae (Gamma[1-ep]/Exp[ep EulerGamma])^2 * (Exp[ep EulerGamma]/Gamma[1-ep]) * 
	HwQCDEWK/ep FLM[1,2,3,4] + 
	I1QCD[ns][e[1],e[2]] FLVfin[EWK,1,2,3,4] + I1EWK[ns][e[1],e[2]]FLVfin[QCD,1,2,3,4] + 
	FLVVpV2fin[1,2,3,4];


(* ::Subsection::Closed:: *)
(*CV*)


P0xNLOQCD = conv[PAP0[qq],NLOQCD[ns]]+conv[NLOQCD[ns],PAP0[qq]]/.
	{FLM[q,qb][a__]:>FLM[a],FLVfin[q,qb][a__]:>FLVfin[a]};


P0xNLOEWK = conv[PAP0[qq],NLOEWK[ns]]+conv[NLOEWK[ns],PAP0[qq]]/.
	{FLM[q,qb][a__]:>FLM[a],FLVfin[q,qb][a__]:>FLVfin[a]};


cv = aemontwopi/ep P0xNLOQCD + asontwopi/ep P0xNLOEWK + 
(* -- *)
-1/2 asontwopi aemontwopi/ep^2 conv[cAP0[qq],cAP0[qq]][z] * (Cf q1^2 + q1^2 Cf) FLM[z,2,3,4] + 
-1/2 asontwopi aemontwopi/ep^2 conv[cAP0[qq],cAP0[qq]][zb] * (Cf q2^2 + q2^2 Cf) FLM[1,zb,3,4] + 
(* -- *)
+1/2 asontwopi aemontwopi/ep * (q1^2 Cf PAP1[nsMIX,z]FLM[z,2,3,4] + 
	q2^2 Cf PAP1[nsMIX,zb]FLM[1,zb,3,4])+
(* -- *)
-asontwopi aemontwopi/ep^2 PAP0[qq,z] PAP0[qq,zb] * (q1^2 Cf + Cf q2^2) FLM[z,zb,3,4];


(* ::Subsection::Closed:: *)
(*Finite part check*)


rr+rv+vv+cv;
Collect[%,{_FLM,ONLO,_FLVfin},Simplify];
%//.subs;
Collect[%,{_FLM,ONLO,_FLVfin},Series[#,{ep,0,0}]&]//Normal;
Expand[%]//.fixDelta;
tmp=%;
Collect[%/.four->4,{_FLM,ONLO,_FLVfin},#//.Log[a_ b_]:>Log[a]+Log[b]//.Log[a_^b_]:>b Log[a]&];
Collect[%,{ep,qas,qae,_FLM,ONLO,_FLVfin},Expand]//.fixDelta;
res=Collect[%/.fixLog2,{ep,qas,qae,_FLM,ONLO,_FLVfin},Together];
nnlo[ns,ga]=res/.tag[int]->0; (* remove the identical quark interferences *)
nnlo[qqb,ga] = nnlo[ns,ga]/.{
	q1 -> -Qq, q2 -> Qq, q3 -> Ql, q4 -> -Ql,
	FLM[a_,b_,c_,d_]:>FLM[q,qb,l,lb][a,b,c,d],
	FLM[aa_,b_,c_,d_,a5]:>FLM[q,qb,l,lb,a][aa,b,c,d,5],
	FLM[aa_,b_,c_,d_,g5]:>FLM[q,qb,l,lb,g][aa,b,c,d,5],
	FLM[aa_,b_,c_,d_,g5,a6]:>FLM[q,qb,l,lb,g,a][aa,b,c,d,5,6],
	FLVfin[ch_,a_,b_,c_,d_]:>FLVfin[q,qb,l,lb][ch,a,b,c,d],
	FLVfin[ch_,aa_,b_,c_,d_,a5]:>FLVfin[q,qb,l,lb,a][ch,aa,b,c,d,5],
	FLVfin[ch_,aa_,b_,c_,d_,g5]:>FLVfin[q,qb,l,lb,g][ch,aa,b,c,d,5],
	FLVVpV2fin[a_,b_,c_,d_]:>FLVVpV2fin[q,qb,l,lb][a,b,c,d]};


Series[nnlo[ns,ga],{ep,0,-1}]/.fixLog2;
%/.q4->-q1-q2-q3//Simplify


(* ::Subsection:: *)
(*Simplified results*)


Geq[e1_,e2_,e3_,e4_,mu2_,eta13_,eta14_,eta23_,eta24_]:=
	3 Log[eta13 eta24/eta14/eta23] + 2 PolyLog[2,1-eta13] - 2 PolyLog[2,1-eta14]+
	-2PolyLog[2,1-eta23]+2PolyLog[2,1-eta24]+
	Log[e3^2/mu2]Log[e2 eta23/e1/eta13] + Log[e4^2/mu2]Log[e1 eta14/e2/eta24];
Ge2[e3_,e4_,ec_,eta34_]:=13-2/3Pi^2+Log[e3/e4]^2+(3-2Log[e3 e4/ec^2])Log[eta34]+
	2PolyLog[2,1-eta34];


PqqTNLO[z_,e_,mu2_]:=4DD[1,z]-2(1+z)Log[1-z]+(1-z)+Log[4e^2/mu2](2DD[0,z]-(1+z));
PqqAPR[z_]:=2DD[0,z]-(1+z);


(* assume typo in Chiara's notes, tag \[Rule] 2 *)
PqqNNLO[z_,e_]:=DD[0,z](12 Log[4 e^2/mu2]^2 - 16/3 Log[z]^3) + 
	Zeta[3](32 DD[0,z] - 4(7z^2+1)/(z-1))+
	DD[1,z](8 Log[4 e^2/mu2](Log[4 e^2/mu2]+3)-8 Log[z]^2) + 
	DD[2,z](24 Log[4e^2/mu2]-16Log[z])+
	16 DD[3,z] + 2PolyLog[2,1-z]/(1-z)(2(1+z^2)Log[4e^2/mu2]+z^2-2z+1+3(z^2+1)Log[1-z])+
	4PolyLog[2,z](-(1+z)Log[4e^2/mu2]+z-1-(1+z^2)Log[z]/(z-1)-2(1+z)Log[1-z])+
	Log[1-z]^2(-12(1+z)Log[4e^2/mu2]+4(z-1)-(5z^2+13)Log[z]/(z-1))+
	Log[1-z](2Log[z](2(1+z^2)Log[4e^2/mu2]-5(1+z^2)-2z)/(z-1)+4(z+1)Log[z]^2+
	1/3(-12(1+z)Log[4e^2/mu2]^2-24(z+2)Log[4e^2/mu2]+4Pi^2(1+z)-27z+36))+
	-Log[z]/3/(z-1)(-3(1+3z^2)Log[4e^2/mu2]^tag+12(z^2+1+z)Log[4e^2/mu2]+
		8Pi^2(1+z^2)+27z^2-42z+15)+
	Log[z]^2((1+z)Log[4e^2/mu2]-z+2)+
	2/3(-6(2+z)Log[4e^2/mu2]^2+Pi^2(1+z)Log[4e^2/mu2]+18Log[4e^2/mu2]-15z Log[4e^2/mu2]+
		-5Pi^2z + 24z + 5 Pi^2-27)+
	-2(3z^2-5)PolyLog[3,1-z]/(z-1)+4(3z^2+5)PolyLog[3,z]/(z-1)-(13z^2+47)Log[z]^3/6/(z-1)+
	-8(z+1)Log[1-z]^3/.tag->2;


Coefficient[nnlo[qqb,ga],FLM[q,qb,l,lb][z,2,3,4]Qq^2 Cf qas qae];
ttt=%;
ttt/.e[1]->e/.e[2]:>e/.Log[e]:>1/2 Log[e^2]/.Log[e^2]:>Log[4e^2/mu2]-2Log[2]+Log[mu2]//Expand;
Collect[%/.fixLi2,{Log[4e^2/mu2],DD[a_,z_],Log[a_],PolyLog[a_,b_],Zeta[3],Pi},Simplify];
%/.(1+z^2)/(-1+z):>-"[(1+z^2)/(1-z)]";
%/.PolyLog[3,z]:>"[Li3[z]-Zeta[3]]"+Zeta[3];
Collect[%/.fixLi2,{Log[4e^2/mu2],DD[a_,z_],"[1+z^2]/(1-z)]",Log[a_],PolyLog[a_,b_],
	"[Li3[z]-Zeta[3]]",Zeta[3],Pi},Simplify]//InputForm;

myPqqNNLO[z_,e_,mu2_] := 
(-8*Pi^2*(-1 + z))/3 + 2*(-9 + 8*z) + (4*"[Li3[z]-Zeta[3]]"*(5 + 3*z^2))/(-1 + z) + 
 16*DD[3, z] - 8*(1 + z)*Log[1 - z]^3 + (2 - z)*Log[z]^2 + 
 ((15 + 13*z^2)*Log[z]^3)/(6 - 6*z) + 
 Log[1 - z]^2*(4*(-1 + z) + ((-5 + 3*z^2)*Log[z])/(-1 + z)) + 
 Log[(4*e^2)/mu2]^2*(-4*(2 + z) + 12*DD[0, z] + 8*DD[1, z] - 4*(1 + z)*Log[1 - z] + 
   ((1 + 3*z^2)*Log[z])/(-1 + z)) + (6 - 6*z)*PolyLog[2, 1 - z] + 
 Log[z]*(5 + (10*"[(1+z^2)/(1-z)]"*Pi^2)/3 - 9*z - 4*"[(1+z^2)/(1-z)]"*PolyLog[2, 1 - z]) + 
 Log[(4*e^2)/mu2]*(12 - 10*z + 24*DD[1, z] + 24*DD[2, z] - 12*(1 + z)*Log[1 - z]^2 - 
   (4*(1 + z + z^2)*Log[z])/(-1 + z) + (1 + z)*Log[z]^2 + 
   Log[1 - z]*(-8*(2 + z) + (8*z^2*Log[z])/(-1 + z)) - (8*PolyLog[2, 1 - z])/(-1 + z)) + 
 Log[1 - z]*(12 - 9*z - (2*(7 - 2*z + 7*z^2)*Log[z])/(-1 + z) - 
   8*"[(1+z^2)/(1-z)]"*Log[z]^2 + (2*(-7 + z^2)*PolyLog[2, 1 - z])/(-1 + z)) + 
 ((10 - 6*z^2)*PolyLog[3, 1 - z])/(-1 + z) - 16*(1 + z)*Zeta[3] + 32*DD[0, z]*Zeta[3];


FLMzzb = 2*Cf*qae*qas*Qq^2*(-1 + z - 4*DD[1, z] + 
  (1 + z)*Log[s/mu2] - 2*DD[0, z]*Log[s/mu2] + 
  2*(1 + z)*Log[1 - z])*(-1 + zb - 4*DD[1, zb] + 
  (1 + zb)*Log[s/mu2] - 2*DD[0, zb]*Log[s/mu2] + 
  2*(1 + zb)*Log[1 - zb]);


myPqqNNLO[z,e,mu2]-PqqNNLO[z,e]/.fixLi2/."[Li3[z]-Zeta[3]]":>PolyLog[3,z]-Zeta[3]/.
"[(1+z^2)/(1-z)]":>(1+z^2)/(1-z)//Simplify//Expand;
%/.fixDelta//Simplify


testres = ONNLO FLM[q,qb,l,lb,g,a][1,2,3,4,5,6] + FLVVpV2fin[q,qb,l,lb][1,2,3,4] + 
ONLO FLVfin[q,qb,l,lb,a][QCD,1,2,3,4,5]+ ONLO FLVfin[q,qb,l,lb,g][EWK,1,2,3,4,5] + 
(* -- FLVfinQCD below -- *)
qae FLVfin[q,qb,l,lb][QCD,1,2,3,4]*(Qq^2(2/3 Pi^2 + 3 Log[4e^2/mu2]) + 
	Ql Qq Geq[e1,e2,e3,e4,e^2,eta13,eta14,eta23,eta24]+
	Ql^2 Ge2[e3,e4,e,eta34]) + 
(* - *)
qae PqqTNLO[z,e,mu2]Qq^2 FLVfin[q,qb,l,lb][QCD,z,2,3,4] + 
qae PqqTNLO[zb,e,mu2]Qq^2 FLVfin[q,qb,l,lb][QCD,1,zb,3,4] + 
(* -- FLVfinEWK below -- *)
qas FLVfin[q,qb,l,lb][EWK,1,2,3,4]Cf(2/3 Pi^2 + 3 Log[4 e^2/mu2]) + 
(* - *)
qas PqqTNLO[z,e,mu2]Cf FLVfin[q,qb,l,lb][EWK,z,2,3,4] + 
qas PqqTNLO[zb,e,mu2]Cf FLVfin[q,qb,l,lb][EWK,1,zb,3,4] + 
(* -- FLM12 below -- *)
qae qas Cf FLM[q,qb,l,lb][1,2,3,4]*(
	Qq^2(16/45Pi^4+4(Pi^2+8Zeta[3])Log[4e^2/mu2]+(9-4/3Pi^2)Log[4e^2/mu2]^2)+
	Ql Qq Geq[e1,e2,e3,e4,e^2,eta13,eta14,eta23,eta24](2/3Pi^2 + 3 Log[4e^2/mu2])+
	Ql^2 Ge2[e3,e4,e,eta34](2/3Pi^2+3 Log[4e^2/mu2])) + 
(* - *)
qae qas 2 Cf Qq^2 PqqTNLO[z,e,mu2] FLM[q,qb,l,lb][z,zb,3,4]PqqTNLO[zb,e,mu2] + 
(* - *)
qae qas Cf Qq^2 (myPqqNNLO[z,e,mu2]FLM[q,qb,l,lb][z,2,3,4] + 
	myPqqNNLO[zb,e,mu2]FLM[q,qb,l,lb][1,zb,3,4]) +
qae qas Cf Ql Qq PqqTNLO[z,e,mu2](Geq[e1,e2,e3,e4,mu2,eta13,eta14,eta23,eta24] + 
	Log[e^2/mu2] Log[eta13 eta24/eta14/eta23] - 2 Log[e3 eta13/e4/eta14]Log[z])*
	FLM[q,qb,l,lb][z,2,3,4] + 
qae qas Cf Ql Qq PqqTNLO[zb,e,mu2](Geq[e1,e2,e3,e4,mu2,eta13,eta14,eta23,eta24] + 
	Log[e^2/mu2] Log[eta13 eta24/eta14/eta23] + 2 Log[e3 eta23/e4/eta24]Log[zb])*
	FLM[q,qb,l,lb][1,zb,3,4] + 
qae qas Cf Ql^2 Ge2[e3,e4,e,eta34] * (PqqTNLO[z,e,mu2]FLM[q,qb,l,lb][z,2,3,4] + 
	PqqTNLO[zb,e,mu2]FLM[q,qb,l,lb][1,zb,3,4]) + 
(* -- ONLO below -- *)
qae Qq^2 (PqqTNLO[z,e,mu2] + ww[1,1] Log[eta15] PqqAPR[z]) *
	 ONLO FLM[q,qb,l,lb,g][z,2,3,4,5]+
qae Qq^2 (PqqTNLO[zb,e,mu2] + ww[2,2] Log[eta25] PqqAPR[zb]) *
	 ONLO FLM[q,qb,l,lb,g][1,zb,3,4,5]+
(* - *)
qas Cf (PqqTNLO[z,e,mu2] + ww[1,1] Log[eta15] PqqAPR[z]) *
	 ONLO FLM[q,qb,l,lb,a][z,2,3,4,5]+
qas Cf (PqqTNLO[zb,e,mu2] + ww[2,2] Log[eta25] PqqAPR[zb]) *
	 ONLO FLM[q,qb,l,lb,a][1,zb,3,4,5]+
(* - *)
(Qq^2 (2/3 Pi^2 + 3 Log[4e^2/mu2])+Ql Qq Geq[e,e,e3,e4,e^2,eta13,eta14,eta23,eta24] +
	Ql^2 Ge2[e3,e4,e,eta34]) qae ONLO FLM[q,qb,l,lb,g][1,2,3,4,5]+
(* - *)
(2/3 Pi^2 + 3 Log[4e^2/mu2]) qas Cf ONLO FLM[q,qb,l,lb,a][1,2,3,4,5];


testres-nnlo[qqb,ga];
%/.{e[1]->e1,e[2]->e2,e[3]->e3,e[4]->e4};
%/.Log[4e^2/mu2]->2Log[2]+2Log[e]-Log[mu2];
diff=Collect[%,{qae,qas,Cf,Qq,Ql,FLVfin[a___][b___],FLM[a__][b__],DD[a_,b_],
	Log[a_],PolyLog[a_,b_]},Simplify];
diff=Expand[%]/.fixLog2/.fixDelta/.fixLi2/.{e1->e,e2->e,em->e};

diff//PowerExpand//Expand;
%/.{zb->z,"[Li3[z]-Zeta[3]]":>PolyLog[3,z]-Zeta[3],"[(1+z^2)/(1-z)]" ->(1+z^2)/(1-z)}//
	Expand//Apart[#,z]&//Expand//Simplify


(* ::Section::Closed:: *)
(*NS, qqb final state [factor of 2 difference w.r.t Chiara]*)


(* ::Subsection:: *)
(*RR*)


tc = qas qae (4e[1]^2/mu2)^(-2ep) 2^(4 ep) * (
	twoISRab[z,5]) Cf q1^2 FLM[z,2,3,4] + 
	qas qae (4e[2]^2/mu2)^(-2ep) 2^(4 ep) * (
	twoISRab[zb,5]) Cf q2^2 FLM[1,zb,3,4];
	
rr = tc;


(* ::Subsection:: *)
(*CV*)


cv = +1/2 asontwopi aemontwopi/ep * 
	(q1^2 Cf PAP1[nsMIX,z]FLM[z,2,3,4] + q2^2 Cf PAP1[nsMIX,zb]FLM[1,zb,3,4]);


(* ::Subsection:: *)
(*Finite part checks*)


nnlo[ns,qqb] = (rr//.subs) + Coefficient[cv//.subs,tag[int]];
nnlo[ns,qqb] = Series[nnlo[ns,qqb],{ep,0,0}];
nnlo[qqb,qqb] = nnlo[ns,qqb]/.{FLM[a_,b_,c_,d_]:>FLM[q,qb,l,lb][a,b,c,d]};
nnlo[qqb,qqb] = nnlo[qqb,qqb]/.{q1->-Qq,q2->Qq};
Series[%,{ep,0,-1}]//Simplify


(* ::Subsection:: *)
(*Simplified results*)


(* Lmu = Log[mu2/mv2] = Log[mu2/z/s] *)
PqqNNLOqqb[z_] := Cf Qq^2 * (1 + 7*z + (2*Pi^2*(1 + z))/3 + 
((5 - 12*z + 4*z^2)*Log[z]^2)/(2*(-1 + z)) + 
 Log[1 - z]*(4*(-8 + 7*z) + ((20 - 8*z^2)*Log[z])/(-1 + z)) + 
 Log[z]*((6 + 11*z - 27*z^2)/(-1 + z) + 8*(1 + z)*Log[1 + z]) - 
 (2*(-13 + 6*z + z^2)*PolyLog[2, 1 - z])/(-1 + z) + 
 Lmu*(16 - 14*z + (2*(-5 + 2*z^2)*Log[z])/(-1 + z) + 
   pqq*(2*Log[z]^2 + 4*PolyLog[2, 1 - z])) + 8*(1 + z)*PolyLog[2, -z] + 
 pqq*(-Log[z]^3/3 + Log[1 - z]*(5*Log[z]^2 - 8*PolyLog[2, 1 - z]) + 
   Log[z]*((-7*Pi^2)/3 + 10*PolyLog[2, 1 - z] - 8*PolyLog[2, -z]) + 
   12*PolyLog[3, 1 - z] + 16*PolyLog[3, -z] + 18*PolyLog[3, z] - 6*Zeta[3]));


testresqqb = qae qas * (FLM[q,qb,l,lb][z,2,3,4] PqqNNLOqqb[z] +  
FLM[q,qb,l,lb][1,zb,3,4] PqqNNLOqqb[zb]);


nnlo[qqb,qqb]- testresqqb/.{pqq->(1+z^2)/(1-z)}/.zb->z//Normal;
%/.fixLi2;
%/.Lmu->-Log[s/mu2]-Log[z]/.Log[s/mu2]->Log[e^2/mu2]+2Log[2];
%/.{e[1]:>e,e[2]:>e}//Simplify


(* ::Section::Closed:: *)
(*qq -> qq [factor of 4 difference w.r.t. Chiara]*)


(* ::Subsection:: *)
(*RR*)


tc = qas qae (4e[1]^2/mu2)^(-2ep) 2^(4 ep) * (
	twoISRab[z,6]) Cf q1^2 FLM[qb,q][z,2,3,4] + 
	qas qae (4e[2]^2/mu2)^(-2ep) 2^(4 ep) * (
	twoISRab[zb,6]) Cf q2^2 FLM[q,qb][1,zb,3,4];
	
rr = tc;


(* ::Subsection::Closed:: *)
(*CV*)


cv = +1/2 asontwopi aemontwopi/ep * 
	(q1^2 Cf PAP1[qqbVMix,z]FLM[qb,q][z,2,3,4] + q2^2 Cf PAP1[qqbVMix,zb]FLM[q,qb][1,zb,3,4]);


(* ::Subsection:: *)
(*Finite part checks*)


nnlo[ns,qq] = rr + cv//.subs/.{q1^2->Qq^2,q2^2->Qq^2};
nnlo[ns,qq] = Series[nnlo[ns,qq],{ep,0,0}];
nnlo[qq,qq] = nnlo[ns,qq]/.
	{FLM[q,qb][a_,b_,c_,d_]:>FLM[q,qb,l,lb][a,b,c,d], FLM[qb,q][a___]:>FLM[qb,q,l,lb][a]};
Series[%,{ep,0,-1}]//Simplify;
Simplify[%]


(* ::Subsection:: *)
(*Simplified result*)


PqqNNLOqq[z_] = Cf Qq^2 (Pi^2*(-1 - z) - 15*(-1 + z) + 4*Log[z]^2 + 
 Log[1 - z]*(-16*(-1 + z) + 8*(1 + z)*Log[z]) + 
 Log[z]*(11 + 19*z - 12*(1 + z)*Log[1 + z]) + 
 4*(3 + z)*PolyLog[2, 1 - z] - 12*(1 + z)*PolyLog[2, -z] + 
 Lmu*(8*(-1 + z) - 4*(1 + z)*Log[z] + 
   pqqmx*((2*Pi^2)/3 - 2*Log[z]^2 + 8*Log[z]*Log[1 + z] + 
     8*PolyLog[2, -z])) + pqqmx*(Log[z]^3/3 - 2*Pi^2*Log[1 + z] + 
   6*Log[z]^2*Log[1 + z] + 4*Log[1 + z]^3 + 
   Log[1 - z]*((-4*Pi^2)/3 - 4*Log[z]^2 - 16*Log[z]*Log[1 + z] - 
     16*PolyLog[2, -z]) + Log[z]*((8*Pi^2)/3 - 12*Log[1 + z]^2 - 
     4*PolyLog[2, 1 - z] + 12*PolyLog[2, -z]) + 16*PolyLog[3, 1 - z] - 
   36*PolyLog[3, -z] - 16*PolyLog[3, z] - 24*PolyLog[3, z/(1 + z)] - 
   8*PolyLog[3, 1 - z^2] + 10*Zeta[3]));


testresqq = qae qas * (FLM[qb,q,l,lb][z,2,3,4] PqqNNLOqq[z] +  
FLM[q,qb,l,lb][1,zb,3,4] PqqNNLOqq[zb]);


nnlo[qq,qq]- testresqq/.{pqqmx->(1+z^2)/(1+z)}/.zb->z//Normal;
%/.fixLi2;
%/.Lmu->-Log[s/mu2]-Log[z]/.Log[s/mu2]->Log[e^2/mu2]+2Log[2];
%/.{e[1]:>e,e[2]:>e}//Simplify


(* ::Section::Closed:: *)
(*g qb \[Checkmark]*)


(* g[1] qb[2] -> l[3] lb[4] qb[5] \[Gamma][6] *)


(* ::Subsection::Closed:: *)
(*RR*)


ss = aSaFull[gqb] ONLO FLM[g,qb][1,2,3,4,5] +
	 aSaFull[ns] * qas * tr * (Ci[qg,1,z] FLM[q,qb][z,2,3,4]);


(* single collinear: 1 = w65 + w62 + w63 + w64 -> Ct51[thC w65 + w62 + w63 + w64] + 
	Ct65[thB + thD]w65 + Ct62 w62 + Ct63 w63 + Ct64 w64,
   all operators add on the phase space. Then use w62+w63+w64 = 1 - w65 *)

(* here Nc is a prefactor that comes in sectors b and d *)

(* the connection between new and old notation is that gatqq = Nc ga[qq,R,2,2] *)

(* note that for the regulation, I do not act on the phase space. 
   But then the b/d term (eta51/(1-eta51))^(-ep) -> x^-ep 
   and x^(-ep) * [x(1-x)]^-ep/x -> Nb *)
 
sc = qas tr CiDi[qg,1,z] * ONLO FLM[q,qb][z,2,3,4,a5] * (1-ww[1,1](1-(eta15/2)^(-ep))) + 
	 qae q5^2 CfDf[qq,5] * ONLO FLM[g,qb][1,2,3,4,5] * Nc *(eta15/4/(1-eta15))^(-ep) + 
	 qae q2^2 CiDi[qq,2,zb] * ONLO FLM[g,qb][1,zb,3,4,5] + 
	 qae * (q3^2 CfDf[qq,3] + q4^2 CfDf[qq,4]) * ONLO FLM[g,qb][1,2,3,4,5] + 
(* -- regulation below, in the same order of the comment in the first line *)
	 qas qae tr q1^2 Nb/2 (1/2)^(-ep) times[qq,qg,1,z] FLM[q,qb][z,2,3,4] + 
	 qas qae tr q2^2 CiDi[qg,1,z] Ci[qq,2,zb] FLM[q,qb][z,zb,3,4] + 
	 qas qae tr CiDi[qg,1,z] * (q3^2 Cf[qq,3] + q4^2 Cf[qq,4]) FLM[q,qb][z,2,3,4] + 
	 (* -- *)
	 qae qas tr q1^2 (4^ep Nb/2) Nc timesIF[qgI,qqF,1,z] FLM[q,qb][z,2,3,4] + 
	 (* -- *)
	 qae qas tr Ci[qg,1,z] q2^2 CiDi[qq,2,zb] FLM[q,qb][z,zb,3,4] + 
	 qae qas tr Ci[qg,1,z] * (q3^2 CfDf[qq,3] + q4^2 CfDf[qq,4]) FLM[q,qb][z,2,3,4];


dc = -qas tr CiDi[qg,1,z] * qae * (q2^2 CiDi[qq,2,zb] FLM[q,qb][z,zb,3,4] + 
								   q3^2 CfDf[qq,3] FLM[q,qb][z,2,3,4] + 
								   q4^2 CfDf[qq,4] FLM[q,qb][z,2,3,4]);


(* max does not include the (1-ep) averaging factors *)
tc = qas qae/ep tr q1^2(4e[1]^2/mu2)^(-2ep) 2^(4 ep)*((-ep) ISab[z,9,e[1]]/(1-ep)) FLM[q,qb][z,2,3,4];


onnlo = ONNLO FLM[g,qb,l,lb,qb,a][1,2,3,4,5,6];


rr = ss + sc + dc + tc + onnlo;


(* ::Subsection::Closed:: *)
(*RV*)


(* convenient structures *)
myewkz = I1EWK[ns][z e[1], e[2]] FLM[q,qb][z,2,3,4] + FLVfin[q,qb][EWK,z,2,3,4];


rv = ONLO (I1EWK[gqb] FLM[g,qb][1,2,3,4,5] + FLVfin[g,qb][EWK,1,2,3,4,5]) + 
(* -- *)
qas tr Ci[qg,1,z] myewkz +
Nb/2 qae qas * tr * q1^2 CiRVDi[qg,1,z] FLM[q,qb][z,2,3,4];


(* ::Subsection::Closed:: *)
(*CV*)


P0xNLOQCD = conv[NLO[gqb],PAP0[qq]];
P0xNLOEWK = conv[PAP0[qg],NLOEWK[ns]];


cv = aemontwopi/ep P0xNLOQCD + asontwopi/ep P0xNLOEWK + 
(* -- *)
-1/2 asontwopi aemontwopi/ep^2 conv[cAP0[qq],cAP0[qg]][z] * (tr q1^2) FLM[q,qb][z,2,3,4] + 
(* -- *)
+1/2 asontwopi aemontwopi/ep * (q1^2 tr PAP1[qgAB,z]FLM[q,qb][z,2,3,4])+
(* -- *)
- asontwopi aemontwopi/ep^2 PAP0[qg,z] PAP0[qq,zb] * (tr q2^2) FLM[q,qb][z,zb,3,4];


(* ::Subsection:: *)
(*Finite part check*)


rr+rv+cv;
Collect[%,{_FLM,ONLO,_FLVfin}];
%//.subs;
ttt=Series[%,{ep,0,0}];
Normal[%]//Expand;
ttt=%//.fixDeltaCH/.fixLog2/.fixLi2/.fixLogz;
ONLObit = Coefficient[ttt,ONLO] ONLO;
rest = ttt/.ONLO:>0;
ONLObit=Collect[ONLObit/.four->4/.{q3->Ql,q4->-Ql,q2->Qq,q5->-Qq}//PowerExpand,
	{ep,qae,qam,FLM[a__][b___]},Simplify];
rest=Collect[rest/.four->4/.{q3->Ql,q4->-Ql,q1->-Qq,q2->Qq}//PowerExpand,
	{ep,qae,qam,FLM[a__][b___],Log[a_],PolyLog[a_,b_],Zeta[a_],Pi},Simplify];
nnlo[gqb] = ONLObit + rest;
nnlo[gqb] = nnlo[gqb]/.{
	FLM[g,qb][a__]:>FLM[g,qb,l,lb,qb][a],
	FLM[q,qb][a_,b_,c_,d_]:>FLM[q,qb,l,lb][a,b,c,d],
	FLM[q,qb][aa_,bb_,cc_,dd_,a5]:>FLM[q,qb,l,lb,a][aa,bb,cc,dd,5],
	FLM[q,qb][aa_,bb_,cc_,dd_,g5]:>FLM[q,qb,l,lb,g][aa,bb,cc,dd,5],
	FLVfin[q,qb][ch_,aa_,bb_,cc_,dd_]:>FLVfin[q,qb,l,lb][ch,aa,bb,cc,dd],
	FLVfin[g,qb][ch_,aa_,bb_,cc_,dd_,5]:>FLVfin[g,qb,l,lb,qb][ch,aa,bb,cc,dd,5]};


(* ::Subsection:: *)
(*Simplified results*)


Gq2b[e5_,e_,mu2_,eta15_,eta25_] := 13/2 - Pi^2 + Log[e5/e]^2 + 
	2*(3/2 - Log[e5/e])*Log[4*eta25] + 3/2*Log[e^2/mu2] + 
	2*PolyLog[2, 1 - eta25] + 
	(2Log[e5/e]-3/2)Log[eta15/(1-eta15)];


testresgqb = ONNLO FLM[g,qb,l,lb,qb,a][1,2,3,4,5,6] + 
(* -- ONLO below -- *)
ONLO qae * (
	Ql^2 Ge2[e3,e4,e,eta34] + 
	Ql Qq (Geq[e,e,e3,e4,e^2,eta35,eta45,eta23,eta24]-2 Log[e5/e]Log[e3 eta35/e4/eta45]) +
	Qq^2 Gq2b[e5,e,mu2,eta15,eta25]) FLM[g,qb,l,lb,qb][1,2,3,4,5] + 
(* -- *)
qae Qq^2 PqqTNLO[zb,e,mu2] ONLO FLM[g,qb,l,lb,qb][1,zb,3,4,5];


(* ::Section::Closed:: *)
(*\[Gamma] qb \[Checkmark]*)


(* \[Gamma][1] qb[2] -> l[3] lb[4] g[5] qb[6] *)


(* ::Subsection::Closed:: *)
(*RR*)


(* for the aa bit: the soft gives an extra eta factor -> Gamma[1-ep]^2/Gamma[1-2ep] replaced by Nb/2 *)
(* also, in the eta ->0 limit, sij[eta]->Gamma[1+ep]Gamma[1-ep]sij[1] *)
ss = aSgFull[aqb] ONLO FLM[a,qb,l,lb,qb][1,2,3,4,5] +
	 aSgFull[ns] * qae * xn q1^2 * (Ci[qg,1,z] FLM[q,qb,l,lb][z,2,3,4])+
	 aSgFull[ns] * qae * q2^2 * Gamma[1+ep]Gamma[1-ep](Nb/2 CiDi[gq,2,zb])FLM[a,a,l,lb][1,zb,3,4];


(* single collinear: 
	1 = w5161 + w5262 + w5162 + w5261 \[Rule] 
	w5161(thB C56 + thC C61 + thD C65) + 
	w5262(thA C52 + thB C56 + thD C65) +
	w5162 + w5261(C52 + C61) \[Rule] 
	C52 (thA w62 + w61 \[Rule] C52 + (thA-1) w62 C52 +
	C61 (thC w51 + w52 \[Rule] C61 + (thC-1) w51 C61 +
	C65 [thB + thD] (w5161 + w5262)
*)

(* here Nc is a prefactor that comes in sectors b and d *)

(* the connection between new and old notation is that gatqq = Nc ga[qq,R,2,2] *)

(* note that for the regulation, I do not act on the phase space. 
   But then the b/d term (eta51/(1-eta51))^(-ep) -> x^-ep 
   and x^(-ep) * [x(1-x)]^-ep/x -> Nb *)

sc = qas Cf CiDi[qq,2,zb] * ONLO FLM[a,qb,l,lb,qb][1,zb,3,4,5] * 
		(1-w[2,2](1-(eta25/2)^(-ep))) + 
	 qae xn q1^2 CiDi[qg,1,z] * ONLO FLM[q,qb,l,lb,g][z,2,3,4,5] * 
	    (1-w[1,1](1-(eta15/2)^(-ep))) +
	 qas * Cf * CfDf[qq,5] * ONLO FLM[a,qb,l,lb,qb][1,2,3,4,5] * 
	    Nc * ((eta15/4/(1-eta15))^(-ep) w56[1,1] + 
	          (eta25/4/(1-eta25))^(-ep) w56[2,2])+
(* -- regulation below *)
	 qas*Cf*qae*(CiDi[qq,2,zb] xn q1^2 * Ci[qg,1,z] * FLM[q,qb,l,lb][z,zb,3,4]+
	  q2^2 * Nb/2 (1/2)^(-ep) times[gq,qq,2,zb] * FLM[a,a,l,lb][1,zb,3,4]) +
	 qae*xn*q1^2*qas*Cf*(CiDi[qg,1,z]Ci[qq,2,zb]FLM[q,qb,l,lb][z,zb,3,4] + 
	  Nb/2 (1/2)^(-ep) times[qq,qg,1,z] FLM[q,qb,l,lb][z,2,3,4]) + 
	 qas * qae * Cf (4^ep Nb/2) Nc * 
	  (xn q1^2 timesIF[qgI,qqF,1,z] FLM[q,qb,l,lb][z,2,3,4] + 
	      q2^2 timesIF[gqI,qqF,2,zb] FLM[a,a,l,lb][1,zb,3,4]);


dc = -qae xn q1^2 CiDi[qg,1,z] * qas * Cf * CiDi[qq,2,zb] FLM[q,qb,l,lb][z,zb,3,4];


(* max does not include the (1-ep) averaging factors *)
tc = qas qae/ep Cf*xn*q1^2(4e[1]^2/mu2)^(-2ep) 2^(4 ep)*((-ep) ISab[z,9,e[1]]/(1-ep)) * 
	FLM[q,qb,l,lb][z,2,3,4] + 
	 qas qae/ep q2^2 Cf (4e[2]^2/mu2)^(-2ep) 2^(4 ep)*((-ep) ISab[zb,8,e[2]]*(1-ep)) * 
	FLM[a,a,l,lb][1,zb,3,4];


onnlo = ONNLO FLM[a,qb,l,lb,g,qb][1,2,3,4,5,6];


rr = ss + sc + dc + tc + onnlo;


(* ::Subsection::Closed:: *)
(*RV*)


(* convenient structures *)
myqcdz = I1QCD[ns][z e[1], e[2]] FLM[q,qb,l,lb][z,2,3,4] + FLVfin[q,qb,l,lb][QCD,z,2,3,4];


rv = ONLO (I1QCD[aqb] FLM[a,qb,l,lb,qb][1,2,3,4,5] + FLVfin[a,qb,l,lb,qb][QCD,1,2,3,4,5]) + 
(* -- *)
qae xn q1^2 Ci[qg,1,z] myqcdz +
Nb/2 qae qas * xn q1^2 * Cf CiRVDi[qg,1,z] FLM[q,qb,l,lb][z,2,3,4] + 
Nb/2 qae qas * q2^2 * Cf CiRVDi[gq,2,zb] FLM[a,a,l,lb][1,zb,3,4];


(* ::Subsection::Closed:: *)
(*CV*)


P0xNLOEWK = conv[NLO[aqb],PAP0[qq]];
P0xNLOQCD = conv[PAP0[qg],NLOQCD[ns]];


cv = aemontwopi/ep P0xNLOQCD + asontwopi/ep P0xNLOEWK + 
(* -- *)
-1/2 asontwopi aemontwopi/ep^2 conv[cAP0[qq],cAP0[qg]][z] * (xn q1^2 Cf) FLM[q,qb,l,lb][z,2,3,4] + 
(* -- *)
-1/2 asontwopi aemontwopi/ep^2 conv[cAP0[gq],cAP0[qq]][zb] * (q2^2 Cf) FLM[a,a,l,lb][1,zb,3,4] + 
(* -- *)
+1/2 asontwopi aemontwopi/ep * (xn q1^2 Cf PAP1[qgAB,z]FLM[q,qb,l,lb][z,2,3,4])+
(* -- *)
+1/2 asontwopi aemontwopi/ep * (q2^2 Cf PAP1[gqAB,zb]FLM[a,a,l,lb][1,zb,3,4])+
(* -- *)
- asontwopi aemontwopi/ep^2 PAP0[qg,z] PAP0[qq,zb] * (xn q1^2 Cf) FLM[q,qb,l,lb][z,zb,3,4];


(* ::Subsection::Closed:: *)
(*Finite part check*)


mytmp = rr+rv+cv//.subs;


(* ONLO *)
Coefficient[mytmp/.ONLO:> tag ONLO,tag];
Series[%,{ep,0,-1}]//Normal;
Expand[%]//.fixDeltaCH;
ttt=Collect[%/.four->4//PowerExpand,{ep,ONLO,FLM[a___][b___],FLVfin[a___][b___]},Simplify]
%/.{w[1,1]+w[2,2]->1}


(* rest FLVfin *)
mytmp/.{ONLO :> 0, FLM[a___][b___]:> 0};
Series[%,{ep,0,-1}]//Simplify


mytmp/.{ONLO :>0, FLVfin[a___][b___]:>0};
Series[%,{ep,0,-1}]//Normal//Expand;
ttt=%//.fixDeltaCH;
Cases[Variables[ttt],FLM[a___][b___]]


(* double boosted *)
Coefficient[ttt,FLM[q,qb,l,lb][z,zb,3,4]];
%/.four->4//PowerExpand;
Series[%,{ep,0,-1}]//Simplify


(* single boosted z *)
Coefficient[ttt,FLM[q,qb,l,lb][z,2,3,4]];
%/.four->4/.fixLi2//PowerExpand;
Series[%,{ep,0,-1}]//Simplify;
%/.{t1->t2,t3->t2}//Simplify


(* single boosted zb *)
Coefficient[ttt,FLM[a,a,l,lb][1,zb,3,4]];
%/.four->4/.fixLi2//PowerExpand;
Series[%,{ep,0,-1}]//Simplify;
%/.{t1->t2,t3->t2}//Simplify


nnlo[aqb] = mytmp;
Series[%,{ep,0,0}];
Normal[%/.four->4]//Expand//PowerExpand;
%/.fixDeltaCH/.fixLi2/.fixLogz/.fixLi2;
nnlo[aqb]=Collect[%,{ep,qae,qas,ONLO,FLM[a__][b__],FLVfin[a__][b___]},Simplify]/.{w[1,1]+w[2,2]->1}/.
	{q1^2->Qq^2,q2^2->Qq^2};


(* ::Subsection:: *)
(*Simplified result*)


PqgTNLO[z_,e_,mu2_]:=2(1-2z+2z^2)Log[1-z]+2z(1-z)+Log[4e^2/mu2](1-2z+2z^2);
PqgAPR[z_]:=1-2z+2z^2;


(* simplificatin with e1=e2=em=e *)
PqgNNLO[z_,e_,mu2_]:=-69/4 + (255*z)/4 - 49*z^2 + (Pi^2*(17 - 26*z + 18*z^2))/12 + 
 (11*(1 - 2*z + 2*z^2)*Log[1 - z]^3)/6 + (-3/8 - z/2 + z^2)*Log[z]^2 + 
 (5/4 - (5*z)/2 + (7*z^2)/3)*Log[z]^3 + Log[(4*e^2)/mu2]^2*
  (5/4 - 2*z + 3*z^2 + (1 - 2*z + 2*z^2)*Log[1 - z] + 
   (-1/2 + z - 2*z^2)*Log[z]) + Log[1 - z]^2*
  (-7 + 21*z - 17*z^2 + (5/2 - 5*z + z^2)*Log[z]) - 
 3*PolyLog[2, 1 - z] + Log[(4*e^2)/mu2]*(8 - (41*z)/2 + 15*z^2 + 
   (Pi^2*(1 - 2*z + 2*z^2))/3 + (3 - 6*z + 6*z^2)*Log[1 - z]^2 + 
   (3/2 - 10*z + 10*z^2)*Log[z] + (1/2 - z)*Log[z]^2 + 
   Log[1 - z]*(-1 + 8*z - 4*z^2 - 4*z^2*Log[z]) + 
   (2 - 4*z)*PolyLog[2, 1 - z]) + 
 Log[z]*(1 + (35*z)/4 - 8*z^2 + Pi^2*(3/2 - 3*z + 3*z^2) + 
   (-1 + 2*z - 2*z^2)*PolyLog[2, 1 - z]) + 
 Log[1 - z]*(16 - (109*z)/2 + 44*z^2 + (Pi^2*(1 - 2*z + 2*z^2))/3 + 
   (3 - 14*z + 14*z^2)*Log[z] + (-7/2 + 7*z - 7*z^2)*Log[z]^2 + 
   (5 - 10*z + 2*z^2)*PolyLog[2, 1 - z]) + 
 (-9 + 18*z - 10*z^2)*PolyLog[3, 1 - z] + 
 (-9 + 18*z - 14*z^2)*PolyLog[3, z] + (18 - 36*z + 32*z^2)*Zeta[3];
 
 PaqNNLO[z_,e_,mu2_]:=73/2 - 27/z - (29*z)/4 + (Pi^2*(-38 + 30/z + 17*z))/12 + 
 (13*(2 - 2*z + z^2)*Log[1 - z]^3)/(6*z) + ((-3 + 5*z)*Log[z])/4 + 
 (5/2 + (17*z)/8)*Log[z]^2 + ((2 - z)*Log[z]^3)/12 + 
 Log[1 - z]^2*((58 - 42/z - 27*z)/4 + (4 - 2*z)*Log[z]) + 
 Log[(4*e^2)/mu2]^2*(1 - z/4 + (-2 + 2/z + z)*Log[1 - z] + 
   (1 - z/2)*Log[z]) + (8 + 5*z)*PolyLog[2, 1 - z] + 
 4*z*Log[2]*PolyLog[2, 1 - z] - z*Log[16]*PolyLog[2, 1 - z] + 
 Log[1 - z]*(-23 + 18/z - (3*z)/2 - (5*Pi^2*(2 - 2*z + z^2))/(3*z) + 
   (8 + 5*z)*Log[z] + (8 - 4*z)*PolyLog[2, 1 - z]) + 
 Log[(4*e^2)/mu2]*(-15/2 + 5/z - 2*z - (2*Pi^2*(2 - 2*z + z^2))/(3*z) + 
   (-6 + 6/z + 3*z)*Log[1 - z]^2 + (4 + (5*z)/2)*Log[z] + 
   (1 - z/2)*Log[z]^2 + Log[1 - z]*(10 - 6/z - 4*z + (4 - 2*z)*Log[z]) + 
   (4 - 2*z)*PolyLog[2, 1 - z]) + 4*(-2 + z)*PolyLog[3, 1 - z] + 
 2*(-2 + z)*PolyLog[3, z] + (-12 + 16/z + 6*z)*Zeta[3];


Gqg = 13/2 - Pi^2+3/2 Log[ec^2/mu2] + Log[ec/e5]^2 +
	+(3+2Log[ec/e5])Log[4 eta25] + 2 PolyLog[2,1-eta25]+
	-(3/2+2Log[ec/e5])Log[eta15/(1-eta15)]w56[1,1] + 
	-(3/2+2Log[ec/e5])Log[eta25/(1-eta25)]w56[2,2];


testresaqb = 
	ONNLO FLM[a,qb,l,lb,g,qb][1,2,3,4,5,6]+ONLO FLVfin[a,qb,l,lb,qb][QCD,1,2,3,4,5] + 
	(* -- *)
	qae xn Qq^2 PqgTNLO[z,e,mu2]FLVfin[q,qb,l,lb][QCD,z,2,3,4]+
	(* -- *)
	qae qas xn Qq^2 Cf PqgTNLO[z,e,mu2] FLM[q,qb,l,lb][z,zb,3,4]PqqTNLO[zb,e,mu2]+
	(* - *)
	qae qas xn Cf Qq^2 PqgNNLO[z,e,mu2] FLM[q,qb,l,lb][z,2,3,4]+
	qae qas Cf Qq^2 PaqNNLO[zb,e,mu2]FLM[a,a,l,lb][1,zb,3,4]+
	(* ONLO below *)
	qae Qq^2 xn (PqgTNLO[z,e,mu2]+w[1,1] Log[eta15/2]PqgAPR[z])*
	ONLO FLM[q,qb,l,lb,g][z,2,3,4,5] + 
	(* -- *)
	qas Cf (PqqTNLO[zb,e,mu2] + w[2,2] Log[eta25/2] PqqAPR[zb]) * 
	ONLO FLM[a,qb,l,lb,qb][1,zb,3,4,5] + 
	qas Cf Gqg ONLO FLM[a,qb,l,lb,qb][1,2,3,4,5];


nnlo[aqb]-testresaqb/.w56[2,2]->1-w56[1,1]/.
{e[1]->e,e[2]->e,em->e,ec->e,e[5]->e5}//PowerExpand//Expand


(* ::Section:: *)
(*\[Gamma] g \[Checkmark]*)


(* \[Gamma][1] g[2] -> l[3] lb[4] q[5] qb[6] *)


(* ::Subsection:: *)
(*RR*)


(* single collinear: use the A...D splitting to use Max's result as it is
	1 = w5161 + w5262 + w5162 + w5261 \[Rule] 
	w5161 + w5262[thA + thB + thC + thD] + w5261 + w5162 ->
	C51 [w5161 + w5162] + C52 [w5262 thC + w5261] +
	C61 [w5161 + w5261] + C62 [w5262 thA + w5162] = 
	C51 + C52[1-w[2,2](thC-1)] + C62 + C62[1-w[2,2](thA-1)]
*)

sc = qae q1^2 xn CiDi[qg,1,z] ONLO FLM[qb,g,l,lb,qb][z,2,3,4,5] + 
	 qas tr CiDi[qg,2,zb] ONLO FLM[a,qb,l,lb,qb][1,zb,3,4,5] (1-w[2,2](1-(eta25/2)^(-ep))) + 
	 qae q1^2 xn CiDi[qg,1,z] ONLO FLM[q,g,l,lb,q][z,2,3,4,5] + 
	 qas tr CiDi[qg,2,zb] ONLO FLM[a,q,l,lb,q][1,zb,3,4,5] (1-w[2,2](1-(eta25/2)^(-ep))) +
	 (* -- *)
	qae q1^2 xn CiDi[qg,1,z] qas tr Ci[qg,2,zb] FLM[qb,q,l,lb][z,zb,3,4] + 
	qas tr (CiDi[qg,2,zb] qae xn q1^2 Ci[qg,1,z] FLM[q,qb,l,lb][z,zb,3,4] + 
	 Nb/2 (1/2)^(-ep) qae q2^2 times[gq,qg,2,zb] FLM[a,a,l,lb][1,zb,3,4]) + 
	qae q1^2 xn CiDi[qg,1,z] qas tr Ci[qg,2,zb]FLM[q,qb,l,lb][z,zb,3,4] +
	qas tr (CiDi[qg,2,zb] qae xn q1^2 Ci[qg,1,z]FLM[qb,q,l,lb][z,zb,3,4] + 
	 Nb/2 (1/2)^(-ep) qae q2^2 times[gq,qg,2,zb]FLM[a,a,l,lb][1,zb,3,4]);


dc = -qas qae xn q1^2 tr CiDi[qg,1,z] CiDi[qg,2,zb] (FLM[q,qb,l,lb][z,zb,3,4] + FLM[qb,q,l,lb][z,zb,3,4]);


(* max does not include the (1-ep) averaging factors, which here are anyhow not needed *)
(* careful with the sign, no flip channel -> +1 *)
tc = qas qae/ep tr q1^2(4e[2]^2/mu2)^(-2ep) 2^(4 ep)*((+ep) ISab[zb,7]) FLM[a,a,l,lb][1,zb,3,4];


onnlo = ONNLO FLM[a,g,l,lb,q,qb][1,2,3,4,5,6];


rr = sc + dc + tc + onnlo;


(* ::Subsection:: *)
(*CV*)


P0xNLOEWK = conv[NLO[aqb],PAP0[qg]] + conv[NLO[aq],PAP0[qg]];
P0xNLOQCD = conv[PAP0[qg],NLO[qbg]] + conv[PAP0[qg],NLO[qg]];


cv = aemontwopi/ep P0xNLOQCD + asontwopi/ep P0xNLOEWK + 
(* -- *)
(* -- extra 2: qqb + qbq *)
-1/2 asontwopi aemontwopi/ep^2 (2*conv[cAP0[gq],cAP0[qg]][zb]) * (tr q1^2) FLM[a,a,l,lb][1,zb,3,4] + 
(* -- *)
+1/2 asontwopi aemontwopi/ep * (q1^2 tr PAP1[agAB,zb]FLM[a,a,l,lb][1,zb,3,4])+
(* -- *)
- asontwopi aemontwopi/ep^2 PAP0[qg,z] PAP0[qg,zb] * (xn q1^2 tr) (FLM[q,qb,l,lb][z,zb,3,4] + 
	FLM[qb,q,l,lb][z,zb,3,4]);


(* ::Subsection::Closed:: *)
(*Finite part check*)


mytmp = rr+cv/.{q1^2->Qq^2,q2^2->Qq^2}//.subs;


Series[mytmp,{ep,0,0}];
Normal[%]/.four->4//PowerExpand;
nnlo[ag] = Collect[%,FLM[a___][b___],Simplify];


(* ::Subsection:: *)
(*Simplified result*)


Paq2NNLO[z_,e_,mu2_]:=
(206 + 1056*z - 840*z^2 - 422*z^3)/(27*z) + (4*Pi^2*(-2 - 3*z - 3*z^2 + 2*z^3))/(9*z) + 
 ((8*Pi^2*(1 + z))/3 + (4*(-20 - 57*z + 39*z^2 + 38*z^3))/(9*z))*Log[1 - z] + 
 (2 + 6*z)*Log[z] + ((-5 - 11*z)*Log[z]^2)/2 + ((1 + z)*Log[z]^3)/3 + 
 Log[1 - z]^2*((4*(3 + 4/z - 3*z - 4*z^2))/3 - 8*(1 + z)*Log[z]) + 
 Log[(4*e^2)/mu2]^2*(1 + 4/(3*z) - z - (4*z^2)/3 + 2*(1 + z)*Log[z]) + 
 (4 + 12*z - 16*(1 + z)*Log[1 - z])*PolyLog[2, z] + 
 Log[(4*e^2)/mu2]*((4*Pi^2*(1 + z))/3 + (2*(-20 - 57*z + 39*z^2 + 38*z^3))/(9*z) + 
   (4*(3 + 4/z - 3*z - 4*z^2)*Log[1 - z])/3 + (-2 - 6*z)*Log[z] + 2*(1 + z)*Log[z]^2 - 
   8*(1 + z)*PolyLog[2, z]) - 16*(1 + z)*PolyLog[3, 1 - z] - 8*(1 + z)*PolyLog[3, z] + 
 8*(1 + z)*Zeta[3];


testresag = 
	ONNLO FLM[a,g,l,lb,q,qb][1,2,3,4,5,6] + 
(* -- ONLO below -- *)
ONLO qae Qq^2 xn (PqgTNLO[z,e,mu2])*
	(FLM[q,g,l,lb,q][z,2,3,4,5] + FLM[qb,g,l,lb,qb][z,2,3,4,5]) + 
(* - *)
ONLO qas tr (PqgTNLO[zb,e,mu2]+w[2,2] Log[eta25/2]PqgAPR[zb])*
	(FLM[a,q,l,lb,q][1,zb,3,4,5] + FLM[a,qb,l,lb,qb][1,zb,3,4,5]) + 
(* -- FLM12 below -- *)
qae qas tr xn Qq^2 * PqgTNLO[z,e,mu2] PqgTNLO[zb,e,mu2] * 
	(FLM[q,qb,l,lb][z,zb,3,4] + FLM[qb,q,l,lb][z,zb,3,4]) + 
(* - *)
qae qas tr Qq^2 Paq2NNLO[zb,e,mu2] FLM[a,a,l,lb][1,zb,3,4];


nnlo[ag]-  testresag;
%/.{e[1]->e,e[2]->e};
%/.Log[4 e^2/mu2]:>2Log[2]+2Log[e]-Log[mu2];
Coefficient[%,FLM[a,a,l,lb][1,zb,3,4]]//Expand//Simplify
