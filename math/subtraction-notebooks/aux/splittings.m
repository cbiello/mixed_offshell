(* ::Package:: *)

(* Splitting functions for Nested soft collinear subtraction *)

splittings = {

(* LO Altarelli-Parisi splitting functions *)    
PqqAP0[z_] :>  2 D0[z] - (1 + z) + 3/2 delta[1 - z],
PgqAP0[z_] :> z^2+(1-z)^2,
PqgAP0[z_] :> (1+(1-z)^2)/z,
PggAP0[z_] :>  2 CA (D0[z] + 1/z + z(1-z) - 2 ) + \[Beta]0 delta[1 - z],


PgqNLOIS[z_]:> (1-2 z+2 z^2)+(-1+(1-2 z+2 z^2) (1-2 Log[1-z])) \[Epsilon]+(-1+2 Log[1-z]+(1-2 z+2 z^2) (1-2 Log[1-z]+2 Log[1-z]^2)) \[Epsilon]^2+(-1+2 Log[1-z]-2 Log[1-z]^2+(1-2 z+2 z^2) (1-2 Log[1-z]+2 Log[1-z]^2-4/3 Log[1-z]^3)) \[Epsilon]^3+ \[Epsilon]^4 ORD[PqgNLOis,4],
PggNLOIS[z_, L_]:> 2 (-2+1/z+(1-z) z)+2 D0[z]-2 L delta[1-z]+\[Epsilon] (-4 D1[z]+2 L^2 delta[1-z]-4 (-2+1/z+(1-z) z) Log[1-z])+\[Epsilon]^2 (4 D2[z]-4/3 L^3 delta[1-z]+4 (-2+1/z+(1-z) z) Log[1-z]^2)+\[Epsilon]^3 (-((8 D3[z])/3)+2/3 L^4 delta[1-z]-8/3 (-2+1/z+(1-z) z) Log[1-z]^3)+\[Epsilon]^4 ((4 D4[z])/3-4/15 L^5 delta[1-z]+4/3 (-2+1/z+(1-z) z) Log[1-z]^4)+\[Epsilon]^5 ORD[PggNLOis,5],
PqgNLOIS[z_]:> (2-2 z+z^2)/z+(-z-(2 (2-2 z+z^2) Log[1-z])/z) \[Epsilon]+(2 z Log[1-z]+(2 (2-2 z+z^2) Log[1-z]^2)/z) \[Epsilon]^2+(-2 z Log[1-z]^2-(4 (2-2 z+z^2) Log[1-z]^3)/(3 z)) \[Epsilon]^3+(4/3 z Log[1-z]^3+(2 (2-2 z+z^2) Log[1-z]^4)/(3 z)) \[Epsilon]^4+ \[Epsilon]^5 ORD[PqgNLOis,5],

(* splitting functions *)
PqqNLOIS[z_, L_] :> (2 D0[z] - (1 + z) - 2 L delta[1 - z] + ep (2 L^2 delta[1 - z] - 4 D1[z] + 2 (1 + z) Log[1 - z] - (1 - z)) + ep^2 (4 D2[z] - 4/3 L^3 delta[1 - z] + 2 (1 - z) Log[1 - z] - 2 (1 + z) Log[1 - z]^2) + ep^3 (2/3 L^4 delta[1 - z] - 8/3 D3[z] - 2 (1 - z) Log[1 - z]^2 + 4/3 (1 + z) Log[1 - z]^3) + ep^4 ORD[PqqNLOIS,4]),
    
Pqqfin[z_,E1_] :> -(1/\[Epsilon]) (Gamma[1 - \[Epsilon]]^2/
          Gamma[1 - 2 \[Epsilon]] ((4 E1^2)/
           mu2)^-\[Epsilon] (PqqNLOIS[z, 0] + 3/2 delta[1 - z]) -
         Gamma[1 - \[Epsilon]]/Exp[\[Epsilon] EulerGamma]  PqqAP0[z]),
         
Pgqfin[z_,E1_]:> -(1/\[Epsilon])(Gamma[1-\[Epsilon]]^2/Gamma[1-2\[Epsilon]] ((4 E1^2)/mu2)^-\[Epsilon] PgqNLOIS[z]-Gamma[1-\[Epsilon]]/Exp[\[Epsilon] EulerGamma]  PgqAP0[z]),
Pggfin[z_,E1_]:> -(1/\[Epsilon])(Gamma[1-\[Epsilon]]^2/Gamma[1-2\[Epsilon]] ((4 E1^2)/mu2)^-\[Epsilon] (CA PggNLOIS[z,0]+\[Beta]0 delta[1 - z])-Gamma[1-\[Epsilon]]/Exp[\[Epsilon] EulerGamma]  PggAP0[z]),
    
PggNLOFS[L1_]:>(11/6+2 L1)+(67/9-2 L1^2-(2 \[Pi]^2)/3) ep +1/27 (808+36 L1^3-33 \[Pi]^2-432 Zeta[3]) ep^2 -2/405 (-24280+135 L1^4+1005 \[Pi]^2+81 \[Pi]^4+5940 Zeta[3]) ep^3 + ep^4 ORD[PggNLOFS,4],
PqgNLOFS[L1_]:> -(4/3)-(46 \[Epsilon])/9+1/27 (-553+24 \[Pi]^2) \[Epsilon]^2+4/81 (-1660+69 \[Pi]^2+432 Zeta[3]) \[Epsilon]^3+(2 (-199240+8295 \[Pi]^2+324 \[Pi]^4+49680 Zeta[3]) \[Epsilon]^4)/1215 + \[Epsilon]^5 ORD[PqgNLOFS,4],

(* energy-ordered splittings *)

(* initial state *)
PqqNLOISeo[z_, L_] :> ((-1 + 2 D0[z] - z - 2 L delta[1 - z]) + (-1 + z + 4 L^2 delta[1 - z] - 8 D1[z] - 4 (-1 - z) Log[1 - z]) ep + (-(16/3) L^3 delta[1 - z] - 4 (-1 + z) Log[1 - z] + 16 D2[z] + 8 (-1 - z) Log[1 - z]^2) ep^2 + (16/3 L^4 delta[1 - z] + 8 (-1 + z) Log[1 - z]^2 - 64/3 D3[z] - 32/3 (-1 - z) Log[1 - z]^3) ep^3 + ep^4 ORD[PqqNLOISeo,4]),
PggNLOISeo[z_,L_] :> (2 (-2+1/z+(1-z) z)+2 D0[z]-2 L delta[1-z])+(-8 D1[z]+4 L^2 delta[1-z]-8 (-2+1/z+(1-z) z) Log[1-z]) \[Epsilon]+(16 D2[z]-16/3 L^3 delta[1-z]+16 (-2+1/z+(1-z) z) Log[1-z]^2) \[Epsilon]^2+(-((64 D3[z])/3)+16/3 L^4 delta[1-z]-64/3 (-2+1/z+(1-z) z) Log[1-z]^3) \[Epsilon]^3+((64 D4[z])/3-64/15 L^5 delta[1-z]+64/3 (-2+1/z+(1-z) z) Log[1-z]^4) \[Epsilon]^4+(-((256 D5[z])/15)+128/45 L^6 delta[1-z]-256/15 (-2+1/z+(1-z) z) Log[1-z]^5) \[Epsilon]^5+ \[Epsilon]^6 ORD[PggNLOiseo,6],
PgqNLOISeo[z_]:> (1+2 (-1+z) z)+(2 (-1+z) z-4 (1+2 (-1+z) z) Log[1-z]) \[Epsilon]+(2 (-1+z) z-8 (-1+z) z Log[1-z]+8 (1+2 (-1+z) z) Log[1-z]^2) \[Epsilon]^2+(2 (-1+z) z-8 (-1+z) z Log[1-z]+16 (-1+z) z Log[1-z]^2-32/3 (1+2 (-1+z) z) Log[1-z]^3) \[Epsilon]^3+(2 (-1+z) z-8 (-1+z) z Log[1-z]+16 (-1+z) z Log[1-z]^2-64/3 (-1+z) z Log[1-z]^3+32/3 (1+2 (-1+z) z) Log[1-z]^4) \[Epsilon]^4+(2 (-1+z) z-8 (-1+z) z Log[1-z]+16 (-1+z) z Log[1-z]^2-64/3 (-1+z) z Log[1-z]^3+64/3 (-1+z) z Log[1-z]^4-128/15 (1+2 (-1+z) z) Log[1-z]^5) \[Epsilon]^5 + \[Epsilon]^6 ORD[PgqNLOISEO] ,

(* final state *)
PggNLOFSeo[L_] :> (11/6 + 2 L) + (157/12 - 4 L^2 - (2 Pi^2)/3) ep + (2039/24 + (16 L^3)/3 - (22 Pi^2)/9 - 24 Zeta[3]) ep^2 + (8471/16 - (16 L^4)/3 - 1/45 Pi^2 (785 + 44 Pi^2) - 88 Zeta[3]) ep^3 + ep^4 ORD[PggNLOFSeo,4],
PqqNLOFSeo[L_] :> 17/6+2 L+1/3 (55-12 L^2-2 \[Pi]^2) \[Epsilon]+73 \[Epsilon]^2+(16 L^3 \[Epsilon]^2)/3-(34 \[Pi]^2 \[Epsilon]^2)/9+(4153 \[Epsilon]^3)/9-(16 L^4 \[Epsilon]^3)/3-(220 \[Pi]^2 \[Epsilon]^3)/9-(44 \[Pi]^4 \[Epsilon]^3)/45+(24220 \[Epsilon]^4)/9+(64 L^5 \[Epsilon]^4)/15-(292 \[Pi]^2 \[Epsilon]^4)/3-(748 \[Pi]^4 \[Epsilon]^4)/135+251/18 \[Epsilon]^2 (3+17 \[Epsilon]+(110-4 \[Pi]^2) \[Epsilon]^2)+1/2 \[Epsilon]^2 (-3-17 \[Epsilon]+2 (-55+2 \[Pi]^2) \[Epsilon]^2)+\[Epsilon]^5 ORD[PqqNLOfseo,5]+8 \[Epsilon]^2 (-3-17 \[Epsilon]+2 (-55+2 \[Pi]^2) \[Epsilon]^2) Zeta[3]-672 \[Epsilon]^4 Zeta[5],   
PqgNLOFSeo[L_] :> -(4/3)-(22 \[Epsilon])/3+(-(127/3)+(16 \[Pi]^2)/9) \[Epsilon]^2+(-249+(88 \[Pi]^2)/9+64 Zeta[3]) \[Epsilon]^3+(-1479+(508 \[Pi]^2)/9+(352 \[Pi]^4)/135+352 Zeta[3]) \[Epsilon]^4+ ep^5 ORD[PqgNLOFSeo,5],
        
(* convolutions *)
PqqNLOCV[z_] :>
      6 D0[z] + 8 D1[z] + (9/4 - (2 \[Pi]^2)/3) delta[1 - z] + \[Epsilon]*(4/3 \[Pi]^2 D0[z] - 6 D1[z] - 12 D2[z] -
     8 Zeta[3] delta[1 - z]) + \[Epsilon]^2*(16 Zeta[3] D0[z] -
     1/3 8 \[Pi]^2 D1[z] + 6 D2[z] + (32 D3[z])/3 -
     8/45 \[Pi]^4 delta[1 - z]) + (-(((3 z^2 + 1) Log[z])/(1 - z)) - z - 4 (z + 1) Log[1 - z] -5) + \[Epsilon]*((z + 1) (2 PolyLog[2, z] + 6 Log[1 - z]^2 -
     3 Log[z]^2 - \[Pi]^2) - 3/2 (1 - z) + (
     4 Log[z]^2 - 6 Log[z])/(1 - z) + (z + 5) Log[1 - z] +
     2 (z + 3) Log[
     z]) + \[Epsilon]^2*((z + 1) (-4 PolyLog[3, 1 - z] -
     4 PolyLog[3, z] + (2 \[Pi]^2 - 4 PolyLog[2, z]) Log[1 - z] -
     16/3 Log[1 - z]^3 - (2 Log[z] + 3) Log[1 - z]^2 +
     2/3 \[Pi]^2 Log[z] - 4 Zeta[3]) - (2 (1 + 3 z^2))/(3 (1 - z))
     Log[z]^3 + (6 - 2 z) PolyLog[2, z] + \[Pi]^2 (z - 5/3) -
     2 (z + 3) Log[z]^2 + (6 Log[z]^2)/(1 - z) -
     2 Log[z] + (1 - z) (3 Log[1 - z] - 2)) + \[Epsilon]^3*
    ORD[PqqNLOCV,3],
    

(* PqqRV2 from 1702.01352 Eq.7.7 and A.11 *)
PqqRV2[z_, L1_] :> -CA CF (1/ep^2 (delta[1 - z] L1 - D0[z]) + 2/ep (2 D1[z] - delta[1 - z] L1^2) + (-(1/3) delta[1 - z] L1 (Pi^2 - 8 L1^2) + Pi^2/3 D0[z] - 8 D2[z]) + (2/3 delta[1 - z] L1^2 (Pi^2 - 4 L1^2) - 4/3 Pi^2 D1[z] + 32/3 D3[z]) ep + (1 + z)/(2 ep^2) +
1/(2 ep) ((1 - z) - 4 (1 + z) Log[1 - z]) + (-(z/2) - 2 (1 - z) Log[1 - z] + (1 + z) (4 Log[1 - z]^2 - Pi^2/6) - 
((1 + z^2) PolyLog[2, 1 - z])/(1 - z)) + (3/2 z Log[1 - z] + (1 + z)/6 (-3 + 4 Pi^2 Log[1 - z] - 32 Log[1 - z]^3) + 
(1 - z) (4 Log[1 - z]^2 - Pi^2/6 + PolyLog[2, 1 - z]) + 3 (1 + z^2)/(1 - z) Log[1 - z] PolyLog[2, 1 - z]) ep) - 
CF^2 (-(1/ep) (((1 + z^2) Log[z])/(1 - z)) + (z/2 + (1 - z) Log[z] + 
(1 + z^2)/(1 - z) (3 Log[1 - z] Log[z] + PolyLog[2, 1 - z])) + 
(1/2 (1 + z - 3 z Log[1 - z]) - (1 - z) (3 Log[1 - z] Log[z] + PolyLog[2, 1 - z]) - 
1/2 ((1 + z^2)/(1 - z)) (9 Log[1 - z]^2 Log[z] + 6 Log[1 - z] PolyLog[2, 1 - z] - 2 PolyLog[3, 1 - z])) ep) + ep^2 * ORD[PqqRV2,2],
    
PggNLOFS42[L3_] :> (11/6 +2 L3) + (37/4 - 2 L3^2 - (4 \[Pi]^2)/3) ep + (1201/24 + (4 L3^3)/
          3 - (22 \[Pi]^2)/9 - 48 Zeta[3]) ep^2 + (4537/16 - (2 L3^4)/
          3 - (37 \[Pi]^2)/3 - (88 \[Pi]^4)/45 -
          88 Zeta[3]) ep^3 + (52827/32 + (4 L3^5)/15 - (484 \[Pi]^4)/
          135 - 444 Zeta[3] + \[Pi]^2 (-(1201/18) + 64 Zeta[3]) -
          1344 Zeta[5]) ep^4 + ORD[PggNloFs42,5] ep^5,
          
PqqNLOFS42[L3_] :>  (17/6+2 L3)+1/3 (49-6 L3^2-4 \[Pi]^2) \[Epsilon]+2/9 (429+6 L3^3-17 \[Pi]^2-216 Zeta[3]) \[Epsilon]^2-2/45 (-12690+15 L3^4+490 \[Pi]^2+44 \[Pi]^4+3060 Zeta[3]) \[Epsilon]^3+4/135 (113400+9 L3^5-4290 \[Pi]^2-187 \[Pi]^4-26460 Zeta[3]+2160 \[Pi]^2 Zeta[3]-45360 Zeta[5]) \[Epsilon]^4 + ORD[PqqNloFs42,5] ep^5,  

PggNLOFSRV[L_] :> 11/6 + 2 L + (67/9 - 2 L^2 - (2 Pi^2)/3) ep + ep^2 (808/27 + (4 L^3)/3 - (11 Pi^2)/9 - 16 Zeta[3]) + ep^3 (9712/81 - (2 L^4)/3 - (134 Pi^2)/27 - (2 Pi^4)/5 - (88 Zeta[3])/3) +  ep^4 * ORD[PggNLOFSRV,4],

\[Gamma]1loopggg :> 11/6 CA^2 1/ep^2 + CA^2/ep (134/9 - (4 Pi^2)/3) + CA^2 (1013/9 - (44 Pi^2)/9 - 60 Zeta[3]) - (CA Nf)/6 + ep (CA^2 (14635/18 - (335 Pi^2)/9 - (106 Pi^4)/45 - 550/3 Zeta[3]) - 31/18 CA Nf) + ORD[\[Gamma]1Loopggg,2]  ep^2,

\[Gamma]1loopqqg :> (17 CA CF)/(6 \[Epsilon]^2)+(5/18 CF (66 CA+17 CF)-2/3 CF (CA+CF) \[Pi]^2)/\[Epsilon]-1/18 CF (6 CA (-347+17 \[Pi]^2+48 Zeta[3])+CF (-1041+34 \[Pi]^2+792 Zeta[3]))-1/90 (CF (5 CF (-9521+366 \[Pi]^2+42 \[Pi]^4+2244 Zeta[3])+CA (-64995+3360 \[Pi]^2+2 \[Pi]^4+14280 Zeta[3]))) \[Epsilon] + ORD[\[Gamma]1Loopqqg,2]  ep^2,

\[Gamma]1loopqgq :> -((4 (CA CF))/(3 \[Epsilon]^2))-(5 (CF (15 CA+4 CF)))/(9 \[Epsilon])+(-(1/6) CF (316 CA+161 CF)+8/9 CF (3 CA+CF) \[Pi]^2)+1/9 CF (3 CA (-995+51 \[Pi]^2+224 Zeta[3])+CF (-2173+84 \[Pi]^2+528 Zeta[3])) \[Epsilon] + ORD[\[Gamma]1Loopqgq,2]  ep^2,     
PqqISminus2ep[z_] :> - (1/\[Epsilon]) delta[1 - z] + (-1 - z +
         2 D0[z]) + (-1 + z - 4 D1[z] -
          2 (-1 - z) Log[1 - z]) \[Epsilon] + (4 D2[z] -
          2 (-1 + z) Log[1 - z] +
          2 (-1 - z) Log[1 - z]^2) \[Epsilon]^2 + (-(8/3) D3[z] +
          2 (-1 + z) Log[1 - z]^2 -
          4/3 (-1 - z) Log[1 - z]^3) \[Epsilon]^3 + (4/3 D4[z] -
          4/3 (-1 + z) Log[1 - z]^3 +
          2/3 (-1 - z) Log[1 - z]^4) \[Epsilon]^4 + ORD[PqqISminus2ep,5] \[Epsilon]^5,
    
PqqISminus4ep[z_] :> - (1/(2 \[Epsilon])) delta[1 - z] + (-1 - z +
         2 D0[z]) + (-1 + z - 8 D1[z] -
          4 (-1 - z) Log[1 - z]) \[Epsilon] + (16 D2[z] -
          4 (-1 + z) Log[1 - z] +
          8 (-1 - z) Log[1 - z]^2) \[Epsilon]^2 + (-(64/3) D3[z] +
          8 (-1 + z) Log[1 - z]^2 -
          32/3 (-1 - z) Log[1 - z]^3) \[Epsilon]^3 + (64/3 D4[z] -
          32/3 (-1 + z) Log[1 - z]^3 +
          32/3 (-1 - z) Log[1 - z]^4) \[Epsilon]^4 + ORD[PqqISminus4ep,5] \[Epsilon]^5,
          
          
PggISminus2ep[z_]:> 2 (-2+1/z+(1-z) z)+2 D0[z]-delta[1-z]/\[Epsilon]+\[Epsilon] (-4 D1[z]-4 (-2+1/z+(1-z) z) Log[1-z])+\[Epsilon]^2 (4 D2[z]+4 (-2+1/z+(1-z) z) Log[1-z]^2)+\[Epsilon]^3 (-((8 D3[z])/3)-8/3 (-2+1/z+(1-z) z) Log[1-z]^3)+\[Epsilon]^4 ((4 D4[z])/3+4/3 (-2+1/z+(1-z) z) Log[1-z]^4)+ ORD[PggISminus2ep,5] \[Epsilon]^5,
PggISminus4ep[z_]:> 2 (-2+1/z+(1-z) z)+2 D0[z]-delta[1-z]/(2 \[Epsilon])+\[Epsilon] (-8 D1[z]-8 (-2+1/z+(1-z) z) Log[1-z])+\[Epsilon]^2 (16 D2[z]+16 (-2+1/z+(1-z) z) Log[1-z]^2)+\[Epsilon]^3 (-((64 D3[z])/3)-64/3 (-2+1/z+(1-z) z) Log[1-z]^3)+\[Epsilon]^4 ((64 D4[z])/3+64/3 (-2+1/z+(1-z) z) Log[1-z]^4)+ ORD[PggISminus4ep,5] \[Epsilon]^5,
PqqNLOFS[L_] :>17/6+2 L+(209/18-2 L^2-(2 \[Pi]^2)/3) \[Epsilon]+\[Epsilon]^3 (15064/81-(2 L^4)/3-(209 \[Pi]^2)/27-(2 \[Pi]^4)/5-(136 Zeta[3])/3)+\[Epsilon]^2 (1255/27+(4 L^3)/3-(17 \[Pi]^2)/9-16 Zeta[3])+\[Epsilon]^4 (180784/243+(4 L^5)/15-(2510 \[Pi]^2)/81-(17 \[Pi]^4)/15-(1672 Zeta[3])/9+32/3 \[Pi]^2 Zeta[3]-192 Zeta[5])+\[Epsilon]^5 ORD[PqqNloFs,5],
          
(* gamma functions with Z factor *)    
\[Gamma]22zggg :> 11/6 + 1/9 (67 - 6 \[Pi]^2) \[Epsilon] + (808/27 - (11 \[Pi]^2)/9 -
          16 Zeta[3]) \[Epsilon]^2 -
       2/405 (-24280 + 1005 \[Pi]^2 + 81 \[Pi]^4 +
          5940 Zeta[3]) \[Epsilon]^3 + ORD[\[Gamma]22Zggg,4] \[Epsilon]^4,
    
\[Gamma]22zgqq :> -(1/3) - (23 \[Epsilon])/18 + (-(281/54) + (2 \[Pi]^2)/9) \[Epsilon]^2 +
    1/162 (-3419 + 138 \[Pi]^2 + 864 Zeta[3]) \[Epsilon]^3 + (-(41297/486) + (
       281 \[Pi]^2)/81 + (2 \[Pi]^4)/15 + (184 Zeta[3])/9) \[Epsilon]^4 + ORD[\[Gamma]22Zgqq,5] \[Epsilon]^5,
    
\[Gamma]24zggg :> 11/6 + (3/4 (157 - 8 \[Pi]^2) +
          2/3 (-157 + 8 \[Pi]^2)) \[Epsilon] + (2039/24 - (22 \[Pi]^2)/
          9 - 24 Zeta[3]) \[Epsilon]^2 + (8471/16 - (157 \[Pi]^2)/9 - (
          44 \[Pi]^4)/45 - 88 Zeta[3]) \[Epsilon]^3 +
        ORD[\[Gamma]24Zggg,4] \[Epsilon]^4,
        
\[Gamma]22zqqg :> 17/6+(209/18-(2 \[Pi]^2)/3) \[Epsilon]+\[Epsilon]^3 (15064/81-(209 \[Pi]^2)/27-(2 \[Pi]^4)/5-(136 Zeta[3])/3)+\[Epsilon]^2 (1255/27-(17 \[Pi]^2)/9-16 Zeta[3])+\[Epsilon]^4 (180784/243-(2510 \[Pi]^2)/81-(17 \[Pi]^4)/15-(1672 Zeta[3])/9+32/3 \[Pi]^2 Zeta[3]-192 Zeta[5]) +\[Epsilon]^5 ORD[\[Gamma]22zQQG,6] ,
\[Gamma]24zqqg :>17/6+(55-2 \[Pi]^2+2/3 (-55+2 \[Pi]^2)) \[Epsilon]+(340/3-(34 \[Pi]^2)/9-24 Zeta[3]) \[Epsilon]^2+(690-(220 \[Pi]^2)/9-(44 \[Pi]^4)/45-136 Zeta[3]) \[Epsilon]^3+(23725/9-(286 \[Pi]^2)/3-(748 \[Pi]^4)/135-251/9 (-55+2 \[Pi]^2)-2 (440-16 \[Pi]^2) Zeta[3]-672 Zeta[5]) \[Epsilon]^4+(57007/3-(16306 \[Pi]^2)/27-(968 \[Pi]^4)/27-(1936 \[Pi]^6)/945-251/27 (-657+34 \[Pi]^2)+(-1936-2 (1752-(272 \[Pi]^2)/3)) Zeta[3]+576 Zeta[3]^2-3808 Zeta[5]) \[Epsilon]^5+\[Epsilon]^6 ORD[\[Gamma]24zQQG,6] ,
\[Gamma]04zggg :> 11/6+203/18 ep +2909/54 ep^2 +38015/162 ep^3 +476141/486 ep^4 +5839463/1458 ep^5 + ORD[\[Gamma]04Zggg,6] ep^6,
\[Gamma]22zqgq :> -(4/3)-(46 \[Epsilon])/9+(-(553/27)+(8 \[Pi]^2)/9) \[Epsilon]^2+\[Epsilon]^3 (-(6640/81)+(92 \[Pi]^2)/27+(64 Zeta[3])/3)+\[Epsilon]^4 (-(79696/243)+(1106 \[Pi]^2)/81+(8 \[Pi]^4)/15+(736 Zeta[3])/9)+ ORD[\[Gamma]22Zqgq,5] \[Epsilon]^5,
\[Gamma]24zgqq :> -(1/3)-(29 \[Epsilon])/12+(-(131/8)+(4 \[Pi]^2)/9) \[Epsilon]^2+ \[Epsilon]^3 (-(1677/16)+(29 \[Pi]^2)/9+16 Zeta[3])+ \[Epsilon]^4 (-(20799/32)+(131 \[Pi]^2)/6+(88 \[Pi]^4)/135+116 Zeta[3])+ ORD[\[Gamma]24Zgqqb,5] \[Epsilon]^5,
\[Gamma]42zqqg :>17/6+(49-4 \[Pi]^2+2/3 (-49+4 \[Pi]^2)) \[Epsilon]+(286/3-(34 \[Pi]^2)/9-48 Zeta[3]) \[Epsilon]^2+(564-(196 \[Pi]^2)/9-(88 \[Pi]^4)/45-136 Zeta[3]) \[Epsilon]^3+(3360-(1144 \[Pi]^2)/9-(748 \[Pi]^4)/135+(-784+64 \[Pi]^2) Zeta[3]-1344 Zeta[5]) \[Epsilon]^4+ ORD[\[Gamma]42Zqqg,5] \[Epsilon]^5,
\[Gamma]42zqgq :> -(4/3)-(25 \[Epsilon])/3+1/9 (-471+16 \[Pi]^2) \[Epsilon]^2+1/9 \[Epsilon]^3 (-2889+100 \[Pi]^2+576 Zeta[3])+1/135 \[Epsilon]^4 (-262845+9420 \[Pi]^2+352 \[Pi]^4+54000 Zeta[3])+ ORD[\[Gamma]42Zqgq,5] \[Epsilon]^5,
\[Gamma]42zgqq :> -(1/3)-(19 \[Epsilon])/12+(-(69/8)+(4 \[Pi]^2)/9) \[Epsilon]^2+\[Epsilon]^3 (-(787/16)+(19 \[Pi]^2)/9+16 Zeta[3])+\[Epsilon]^4 (-(9185/32)+(23 \[Pi]^2)/2+(88 \[Pi]^4)/135+76 Zeta[3]),
\[Gamma]24zqgq :> -(4/3)-(22 \[Epsilon])/3+(5/9 (381-16 \[Pi]^2)+2/3 (-381+16 \[Pi]^2)) \[Epsilon]^2+(-249+(88 \[Pi]^2)/9+64 Zeta[3]) \[Epsilon]^3+(-1479+(508 \[Pi]^2)/9+(352 \[Pi]^4)/135+352 Zeta[3]) \[Epsilon]^4+ ORD[\[Gamma]24Zqgq,5] \[Epsilon]^5,

(* gamma functions without Z factor *) 
\[Gamma]22qqg :> 3/2+1/6 (39-4 \[Pi]^2) \[Epsilon]+(26-\[Pi]^2-16 Zeta[3]) \[Epsilon]^2+(104-(13 \[Pi]^2)/3-(2 \[Pi]^4)/5-24 Zeta[3]) \[Epsilon]^3+(-((3 \[Pi]^4)/5)+4/3 \[Pi]^2 (-13+8 Zeta[3])-8 (-52+13 Zeta[3]+24 Zeta[5])) \[Epsilon]^4+\[Epsilon]^5 ORD[\[Gamma]22QQG,5],
\[Gamma]22gqq  :> -(2/3)-(23 \[Epsilon])/9+1/27 (-281+12 \[Pi]^2) \[Epsilon]^2+\[Epsilon]^4 (-(41297/243)+(562 \[Pi]^2)/81+(4 \[Pi]^4)/15+(368 Zeta[3])/9)+1/81 \[Epsilon]^3 (-3419+138 \[Pi]^2+864 Zeta[3])+\[Epsilon]^4 ORD[\[Gamma]22GQQ,4],
\[Gamma]22ggg :> 11/3+2/9 (67-6 \[Pi]^2) \[Epsilon]-2/27 (-808+33 \[Pi]^2+432 Zeta[3]) \[Epsilon]^2-4/405 (-24280+1005 \[Pi]^2+81 \[Pi]^4+5940 Zeta[3]) \[Epsilon]^3+(233216/243-(3232 \[Pi]^2)/81-(22 \[Pi]^4)/15+32/9 (-67+6 \[Pi]^2) Zeta[3]-384 Zeta[5]) \[Epsilon]^4+\[Epsilon]^5 ORD[\[Gamma]22GGG,5],
\[Gamma]22qgq :> 3/2+1/6 (39-4 \[Pi]^2) \[Epsilon]+(26-\[Pi]^2-16 Zeta[3]) \[Epsilon]^2+(104-(13 \[Pi]^2)/3-(2 \[Pi]^4)/5-24 Zeta[3]) \[Epsilon]^3+(-((3 \[Pi]^4)/5)+4/3 \[Pi]^2 (-13+8 Zeta[3])-8 (-52+13 Zeta[3]+24 Zeta[5])) \[Epsilon]^4,

(* splittings for sectors B and D *)  
\[Gamma]perp22ggg :> 11/3+(131/9-(4 \[Pi]^2)/3) ep+ep^2 (1604/27-(22 \[Pi]^2)/9-32 Zeta[3])+ep^3 ORD[\[Gamma]Perp22ggg,3],
\[Gamma]perptilde22ggg :> -(2/3)-14/9 ep +(-(164/27)+(4 \[Pi]^2)/9) ep^2+ep^3 ORD[\[Gamma]Perptilde22ggg,3],
PqqSpinISminus4ep[z_]:> 1/2 (-3-z)+2 (D0[z]-4 ep D1[z]+8 ep^2 D2[z]-32/3 ep^3 D3[z]+32/3 ep^4 D4[z]-delta[1-z]/(4 ep))+2 ep (3+z) Log[1-z]-4 ep^2 (3+z) Log[1-z]^2+16/3 ep^3 (3+z) Log[1-z]^3-16/3 ep^4 (3+z) Log[1-z]^4+ ep^5 ORD[PqqspinISminus4ep,5],deltaperpg :> 13/36-Log[2]/3+ep ORD[deltaPerpg,1],
deltag :> -(131/72)+\[Pi]^2/6+(11 Log[2])/6+ ep (-(1541/216)+(11 \[Pi]^2)/18-(67 Log[2])/9+2/3 \[Pi]^2 Log[2]+(22 Log[2]^2)/3+4 Zeta[3])+ep^2 (-(9607/324)+(125 \[Pi]^2)/216+(7 \[Pi]^4)/45-(1487 Log[2])/54+55/18 \[Pi]^2 Log[2]-(137 Log[2]^2)/9+4/3 \[Pi]^2 Log[2]^2+(44 Log[2]^3)/3+(77 Zeta[3])/6+16 Log[2] Zeta[3])+ ep^3 ORD[\[Delta]g,3],
\[Gamma]24rzggg :> 29/12+(371/24-(2 \[Pi]^2)/3) ep+(1559/16-(29 \[Pi]^2)/9-24 Zeta[3]) ep^2+(19221/32-1/90 \[Pi]^2 (1855+88 \[Pi]^2)-116 Zeta[3]) ep^3+ep^4 ORD[\[Gamma]24Rzggg,5],

(* splittings for sectors B and D gq channel *) 
\[Gamma]perptilde22gqq :> 2/3+(20 \[Epsilon])/9-4/27 (-56+3 \[Pi]^2) \[Epsilon]^2-8/81 \[Epsilon]^3 (-328+15 \[Pi]^2+108 Zeta[3])-(4 \[Epsilon]^4 (1680 \[Pi]^2+81 \[Pi]^4+80 (-488+135 Zeta[3])))/1215+ep^4 ORD[\[Gamma]Perptilde22gqq,5],
\[Gamma]perptilde22Rgqq :> 1/3+(7 \[Epsilon])/9+1/27 (73-6 \[Pi]^2) \[Epsilon]^2+(\[Epsilon]^4 (51605-2190 \[Pi]^2-162 \[Pi]^4-15120 Zeta[3]))/1215+1/81 \[Epsilon]^3 (859-42 \[Pi]^2-432 Zeta[3])+ep^4 ORD[\[Gamma]Perptilde22Rgqq,5],
deltaq :> 2^(4\[Epsilon]) (23/72-Log[2]/3+\[Epsilon](103/108-\[Pi]^2/9+35/36 Log[2]-(2Log[2]^2)/3)+\[Epsilon]^2 (373/81-5/216 \[Pi]^2+67/54 Log[2]-\[Pi]^2/3 Log[2]+47/36 Log[2]^2-2/3 Log[2]^3-7/3 Zeta[3]))+\[Epsilon]^3 ORD[deltaQ,3],
deltaperpq :> -13/36+Log[2]/3 + \[Epsilon] ORD[deltaQ,1],
PggSpinISminus4ep[z_]:>(1-3 z+z^2-z^3)/z+2 (D0[z]-4 \[Epsilon] D1[z]+8 \[Epsilon]^2 D2[z]-32/3 \[Epsilon]^3 D3[z]+32/3 \[Epsilon]^4 D4[z]-128/15 \[Epsilon]^5 D5[z])-delta[1-z]/(2 \[Epsilon])+2 \[Epsilon] ((1-z+z^2-z^3)/(2 z)-(2 (1-3 z+z^2-z^3) Log[1-z])/z)+2 \[Epsilon]^2 ((1-z+z^2-z^3)/(2 z)-(2 (1-z+z^2-z^3) Log[1-z])/z+(4 (1-3 z+z^2-z^3) Log[1-z]^2)/z)+2 \[Epsilon]^3 ((1-z+z^2-z^3)/(2 z)-(2 (1-z+z^2-z^3) Log[1-z])/z+(4 (1-z+z^2-z^3) Log[1-z]^2)/z-(16 (1-3 z+z^2-z^3) Log[1-z]^3)/(3 z))+2 \[Epsilon]^4 ((1-z+z^2-z^3)/(2 z)-(2 (1-z+z^2-z^3) Log[1-z])/z+(4 (1-z+z^2-z^3) Log[1-z]^2)/z-(16 (1-z+z^2-z^3) Log[1-z]^3)/(3 z)+(16 (1-3 z+z^2-z^3) Log[1-z]^4)/(3 z)),
\[Gamma]24zqqgSpin :> 35/12+1/12 (223-8 \[Pi]^2) \[Epsilon]+(459/4-(35 \[Pi]^2)/9-24 Zeta[3]) \[Epsilon]^2+(2793/4-(223 \[Pi]^2)/9-(44 \[Pi]^4)/45-140 Zeta[3]) \[Epsilon]^3+(16875/4-(154 \[Pi]^4)/27-892 Zeta[3]+\[Pi]^2 (-153+32 Zeta[3])-672 Zeta[5]) \[Epsilon]^4 + \[Epsilon]^5 ORD[\[Gamma]24ZqqgSpin, 5],
\[Gamma]perptilde22Rggg :> -(1/3)-(4 \[Epsilon])/9+2/27 (-26+3 \[Pi]^2) \[Epsilon]^2+8/81 (-80+3 \[Pi]^2+54 Zeta[3]) \[Epsilon]^3+(2 (-19360+780 \[Pi]^2+81 \[Pi]^4+4320 Zeta[3]) \[Epsilon]^4)/1215+ep^5 ORD[\[Gamma]Perptilde22Rggg,5],

(*convolutions*)
Pqq22timesPqq02[z_] :>
      delta[1 - z]/\[Epsilon]^2 + 1/\[Epsilon] 2 (1 + z - 2 D0[z]) +
       16 D1[z] - (2 \[Pi]^2)/3 delta[
         1 - z] - (1 + z) (8 Log[1 - z] -
          Log[z]) + \[Epsilon] ((8 \[Pi]^2 D0[z])/3 -
          32 D2[z] - (1 + z) (2 \[Pi]^2 -
             16 Log[1 - z]^2 + (Log[z] - 2) Log[z] - 4 PolyLog[2, z]) -
          16 delta[1 - z] Zeta[3])+ \[Epsilon]^2 (-(32/3) \[Pi]^2 D1[z]+128/3 D3[z]-2/5 \[Pi]^4 delta[1-z]+1/3 (-30-4 \[Pi]^2+30 z-4 \[Pi]^2 z+24 \[Pi]^2 Log[1-z]+24 \[Pi]^2 z Log[1-z]-64 Log[1-z]^3-64 z Log[1-z]^3-15 Log[z]+2 \[Pi]^2 Log[z]-15 z Log[z]+2 \[Pi]^2 z Log[z]-24 Log[1-z]^2 Log[z]-24 z Log[1-z]^2 Log[z]-6 Log[z]^2-6 z Log[z]^2+2 Log[z]^3+2 z Log[z]^3-24 (1+z) (-1+2 Log[1-z]) PolyLog[2,z]-48 (1+z) PolyLog[3,1-z]-24 PolyLog[3,z]-24 z PolyLog[3,z])-24 (1+z) Zeta[3]+64 D0[z] Zeta[3]) + ORD[Pqq22timesPqq02,3] \[Epsilon]^3,

Pgq22timesPqg02[z_] :> (1+4/(3 z)-z-(4 z^2)/3+2 (1+z) Log[z])+1/(9 z) (46-3 (1+4 \[Pi]^2) z+(3-12 \[Pi]^2) z^2-46 z^3+12 (-1+z) (4+7 z+4 z^2) Log[1-z]-18 z (1+z) (Log[z]+Log[z]^2-4 PolyLog[2,z])) \[Epsilon]+1/(27 z) 2 (281-12 \[Pi]^2+21 z+9 \[Pi]^2 z-21 z^2+27 \[Pi]^2 z^2-281 z^3+12 \[Pi]^2 z^3+6 (-46+3 (1+4 \[Pi]^2) z+3 (-1+4 \[Pi]^2) z^2+46 z^3) Log[1-z]-36 (-1+z) (4+7 z+4 z^2) Log[1-z]^2+9 (3+2 \[Pi]^2) z (1+z) Log[z]-9 z (1+z) (24 Log[1-z]^2 Log[z]-3 Log[z]^2-2 Log[z]^3+12 (1+4 Log[1-z]) PolyLog[2,z]+48 PolyLog[3,1-z]+24 PolyLog[3,z])+216 z Zeta[3]+216 z^2 Zeta[3]) \[Epsilon]^2+ ORD[Pgq22TimesPqg02,3] \[Epsilon]^3, 
Pqg22timesPgq02[z_] :> (1+4/(3 z)-z-(4 z^2)/3+2 (1+z) Log[z])-((2+z+4 \[Pi]^2 z-z^2+4 \[Pi]^2 z^2-2 z^3-4 (-1+z) (4+7 z+4 z^2) Log[1-z]+8 (1+z)^3 Log[z]+6 z (1+z) (Log[z]^2-4 PolyLog[2,z])) \[Epsilon])/(3 z)+1/(9 z) 2 (9+4 \[Pi]^2+51 z+21 \[Pi]^2 z-51 z^2+27 \[Pi]^2 z^2-9 z^3+12 \[Pi]^2 z^3+6 (2+z+4 \[Pi]^2 z+(-1+4 \[Pi]^2) z^2-2 z^3) Log[1-z]-12 (-1+z) (4+7 z+4 z^2) Log[1-z]^2+6 (1+z) (1+(8+\[Pi]^2) z+z^2) Log[z]+12 (1+z)^3 (Log[z]^2-4 PolyLog[2,z])-6 z (1+z) (12 Log[1-z]^2 Log[z]-Log[z]^3+24 Log[1-z] PolyLog[2,z]+24 PolyLog[3,1-z]+12 PolyLog[3,z])+72 z Zeta[3]+72 z^2 Zeta[3]) \[Epsilon]^2+ ORD[Pqg22TimesPgq02,3] \[Epsilon]^3,
Pgg22timesPgg02[z_] :> (4 (-1+z) (11+2 z+11 z^2))/(3 z)+16 D1[z]+delta[1-z]/\[Epsilon]^2-(16 (-1+2 z-z^2+z^3) Log[1-z])/z-8 (1+z) Log[z]+((4 (-1+2 z-z^2+z^3))/z-4 D0[z]-2 delta[1-z] Log[z])/\[Epsilon]+delta[1-z] (-((2 \[Pi]^2)/3)+2 Log[z]^2)+\[Epsilon] (-((4 (67-6 \[Pi]^2-69 z+69 z^2-18 \[Pi]^2 z^2-67 z^3+6 \[Pi]^2 z^3))/(9 z))+8/3 \[Pi]^2 D0[z]-32 D2[z]-(16 (-1+z) (11+2 z+11 z^2) Log[1-z])/(3 z)+(32 (-1+2 z-z^2+z^3) Log[1-z]^2)/z+(4 (1+z) (11+4 z+11 z^2) Log[z])/(3 z)+8 (1+z) (Log[z]^2-4 PolyLog[2,z])+delta[1-z] (-(4/3) Log[z]^3-16 Zeta[3]))+\[Epsilon]^2 (-(32/3) \[Pi]^2 D1[z]+(128 D3[z])/3+(16 (67-6 \[Pi]^2-69 z+69 z^2-18 \[Pi]^2 z^2-67 z^3+6 \[Pi]^2 z^3) Log[1-z])/(9 z)+(32 (-1+z) (11+2 z+11 z^2) Log[1-z]^2)/(3 z)-(128 (-1+2 z-z^2+z^3) Log[1-z]^3)/(3 z)-16/3 (3+\[Pi]^2) (1+z) Log[z]+delta[1-z] (-((2 \[Pi]^4)/5)+(2 Log[z]^4)/3)-(4 (1+z) (11+4 z+11 z^2) (Log[z]^2-4 PolyLog[2,z]))/(3 z)+16/3 (1+z) (12 Log[1-z]^2 Log[z]-Log[z]^3+24 Log[1-z] PolyLog[2,z]+24 PolyLog[3,1-z]+12 PolyLog[3,z])+64 D0[z] Zeta[3]-(8 (413-339 z+72 \[Pi]^2 z+339 z^2+18 \[Pi]^2 z^2-413 z^3+66 \[Pi]^2 z^3-216 Zeta[3]+648 z Zeta[3]+216 z^3 Zeta[3]))/(27 z))+ ORD[Pgg22TimesPgg02,3] \[Epsilon]^3,
Pgq22timesPqq02[z_] :> -(-1+z) (-2+5 z)+(-1+2 z-2 z^2)/\[Epsilon]+4 (1-2 z+2 z^2) Log[1-z]+(1-2 z) Log[z]+\[Epsilon] (1/3 (-9+63 z-54 z^2+4 \[Pi]^2 z^2)+4 (-1+z) (-2+5 z) Log[1-z]-8 (1-2 z+2 z^2) Log[1-z]^2+(3+2 z) Log[z]+(-1+2 z) (Log[z]^2-4 PolyLog[2,z]))+\[Epsilon]^2 (-(4/3) (-9+63 z-54 z^2+4 \[Pi]^2 z^2) Log[1-z]-8 (-1+z) (-2+5 z) Log[1-z]^2+32/3 (1-2 z+2 z^2) Log[1-z]^3-2/3 (3-\[Pi]^2+3 z+2 \[Pi]^2 z) Log[z]+(-3-2 z) Log[z]^2+4 (3+2 z) PolyLog[2,z]+2/3 (-1+2 z) (12 Log[1-z]^2 Log[z]-Log[z]^3+24 Log[1-z] PolyLog[2,z]+24 PolyLog[3,1-z]+12 PolyLog[3,z])+2/3 (-24-\[Pi]^2+129 z-9 \[Pi]^2 z-105 z^2+5 \[Pi]^2 z^2+36 Zeta[3]-72 z Zeta[3]+48 z^2 Zeta[3]))+ ORD[Pgq22TimesPqq02,3] \[Epsilon]^3,
Pgg22timesPgq02[z_] :>  -(((-1+z) (4+7 z+37 z^2))/(3 z))+(-1+2 z-2 z^2)/\[Epsilon]+4 (1-2 z+2 z^2) Log[1-z]+2 (1+4 z) Log[z]+\[Epsilon] ((2 (-1+z-\[Pi]^2 z+38 z^2-10 \[Pi]^2 z^2-38 z^3+2 \[Pi]^2 z^3))/(3 z)+(4 (-1+z) (4+7 z+37 z^2) Log[1-z])/(3 z)-8 (1-2 z+2 z^2) Log[1-z]^2-(4 (2+6 z+9 z^2+11 z^3) Log[z])/(3 z)-2 (1+4 z) (Log[z]^2-4 PolyLog[2,z]))+\[Epsilon]^2 (-((8 (-1+z-\[Pi]^2 z+38 z^2-10 \[Pi]^2 z^2-38 z^3+2 \[Pi]^2 z^3) Log[1-z])/(3 z))-(8 (-1+z) (4+7 z+37 z^2) Log[1-z]^2)/(3 z)+32/3 (1-2 z+2 z^2) Log[1-z]^3+(4 (1+9 z+\[Pi]^2 z+3 z^2+4 \[Pi]^2 z^2-11 z^3) Log[z])/(3 z)+(4 (2+6 z+9 z^2+11 z^3) (Log[z]^2-4 PolyLog[2,z]))/(3 z)-4/3 (1+4 z) (12 Log[1-z]^2 Log[z]-Log[z]^3+24 Log[1-z] PolyLog[2,z]+24 PolyLog[3,1-z]+12 PolyLog[3,z])+(2 (9+4 \[Pi]^2+69 z+21 \[Pi]^2 z+282 z^2+6 \[Pi]^2 z^2-360 z^3+81 \[Pi]^2 z^3+144 z Zeta[3]+144 z^2 Zeta[3]+144 z^3 Zeta[3]))/(9 z))+ ORD[Pgg22TimesPgq02,3] \[Epsilon]^3,

PAP0timesPAP0[z_] :> 6 D0[z] + 8 D1[z] + (9/4 - (2 \[Pi]^2)/3) delta[1 - z] - ((3 z^2 + 1) Log[z])/(1 - z) - z - 4 (1 + z) Log[1 - z] - 5,

(* see Eq.D9-D12 in DY_2019_Other_Channels.pdf with convention gq \[TwoWayRule] gq *)

PAP0qqtimesPAP0qq[z_] :> 6 D0[z] + 8 D1[z] + (9/4 - (2 \[Pi]^2)/3) delta[1 - z] - ((3 z^2 + 1) Log[z])/(1 - z) - z - 4 (1 + z) Log[1 - z] - 5,

 PAP0gqtimesPAP0qq[z_]:> 2z -1/2 +2(1-2z+2z^2)*Log[1-z]-(1-2z+4z^2)*Log[z],   
 
 PAP0ggtimesPAP0gq[z_]:> \[Beta]0*(2z^2-2z+1)+CA(-31/3z^2 + (4z^2-4z+2)Log[1-z]+8z+4/3/z+(8z+2)Log[z]+1), 
 
 PAP0gqtimesPAP0qg[z_]:> 1 + 4/3/z - z - 4/3*z^2 + 2(1+z)Log[z] , 
 
 (* convolution computed with Federica's package *)
 PAP0ggtimesPAP0gg[z_]:> -((4 CA (-1+2 z-z^2+z^3) \[Beta]0)/z)+4 CA \[Beta]0 D0[z]+8 CA^2 D1[z]+(-(2/3) CA^2 \[Pi]^2+\[Beta]0^2) delta[1-z]+CA^2 ((4 (-1+z) (11+2 z+11 z^2))/(3 z)-(8 (-1+2 z-z^2+z^3) Log[1-z])/z+(4 (1+3 z^2-4 z^3+z^4) Log[z])/((-1+z) z)) , 
 PAP0qgtimesPAP0gg[z_]:> (\[Beta]0 CF (2-2 z+z^2))/z+CA CF (((-1+z) (31+7 z+4 z^2))/(3 z)+(2 (2-2 z+z^2) Log[1-z])/z-(4 (1+z+z^2) Log[z])/z) , 
(*one loop splitting functions*)
PqqAP1[z_] :> (* CF^2 *) ( 3 - 2*z - 2*(1 + z^2)/(1 - z)*Log[1 - z]*Log[z] +
          2*Log[z] + (1 + 3*z^2)/(2*(1 - z))*Log[z]^2 +
          2*(1 + z^2)/(1 - z)*PolyLog[2, 1 - z] +
          delta[1 - z]*(3/8 - \[Pi]^2/2 + 6*Zeta[3])) +
       0* CA CF ((3 \[Pi]^2 (1 + z) - 124 z - 19)/
          18 + (67/9 - \[Pi]^2/3) D0[z] + (2 + 11 z^2)/(6 (1 - z))
            Log[z] - (1 + z^2)/(1 - z) PolyLog[2, 1 - z] +
          delta[1 - z] (17/24 + 11/18 \[Pi]^2 - 3 Zeta[3])),
          
PgqAP1[z_] :> 1/2 CF TR (4 - 9*z - (1-4*z)*Log[z] - (1 - 2*z)*Log[z]^2 + 4*Log[1-z] + (2*Log[(1-z)/z]^2 - 4*Log[(1-z)/z] - 2/3*\[Pi]^2 +10)*(z^2 + (1-z)^2)) + 1/2*CA*TR*(182/9 + 14/9*z + 40/9/z + (136/3*z - 38/3)*Log[z] - 4*Log[1-z] - (2 + 8*z)*Log[z]^2 + 2*(z^2 + (1+z)^2)*(-2*PolyLog[2, -z] + 1/2*Log[z]^2 - 2*Log[z]*Log[1+z] - \[Pi]^2/6) + (-Log[z]^2 + 44/3*Log[z] - 2*Log[1-z]^2 + 4*Log[1-z] + \[Pi]^2/3 - 218/9 )*(z^2 + (1-z)^2 ) ), 
   
PAP0timesPqqfin[z_,E1_] :> -(1/\[Epsilon]) (((4 E1^2)/mu2)^-\[Epsilon]  (
          Gamma[1 - \[Epsilon]]^2 PqqNLOCV[z])/Gamma[1 - 2 \[Epsilon]] - Exp[-\[Epsilon] EulerGamma] Gamma[1 - \[Epsilon]] PAP0timesPAP0[z]),
   
(* Eq. C9 in 2009.10386 *) 
PgqNLOCV[z_]:> -3*z^2+2*(2*z^2-2*z+1)Log[1-z]+(-4z^2+2z-1)Log[z]+5z-2+\[Epsilon](2(1-2z)PolyLog[2,z]+1/3(\[Pi]^2(2z-1)-3(8z^2-9z+1))-2(2z^2-2z+1)Log[1-z]^2+(4z^2-2z+1)Log[z]^2+2(3z^2-5z+2)Log[1-z]+(2z+3)Log[z])+\[Epsilon]^2(4(2z-1)PolyLog[3,1-z]+4(2z-1)PolyLog[3,z]+PolyLog[2,z](2(2z+3)+4(2z-1)Log[1-z])+1/3(-12(4z^2+z(2 Zeta[3]-5)-Zeta[3]+1)-\[Pi]^2(2z+3))+4/3(2z^2-2z+1)Log[1-z]^3-2/3(4z^2-2z+1)Log[z]^3+(2(2z-1)Log[z]-2(3z^2-5z+2))Log[1-z]^2+2/3(24z^2-2\[Pi]^2z-27z+\[Pi]^2+3)Log[1-z]-(2z+3)Log[z]^2+2/3(\[Pi]^2(1-2z)-3)Log[z])+ \[Epsilon]^3 ORD[PgqCV,3],
              
(* Eq. C7 in 2009.10386 *)          
PAP0gqtimesPqqfin[z_,E1_] :> -(1/\[Epsilon]) (((4 E1^2)/mu2)^-\[Epsilon] Gamma[1 - \[Epsilon]]^2/Gamma[1 - 2 \[Epsilon]] (3/2 z^(-2 \[Epsilon])(z^2+(1-z)^2) + PgqNLOCV[z])- Exp[-\[Epsilon] EulerGamma] Gamma[1 - \[Epsilon]] PAP0gqtimesPAP0qq[z]),

(* convolutions Eq. C8 in 2009.10386 *)
PqqNLOIStimesPgqNLOIS[z_, E1_]:> -2+5z-3z^2+2(2z^2-2z+1)Log[1-z]-(4z^2-2z+1)Log[z]+\[Epsilon](4(1-2z)PolyLog[2,z]+1/3 (4(\[Pi]^2-12)z^2+57z-9)-6(2z^2-2z+1)Log[1-z]^2+(4z^2-2z+1)Log[z]^2+(8(2z^2-3z+1)+4(2z^2-2z+1)Log[z])Log[1-z]+(3+6z-4z^2)Log[z])+\[Epsilon]^2 (16(2z-1)PolyLog[3,1-z]+8(2z-1)PolyLog[3,z]+PolyLog[2,z](4(2z+3)+16(2z-1)Log[1-z])+8(3-6z+4z^2)Zeta[3]-(2\[Pi]^2)/3 (1+9z-5z^2)-16+84z-68z^2+28/3 (2z^2-2z+1)Log[1-z]^3-2/3 (4z^2-2z+1)Log[z]^3-4((9z^2-13z+4)+(2z^2-6z+3)Log[z])Log[1-z]^2+(-(4/3)((4\[Pi]^2-51)z^2+60z-9)-4(2z^2-2z+1)Log[z]^2-8z(1-z)Log[z])Log[1-z]+(4z^2-6z-3)Log[z]^2-2(1-z+2z^2-\[Pi]^2/3 (1-2z))Log[z])-1/\[Epsilon] ((1-2 z+2 z^2)+(-1+(1-2 z+2 z^2) (1-2 Log[1-z])) \[Epsilon]+(-1+2 Log[1-z]+(1-2 z+2 z^2) (1-2 Log[1-z]+2 Log[1-z]^2)) \[Epsilon]^2+(-1+2 Log[1-z]-2 Log[1-z]^2+(1-2 z+2 z^2) (1-2 Log[1-z]+2 Log[1-z]^2-4/3 Log[1-z]^3)) \[Epsilon]^3)(z^(-2\[Epsilon])-(Emax/E1)^(-2\[Epsilon])) + \[Epsilon]^3 ORD[PgqNLOISXPqqNLOIS,3],

PgqNLOIStimesPggNLOIS[z_, L1_]:> -((-4-3 z+6 L1 z-24 z^2-12 L1 z^2+31 z^3+12 L1 z^3)/(3 z))+2 (1-2 z+2 z^2) Log[1-z]+2 (1+4 z) Log[z]+ep ((2 (-1+z+3 L1^2 z-\[Pi]^2 z+35 z^2+6 L1 z^2-6 L1^2 z^2-10 \[Pi]^2 z^2-35 z^3-6 L1 z^3+6 L1^2 z^3+2 \[Pi]^2 z^3))/(3 z)+(4 (-4-3 z+3 L1 z-27 z^2-6 L1 z^2+34 z^3+6 L1 z^3) Log[1-z])/(3 z)-6 (1-2 z+2 z^2) Log[1-z]^2-(4 (2+6 z+9 z^2+11 z^3) Log[z])/(3 z)-2 (1+4 z) (Log[z]^2-4 PolyLog[2,z]))+ep^2 (-((4 (-2+2 z+3 L1^2 z-2 \[Pi]^2 z+73 z^2+6 L1 z^2-6 L1^2 z^2-20 \[Pi]^2 z^2-73 z^3-6 L1 z^3+6 L1^2 z^3+4 \[Pi]^2 z^3) Log[1-z])/(3 z))-(4 (-8-6 z+3 L1 z-57 z^2-6 L1 z^2+71 z^3+6 L1 z^3) Log[1-z]^2)/(3 z)+28/3 (1-2 z+2 z^2) Log[1-z]^3+(4 (1+9 z+\[Pi]^2 z+3 z^2+4 \[Pi]^2 z^2-11 z^3) Log[z])/(3 z)+(4 (2+6 z+9 z^2+11 z^3) (Log[z]^2-4 PolyLog[2,z]))/(3 z)-4/3 (1+4 z) (12 Log[1-z]^2 Log[z]-Log[z]^3+24 Log[1-z] PolyLog[2,z]+24 PolyLog[3,1-z]+12 PolyLog[3,z])+1/(9 z) 2 (9+4 \[Pi]^2+69 z-6 L1^3 z+21 \[Pi]^2 z+273 z^2+18 L1 z^2-18 L1^2 z^2+12 L1^3 z^2+6 \[Pi]^2 z^2-351 z^3-18 L1 z^3+18 L1^2 z^3-12 L1^3 z^3+81 \[Pi]^2 z^3+144 z Zeta[3]+144 z^2 Zeta[3]+144 z^3 Zeta[3]))+ ep^3 ORD[PgqNLOISXPggNLOIS, 3],   

(* \[CapitalGamma]2 functions for PDF renormalisations \[Rule] computed in ingredients_for_PDF_renormalisation *)
\[CapitalGamma]2gg[z_]:> -((\[Beta]0 (2 CA (-2+1/z+(1-z) z+D0[z])+\[Beta]0 delta[1-z]))/(2 \[Epsilon]^2))+1/(2 \[Epsilon]^2) (-((4 CA (-1+2 z-z^2+z^3) \[Beta]0)/z)+4 CA \[Beta]0 D0[z]+8 CA^2 D1[z]+(-(2/3) CA^2 \[Pi]^2+\[Beta]0^2) delta[1-z]+2 CF TR (1+4/(3 z)-z-(4 z^2)/3+2 (1+z) Log[z])+CA^2 ((4 (-1+z) (11+2 z+11 z^2))/(3 z)-(8 (-1+2 z-z^2+z^3) Log[1-z])/z+(4 (1+3 z^2-4 z^3+z^4) Log[z])/((-1+z) z)))+1/(2 \[Epsilon]) (CA Nf TR (2-2 z-20/9 (-2+1/(1-z)+1/z+(1-z) z)+26/9 (-(1/z)+z^2)-4/3 (1+z) Log[z])+CF Nf TR (-16+4/(3 z)+8 z+(20 z^2)/3-(6+10 z) Log[z]-(2+2 z) Log[z]^2)+CA^2 ((27 (1-z))/2+67/9 (-(1/z)+z^2)-(25/3-(11 z)/3+(44 z^2)/3) Log[z]+4 (1+z) Log[z]^2+(-2+1/(1-z)+1/z+(1-z) z) (67/9-\[Pi]^2/3-4 Log[1-z] Log[z]+Log[z]^2)+2 (-2-1/z+1/(1+z)-z (1+z)) (-(\[Pi]^2/6)+Log[z]^2/2-2 Log[z] Log[1+z]-2 PolyLog[2,-z]))+delta[1-z] (-(4/3) CA Nf TR-CF Nf TR+CA^2 (8/3+3 Zeta[3]))),
\[CapitalGamma]2gq[z_]:>-((TR ((1-z)^2+z^2) \[Beta]0)/(2 \[Epsilon]^2))+1/(2 \[Epsilon]^2) (CF TR (-(1/2)+2 z+2 (1-2 z+2 z^2) Log[1-z]-(1-2 z+4 z^2) Log[z])+TR ((1-2 z+2 z^2) \[Beta]0+CA (1+4/(3 z)+8 z-(31 z^2)/3+(2-4 z+4 z^2) Log[1-z]+(2+8 z) Log[z])))+1/(2 \[Epsilon]) (1/2 CF TR (4-9 z+4 Log[1-z]+((1-z)^2+z^2) (10-(2 \[Pi]^2)/3-4 Log[(1-z)/z]+2 Log[(1-z)/z]^2)-(1-4 z) Log[z]-(1-2 z) Log[z]^2)+1/2 CA TR (182/9+40/(9 z)+(14 z)/9-4 Log[1-z]+(-(38/3)+(136 z)/3) Log[z]-(2+8 z) Log[z]^2+((1-z)^2+z^2) (-(218/9)+\[Pi]^2/3+4 Log[1-z]-2 Log[1-z]^2+(44 Log[z])/3-Log[z]^2)+2 (z^2+(1+z)^2) (-(\[Pi]^2/6)+Log[z]^2/2-2 Log[z] Log[1+z]-2 PolyLog[2,-z]))),

\[CapitalGamma]2qq[z_]:>-((CF \[Beta]0 (-1-z+2 D0[z]+3/2 delta[1-z]))/(2 \[Epsilon]^2))+1/(2 \[Epsilon]^2) (CF TR (1+4/(3 z)-z-(4 z^2)/3+2 (1+z) Log[z])+CF^2 (-5-z+6 D0[z]+8 D1[z]+(9/4-(2 \[Pi]^2)/3) delta[1-z]-4 (1+z) Log[1-z]-((1+3 z^2) Log[z])/(1-z)))+1/(2 \[Epsilon]) (-CF Nf TR ((4 (1-z))/3-(10 (1+z))/9+(20 D0[z])/9+(1/6+(2 \[Pi]^2)/9) delta[1-z]+(2 (1+z^2) Log[z])/(3 (1-z)))+1/2 CF TR (-(16/3)+(4 (1-z))/3+40/(9 z)+(40 z)/3-(112 z^2)/9+(2+10 z+(16 z^2)/3) Log[z]-2 (1+z) Log[z]^2)+(CA-2 CF) CF (4-(7 z)/2+((5-2 z^2) Log[z])/(2 (1-z))+((1+z^2) Log[z]^2)/(2 (1-z))+((1+z^2) PolyLog[2,1-z])/(1-z))+CA CF (1/18 (-19-124 z+3 \[Pi]^2 (1+z))+(67/9-\[Pi]^2/3) D0[z]+((2+11 z^2) Log[z])/(6 (1-z))-((1+z^2) PolyLog[2,1-z])/(1-z)+delta[1-z] (17/24+(11 \[Pi]^2)/18-3 Zeta[3]))+CF^2 (3-2 z+2 Log[z]-(2 (1+z^2) Log[1-z] Log[z])/(1-z)+((1+3 z^2) Log[z]^2)/(2 (1-z))+(2 (1+z^2) PolyLog[2,1-z])/(1-z)+delta[1-z] (3/8-\[Pi]^2/2+6 Zeta[3]))),
\[CapitalGamma]2qqp[z_]:> (CF TR (1+4/(3 z)-z-(4 z^2)/3+2 (1+z) Log[z]))/(2 \[Epsilon]^2)+1/(4 \[Epsilon]) CF TR (-(16/3)+(4 (1-z))/3+40/(9 z)+(40 z)/3-(112 z^2)/9+(2+10 z+(16 z^2)/3) Log[z]-2 (1+z) Log[z]^2),
\[CapitalGamma]2qqb[z_]:> (CF TR (1+4/(3 z)-z-(4 z^2)/3+2 (1+z) Log[z]))/(2 \[Epsilon]^2)+1/(2 \[Epsilon]) (1/2 CF TR (-(16/3)+(4 (1-z))/3+40/(9 z)+(40 z)/3-(112 z^2)/9+(2+10 z+(16 z^2)/3) Log[z]-2 (1+z) Log[z]^2)+CF (-(CA/2)+CF) (4 (1-z)+2 (1+z) Log[z]+2 (-1+z+2/(1+z)) (-(\[Pi]^2/6)+Log[z]^2/2-2 Log[z] Log[1+z]-2 PolyLog[2,-z]))),
\[CapitalGamma]2qqbp[z_]:> (CF TR (1+4/(3 z)-z-(4 z^2)/3+2 (1+z) Log[z]))/(2 \[Epsilon]^2)+1/(4 \[Epsilon]) CF TR (-(16/3)+(4 (1-z))/3+40/(9 z)+(40 z)/3-(112 z^2)/9+(2+10 z+(16 z^2)/3) Log[z]-2 (1+z) Log[z]^2),

(* convolutions of AP splitting function and hard-collinear splitting [no PDF included, different with repsect to Pfin] \[Rule] see ingredients_PDF_renormalisation *)

PAP0qgtimesPgqNLO[z_,E1_] :> 1/9 ((6 CF (1+z) (1+(-3+\[Pi]^2) z+2 z^2))/z-6 CF \[Epsilon]-12 CF \[Pi]^2 \[Epsilon]-(18 CF \[Epsilon])/z-(8 CF \[Pi]^2 \[Epsilon])/z+12 CF z \[Epsilon]-6 CF \[Pi]^2 z \[Epsilon]+12 CF z^2 \[Epsilon]-(6 CF (1+z) (1+(-3+\[Pi]^2) z+2 z^2) \[Epsilon] Log[(4 E1^2)/mu2])/z-(6 CF (-1+z) (4+7 z+4 z^2) Log[1-z])/z-(12 CF (1+z) (1+(-3+\[Pi]^2) z+2 z^2) \[Epsilon] Log[1-z])/z+(6 CF (-1+z) (4+7 z+4 z^2) \[Epsilon] Log[(4 E1^2)/mu2] Log[1-z])/z+(6 CF (-1+z) (4+7 z+4 z^2) \[Epsilon] Log[1-z]^2)/z+(6 CF (4+6 z+3 z^2) Log[z])/z-(6 CF (2+2 (3+\[Pi]^2) z+(-3+2 \[Pi]^2) z^2) \[Epsilon] Log[z])/z-(6 CF (4+6 z+3 z^2) \[Epsilon] Log[(4 E1^2)/mu2] Log[z])/z+36 CF (1+z) \[Epsilon] Log[1-z]^2 Log[z]+18 CF (1+z) Log[z]^2-(6 CF (4+6 z+3 z^2) \[Epsilon] Log[z]^2)/z-18 CF (1+z) \[Epsilon] Log[(4 E1^2)/mu2] Log[z]^2-12 CF (1+z) \[Epsilon] Log[z]^3+(3 CF (-4-3 z+3 z^2+4 z^3-6 z (1+z) Log[z]))/(z \[Epsilon])-(CF \[Epsilon] (\[Pi]^2-3 Log[4]^2-6 Log[4] Log[E1^2/mu2]-3 Log[E1^2/mu2]^2) (-4-3 z+3 z^2+4 z^3-6 z (1+z) Log[z]))/(2 z)-(3 CF Log[(4 E1^2)/mu2] (-4-3 z+3 z^2+4 z^3-6 z (1+z) Log[z]))/z+9 \[Epsilon]^3 ORD[CONV,3]+12 CF (-3 (1+z)+(6+4/z+3 z) \[Epsilon]+3 (1+z) \[Epsilon] Log[(4 E1^2)/mu2]+6 (1+z) \[Epsilon] Log[1-z]) PolyLog[2,z]+72 CF (1+z) \[Epsilon] PolyLog[3,1-z]+72 CF (1+z) \[Epsilon] PolyLog[3,z]-72 CF (1+z) \[Epsilon] Zeta[3]),
PAP0qgtimesPgqNLO[z_,E1_] :> 1/2 CF (3-3 z+2 \[Pi]^2 (1+z))+2 CF \[Epsilon]+5/3 CF \[Pi]^2 \[Epsilon]-2 CF z \[Epsilon]-CF \[Pi]^2 z \[Epsilon]-4/3 CF \[Pi]^2 D0[z]+6 CF D1[z]+8/3 CF \[Pi]^2 \[Epsilon] D1[z]+12 CF D2[z]-6 CF \[Epsilon] D2[z]-32/3 CF \[Epsilon] D3[z]-1/2 CF (3-3 z+2 \[Pi]^2 (1+z)) \[Epsilon] Log[(4 E1^2)/mu2]+4/3 CF \[Pi]^2 \[Epsilon] D0[z] Log[(4 E1^2)/mu2]-6 CF \[Epsilon] D1[z] Log[(4 E1^2)/mu2]-12 CF \[Epsilon] D2[z] Log[(4 E1^2)/mu2]-CF (5+z) Log[1-z]-CF (3-3 z+2 \[Pi]^2 (1+z)) \[Epsilon] Log[1-z]+CF (5+z) \[Epsilon] Log[(4 E1^2)/mu2] Log[1-z]-6 CF (1+z) Log[1-z]^2+3 CF (1+z) \[Epsilon] Log[1-z]^2+6 CF (1+z) \[Epsilon] Log[(4 E1^2)/mu2] Log[1-z]^2+16/3 CF (1+z) \[Epsilon] Log[1-z]^3-(2 CF z (2+z) Log[z])/(-1+z)-2/3 CF (-3+\[Pi]^2 (1+z)) \[Epsilon] Log[z]+9/2 CF delta[1-z] Log[z]+(2 CF z (2+z) \[Epsilon] Log[(4 E1^2)/mu2] Log[z])/(-1+z)-9/2 CF \[Epsilon] delta[1-z] Log[(4 E1^2)/mu2] Log[z]+2 CF (1+z) \[Epsilon] Log[1-z]^2 Log[z]+(CF (1+3 z^2) Log[z]^2)/(-1+z)+(2 CF z (2+z) \[Epsilon] Log[z]^2)/(-1+z)-(CF (1+3 z^2) \[Epsilon] Log[(4 E1^2)/mu2] Log[z]^2)/(-1+z)-(2 CF (1+3 z^2) \[Epsilon] Log[z]^3)/(3 (-1+z))+1/\[Epsilon] (5 CF+CF z-6 CF D0[z]-8 CF D1[z]+1/12 CF (-27+8 \[Pi]^2) delta[1-z]+4 CF (1+z) Log[1-z]-(CF (1+3 z^2) Log[z])/(-1+z))+1/90 CF \[Epsilon] delta[1-z] (16 \[Pi]^4-405 Log[z]^2)-1/6 \[Epsilon] (-\[Pi]^2+3 Log[4]^2+6 Log[4] Log[E1^2/mu2]+3 Log[E1^2/mu2]^2) (6 CF D0[z]+8 CF D1[z]+1/12 CF (27-8 \[Pi]^2) delta[1-z]+CF (-5-z-4 (1+z) Log[1-z]+((1+3 z^2) Log[z])/(-1+z)))+Log[(4 E1^2)/mu2] (6 CF D0[z]+8 CF D1[z]+1/12 CF (27-8 \[Pi]^2) delta[1-z]+CF (-5-z-4 (1+z) Log[1-z]+((1+3 z^2) Log[z])/(-1+z)))+2 CF (-1-z+(-3+z) \[Epsilon]+(1+z) \[Epsilon] Log[(4 E1^2)/mu2]+2 (1+z) \[Epsilon] Log[1-z]) PolyLog[2,z]+4 CF (1+z) \[Epsilon] PolyLog[3,1-z]+4 CF (1+z) \[Epsilon] PolyLog[3,z]+4 CF ((1+z) \[Epsilon]-4 \[Epsilon] D0[z]+delta[1-z] (2-2 \[Epsilon] Log[(4 E1^2)/mu2])) Zeta[3] + \[Epsilon]^2 ORD[PAP0qqTimesPqqNLO,2],
PAP0ggtimesPggNLO[z_,E1_] :> (2 CA^2 (67-69 z+69 z^2-67 z^3+6 \[Pi]^2 (-1-3 z^2+z^3)))/(9 z)-(274 CA^2 \[Epsilon])/9+4/3 CA^2 \[Pi]^2 \[Epsilon]+(898 CA^2 \[Epsilon])/(27 z)+(44 CA^2 \[Pi]^2 \[Epsilon])/(9 z)+250/9 CA^2 z \[Epsilon]+16/3 CA^2 \[Pi]^2 z \[Epsilon]-826/27 CA^2 z^2 \[Epsilon]-4/3 CA^2 \[Pi]^2 D0[z]+4 CA \[Beta]0 D1[z]+8/3 CA^2 \[Pi]^2 \[Epsilon] D1[z]+12 CA^2 D2[z]-4 CA \[Beta]0 \[Epsilon] D2[z]-32/3 CA^2 \[Epsilon] D3[z]-(2 CA^2 (67-69 z+69 z^2-67 z^3+6 \[Pi]^2 (-1-3 z^2+z^3)) \[Epsilon] Log[(4 E1^2)/mu2])/(9 z)+4/3 CA^2 \[Pi]^2 \[Epsilon] D0[z] Log[(4 E1^2)/mu2]-4 CA \[Beta]0 \[Epsilon] D1[z] Log[(4 E1^2)/mu2]-12 CA^2 \[Epsilon] D2[z] Log[(4 E1^2)/mu2]+(8 CA^2 (-1+z) (11+2 z+11 z^2) Log[1-z])/(3 z)-(4 CA (-1+2 z-z^2+z^3) \[Beta]0 Log[1-z])/z-(4 CA^2 (67-69 z+69 z^2-67 z^3+6 \[Pi]^2 (-1-3 z^2+z^3)) \[Epsilon] Log[1-z])/(9 z)-(8 CA^2 (-1+z) (11+2 z+11 z^2) \[Epsilon] Log[(4 E1^2)/mu2] Log[1-z])/(3 z)+(4 CA (-1+2 z-z^2+z^3) \[Beta]0 \[Epsilon] Log[(4 E1^2)/mu2] Log[1-z])/z-(12 CA^2 (-1+2 z-z^2+z^3) Log[1-z]^2)/z-(8 CA^2 (-1+z) (11+2 z+11 z^2) \[Epsilon] Log[1-z]^2)/(3 z)+(4 CA (-1+2 z-z^2+z^3) \[Beta]0 \[Epsilon] Log[1-z]^2)/z+(12 CA^2 (-1+2 z-z^2+z^3) \[Epsilon] Log[(4 E1^2)/mu2] Log[1-z]^2)/z+(32 CA^2 (-1+2 z-z^2+z^3) \[Epsilon] Log[1-z]^3)/(3 z)-(4 CA^2 (11+3 z+12 z^2) Log[z])/(3 z)-(4 CA (1-z+z^2)^2 \[Beta]0 Log[z])/((-1+z) z)+8/3 CA^2 (3+2 \[Pi]^2 (1+z)) \[Epsilon] Log[z]+2 \[Beta]0^2 delta[1-z] Log[z]+(4 CA^2 (11+3 z+12 z^2) \[Epsilon] Log[(4 E1^2)/mu2] Log[z])/(3 z)+(4 CA (1-z+z^2)^2 \[Beta]0 \[Epsilon] Log[(4 E1^2)/mu2] Log[z])/((-1+z) z)-2 \[Beta]0^2 \[Epsilon] delta[1-z] Log[(4 E1^2)/mu2] Log[z]-16 CA^2 (1+z) \[Epsilon] Log[1-z]^2 Log[z]+(4 CA^2 (1+3 z^2-4 z^3+z^4) Log[z]^2)/((-1+z) z)+(4 CA^2 (11+3 z+12 z^2) \[Epsilon] Log[z]^2)/(3 z)+(4 CA (1-z+z^2)^2 \[Beta]0 \[Epsilon] Log[z]^2)/((-1+z) z)-(4 CA^2 (1+3 z^2-4 z^3+z^4) \[Epsilon] Log[(4 E1^2)/mu2] Log[z]^2)/((-1+z) z)-(8 CA^2 (1+3 z^2-4 z^3+z^4) \[Epsilon] Log[z]^3)/(3 (-1+z) z)+2/45 \[Epsilon] delta[1-z] (4 CA^2 \[Pi]^4-45 \[Beta]0^2 Log[z]^2)+1/\[Epsilon] ((4 CA (-1+2 z-z^2+z^3) \[Beta]0)/z-4 CA \[Beta]0 D0[z]-8 CA^2 D1[z]-(-(2/3) CA^2 \[Pi]^2+\[Beta]0^2) delta[1-z]+1/(3 (-1+z) z) 4 CA^2 (-(-1+z)^2 (11+2 z+11 z^2)+6 (1-3 z+3 z^2-2 z^3+z^4) Log[1-z]-3 (1+3 z^2-4 z^3+z^4) Log[z]))-1/6 \[Epsilon] (-\[Pi]^2+3 Log[4]^2+6 Log[4] Log[E1^2/mu2]+3 Log[E1^2/mu2]^2) (-((4 CA (-1+2 z-z^2+z^3) \[Beta]0)/z)+4 CA \[Beta]0 D0[z]+8 CA^2 D1[z]+(-(2/3) CA^2 \[Pi]^2+\[Beta]0^2) delta[1-z]+(4 CA^2 ((-1+z) (11+2 z+11 z^2)-6 (-1+2 z-z^2+z^3) Log[1-z]+(3 (1+3 z^2-4 z^3+z^4) Log[z])/(-1+z)))/(3 z))+Log[(4 E1^2)/mu2] (-((4 CA (-1+2 z-z^2+z^3) \[Beta]0)/z)+4 CA \[Beta]0 D0[z]+8 CA^2 D1[z]+(-(2/3) CA^2 \[Pi]^2+\[Beta]0^2) delta[1-z]+(4 CA^2 ((-1+z) (11+2 z+11 z^2)-6 (-1+2 z-z^2+z^3) Log[1-z]+(3 (1+3 z^2-4 z^3+z^4) Log[z])/(-1+z)))/(3 z))+8/3 CA^2 (6+6 z-3 \[Epsilon]-(11 \[Epsilon])/z-12 z \[Epsilon]-6 (1+z) \[Epsilon] Log[(4 E1^2)/mu2]-12 (1+z) \[Epsilon] Log[1-z]) PolyLog[2,z]-32 CA^2 (1+z) \[Epsilon] PolyLog[3,1-z]-32 CA^2 (1+z) \[Epsilon] PolyLog[3,z]+1/z 8 CA^2 (2 (-1+4 z+z^2+z^3) \[Epsilon]-2 z \[Epsilon] D0[z]-z delta[1-z] (-1+\[Epsilon] Log[4]+\[Epsilon] Log[E1^2/mu2])) Zeta[3] + \[Epsilon]^2 ORD[PAP0ggTimesPggNLO,2],
PAP0ggtimesPgqNLO[z_,E1_] :> 1/9 (-((3 CA (-2+7 z+(44-12 \[Pi]^2) z^2+(-49+4 \[Pi]^2) z^3))/z)-18 (-1+z) z \[Beta]0-33 CA \[Epsilon]-12 CA \[Pi]^2 \[Epsilon]-(18 CA \[Epsilon])/z-(8 CA \[Pi]^2 \[Epsilon])/z-204 CA z \[Epsilon]-12 CA \[Pi]^2 z \[Epsilon]+255 CA z^2 \[Epsilon]-12 CA \[Pi]^2 z^2 \[Epsilon]+(3 CA (-2+7 z+(44-12 \[Pi]^2) z^2+(-49+4 \[Pi]^2) z^3) \[Epsilon] Log[(4 E1^2)/mu2])/z+18 (-1+z) z \[Beta]0 \[Epsilon] Log[(4 E1^2)/mu2]-(6 CA (-1+z) (4+7 z+37 z^2) Log[1-z])/z+18 (1-2 z+2 z^2) \[Beta]0 Log[1-z]+(6 CA (-2+7 z+(50-12 \[Pi]^2) z^2+(-55+4 \[Pi]^2) z^3) \[Epsilon] Log[1-z])/z+(6 CA (-1+z) (4+7 z+37 z^2) \[Epsilon] Log[(4 E1^2)/mu2] Log[1-z])/z-18 (1-2 z+2 z^2) \[Beta]0 \[Epsilon] Log[(4 E1^2)/mu2] Log[1-z]+36 CA (1-2 z+2 z^2) Log[1-z]^2+(6 CA (-1+z) (4+7 z+43 z^2) \[Epsilon] Log[1-z]^2)/z-18 (1-2 z+2 z^2) \[Beta]0 \[Epsilon] Log[1-z]^2-36 CA (1-2 z+2 z^2) \[Epsilon] Log[(4 E1^2)/mu2] Log[1-z]^2-36 CA (1-2 z+2 z^2) \[Epsilon] Log[1-z]^3+18 (-1+z) z \[Beta]0 \[Epsilon] (-1+2 Log[1-z])+(12 CA (2+3 z+6 z^2) Log[z])/z-(12 CA (1+(3+\[Pi]^2) z+(-6+4 \[Pi]^2) z^2) \[Epsilon] Log[z])/z-(12 CA (2+3 z+6 z^2) \[Epsilon] Log[(4 E1^2)/mu2] Log[z])/z+36 CA (1+4 z) \[Epsilon] Log[1-z]^2 Log[z]+18 CA (1+4 z) Log[z]^2-(12 CA (2+3 z+6 z^2) \[Epsilon] Log[z]^2)/z-18 CA (1+4 z) \[Epsilon] Log[(4 E1^2)/mu2] Log[z]^2-12 CA (1+4 z) \[Epsilon] Log[z]^3+(9 ((CA (-1+z) (4+7 z+31 z^2))/(3 z)-(1-2 z+2 z^2) \[Beta]0-2 CA (1-2 z+2 z^2) Log[1-z]-2 CA (1+4 z) Log[z]))/\[Epsilon]+1/(2 z) \[Epsilon] (\[Pi]^2-3 Log[4]^2-6 Log[4] Log[E1^2/mu2]-3 Log[E1^2/mu2]^2) (4 CA+3 CA z+24 CA z^2-31 CA z^3+3 z \[Beta]0-6 z^2 \[Beta]0+6 z^3 \[Beta]0+6 CA z (1-2 z+2 z^2) Log[1-z]+6 CA z (1+4 z) Log[z])+9 Log[(4 E1^2)/mu2] ((1-2 z+2 z^2) \[Beta]0+CA (1+4/(3 z)+8 z-(31 z^2)/3+(2-4 z+4 z^2) Log[1-z]+(2+8 z) Log[z]))+12 CA (-3-12 z+6 \[Epsilon]+(4 \[Epsilon])/z+12 z \[Epsilon]+3 (\[Epsilon]+4 z \[Epsilon]) Log[(4 E1^2)/mu2]+6 (\[Epsilon]+4 z \[Epsilon]) Log[1-z]) PolyLog[2,z]+72 CA (1+4 z) \[Epsilon] PolyLog[3,1-z]+72 CA (1+4 z) \[Epsilon] PolyLog[3,z]-144 CA (1+z+z^2) \[Epsilon] Zeta[3])+\[Epsilon]^2 ORD[PAP0ggTimesPgqNLO,2],
PAP0gqtimesPqqNLO[z_,E1_] :> 1/12 TR (4 (3+\[Pi]^2-27 z-2 \[Pi]^2 z+24 z^2)+48 \[Epsilon]+12 \[Pi]^2 \[Epsilon]-240 z \[Epsilon]+8 \[Pi]^2 z \[Epsilon]+192 z^2 \[Epsilon]+4 (\[Pi]^2 (-1+2 z)-3 (1-9 z+8 z^2)) \[Epsilon] Log[(4 E1^2)/mu2]-24 (-1+z) (-2+3 z) Log[1-z]+8 (\[Pi]^2 (-1+2 z)-3 (1-9 z+8 z^2)) \[Epsilon] Log[1-z]+24 (-1+z) (-2+3 z) \[Epsilon] Log[(4 E1^2)/mu2] Log[1-z]+24 (1-2 z+2 z^2) Log[1-z]^2+24 (-1+z) (-2+3 z) \[Epsilon] Log[1-z]^2-24 (1-2 z+2 z^2) \[Epsilon] Log[(4 E1^2)/mu2] Log[1-z]^2-16 (1-2 z+2 z^2) \[Epsilon] Log[1-z]^3+24 z (-4+3 z) Log[z]+8 (3+\[Pi]^2 (-1+2 z)) \[Epsilon] Log[z]-24 z (-4+3 z) \[Epsilon] Log[(4 E1^2)/mu2] Log[z]-24 (-1+2 z) \[Epsilon] Log[1-z]^2 Log[z]-12 (1-2 z+4 z^2) Log[z]^2-24 z (-4+3 z) \[Epsilon] Log[z]^2+12 (1-2 z+4 z^2) \[Epsilon] Log[(4 E1^2)/mu2] Log[z]^2+8 (1-2 z+4 z^2) \[Epsilon] Log[z]^3-(6 (-1+4 z+(4-8 z+8 z^2) Log[1-z]+(-2+4 z-8 z^2) Log[z]))/\[Epsilon]+\[Epsilon] (\[Pi]^2-3 Log[4]^2-6 Log[4] Log[E1^2/mu2]-3 Log[E1^2/mu2]^2) (-1+4 z+(4-8 z+8 z^2) Log[1-z]+(-2+4 z-8 z^2) Log[z])+6 Log[(4 E1^2)/mu2] (-1+4 z+(4-8 z+8 z^2) Log[1-z]+(-2+4 z-8 z^2) Log[z])+24 (-1+2 z-3 \[Epsilon]-2 z \[Epsilon]+(\[Epsilon]-2 z \[Epsilon]) Log[(4 E1^2)/mu2]+(2-4 z) \[Epsilon] Log[1-z]) PolyLog[2,z]-48 (-1+2 z) \[Epsilon] PolyLog[3,1-z]-48 (-1+2 z) \[Epsilon] PolyLog[3,z]+48 (-1+2 z) \[Epsilon] Zeta[3])+\[Epsilon]^2 ORD[PAP0gqTimesPqqNLO,2],
PAP0gqtimesPqgNLO[z_,E1_] :> (TR (-26+3 (5+2 \[Pi]^2) z+(3+6 \[Pi]^2) z^2+8 z^3))/(9 z)+(34 TR \[Epsilon])/9-2/3 \[Pi]^2 TR \[Epsilon]-(178 TR \[Epsilon])/(27 z)+(2 TR z \[Epsilon])/9-2/3 \[Pi]^2 TR z \[Epsilon]+70/27 TR z^2 \[Epsilon]-(TR (-26+3 (5+2 \[Pi]^2) z+(3+6 \[Pi]^2) z^2+8 z^3) \[Epsilon] Log[(4 E1^2)/mu2])/(9 z)-(2 TR (-1+z) (4+7 z+4 z^2) Log[1-z])/(3 z)-(2 TR (-26+3 (5+2 \[Pi]^2) z+(3+6 \[Pi]^2) z^2+8 z^3) \[Epsilon] Log[1-z])/(9 z)+(2 TR (-1+z) (4+7 z+4 z^2) \[Epsilon] Log[(4 E1^2)/mu2] Log[1-z])/(3 z)+(2 TR (-1+z) (4+7 z+4 z^2) \[Epsilon] Log[1-z]^2)/(3 z)+2 TR (1+z) Log[z]-2/3 TR (3+2 \[Pi]^2 (1+z)) \[Epsilon] Log[z]-2 TR (1+z) \[Epsilon] Log[(4 E1^2)/mu2] Log[z]+4 TR (1+z) \[Epsilon] Log[1-z]^2 Log[z]+2 TR (1+z) Log[z]^2-2 TR (1+z) \[Epsilon] Log[z]^2-2 TR (1+z) \[Epsilon] Log[(4 E1^2)/mu2] Log[z]^2-4/3 TR (1+z) \[Epsilon] Log[z]^3+(TR (-3-4/z+3 z+4 z^2-6 (1+z) Log[z]))/(3 \[Epsilon])+1/3 TR Log[(4 E1^2)/mu2] (3+4/z-3 z-4 z^2+6 (1+z) Log[z])-(TR \[Epsilon] (\[Pi]^2-3 Log[4]^2-6 Log[4] Log[E1^2/mu2]-3 Log[E1^2/mu2]^2) (-4-3 z+3 z^2+4 z^3-6 z (1+z) Log[z]))/(18 z)+\[Epsilon]^3 ORD[CONV,3]+4 TR (1+z) (-1+\[Epsilon]+\[Epsilon] Log[(4 E1^2)/mu2]+2 \[Epsilon] Log[1-z]) PolyLog[2,z]+8 TR (1+z) \[Epsilon] PolyLog[3,1-z]+8 TR (1+z) \[Epsilon] PolyLog[3,z]-8 TR (1+z) \[Epsilon] Zeta[3]+\[Epsilon]^2 ORD[PAP0gqTimesPqgNLO,2], 
PAP0qqtimesPqqNLO[z_,E1_] :> 1/2 CF (3-3 z+2 \[Pi]^2 (1+z))+2 CF \[Epsilon]+5/3 CF \[Pi]^2 \[Epsilon]-2 CF z \[Epsilon]-CF \[Pi]^2 z \[Epsilon]-4/3 CF \[Pi]^2 D0[z]+6 CF D1[z]+8/3 CF \[Pi]^2 \[Epsilon] D1[z]+12 CF D2[z]-6 CF \[Epsilon] D2[z]-32/3 CF \[Epsilon] D3[z]-1/2 CF (3-3 z+2 \[Pi]^2 (1+z)) \[Epsilon] Log[(4 E1^2)/mu2]+4/3 CF \[Pi]^2 \[Epsilon] D0[z] Log[(4 E1^2)/mu2]-6 CF \[Epsilon] D1[z] Log[(4 E1^2)/mu2]-12 CF \[Epsilon] D2[z] Log[(4 E1^2)/mu2]-CF (5+z) Log[1-z]-CF (3-3 z+2 \[Pi]^2 (1+z)) \[Epsilon] Log[1-z]+CF (5+z) \[Epsilon] Log[(4 E1^2)/mu2] Log[1-z]-6 CF (1+z) Log[1-z]^2+3 CF (1+z) \[Epsilon] Log[1-z]^2+6 CF (1+z) \[Epsilon] Log[(4 E1^2)/mu2] Log[1-z]^2+16/3 CF (1+z) \[Epsilon] Log[1-z]^3-(2 CF z (2+z) Log[z])/(-1+z)-2/3 CF (-3+\[Pi]^2 (1+z)) \[Epsilon] Log[z]+9/2 CF delta[1-z] Log[z]+(2 CF z (2+z) \[Epsilon] Log[(4 E1^2)/mu2] Log[z])/(-1+z)-9/2 CF \[Epsilon] delta[1-z] Log[(4 E1^2)/mu2] Log[z]+2 CF (1+z) \[Epsilon] Log[1-z]^2 Log[z]+(CF (1+3 z^2) Log[z]^2)/(-1+z)+(2 CF z (2+z) \[Epsilon] Log[z]^2)/(-1+z)-(CF (1+3 z^2) \[Epsilon] Log[(4 E1^2)/mu2] Log[z]^2)/(-1+z)-(2 CF (1+3 z^2) \[Epsilon] Log[z]^3)/(3 (-1+z))+1/\[Epsilon] (5 CF+CF z-6 CF D0[z]-8 CF D1[z]+1/12 CF (-27+8 \[Pi]^2) delta[1-z]+4 CF (1+z) Log[1-z]-(CF (1+3 z^2) Log[z])/(-1+z))+1/90 CF \[Epsilon] delta[1-z] (16 \[Pi]^4-405 Log[z]^2)-1/6 \[Epsilon] (-\[Pi]^2+3 Log[4]^2+6 Log[4] Log[E1^2/mu2]+3 Log[E1^2/mu2]^2) (6 CF D0[z]+8 CF D1[z]+1/12 CF (27-8 \[Pi]^2) delta[1-z]+CF (-5-z-4 (1+z) Log[1-z]+((1+3 z^2) Log[z])/(-1+z)))+Log[(4 E1^2)/mu2] (6 CF D0[z]+8 CF D1[z]+1/12 CF (27-8 \[Pi]^2) delta[1-z]+CF (-5-z-4 (1+z) Log[1-z]+((1+3 z^2) Log[z])/(-1+z)))+2 CF (-1-z+(-3+z) \[Epsilon]+(1+z) \[Epsilon] Log[(4 E1^2)/mu2]+2 (1+z) \[Epsilon] Log[1-z]) PolyLog[2,z]+4 CF (1+z) \[Epsilon] PolyLog[3,1-z]+4 CF (1+z) \[Epsilon] PolyLog[3,z]+4 CF ((1+z) \[Epsilon]-4 \[Epsilon] D0[z]+delta[1-z] (2-2 \[Epsilon] Log[(4 E1^2)/mu2])) Zeta[3] +\[Epsilon]^2 ORD[PAP0gqTimesPqgNLO,2]

} /. ep -> \[Epsilon];
