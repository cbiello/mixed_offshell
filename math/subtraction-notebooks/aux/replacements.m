(* ::Package:: *)

(*Useful replacements *)
 
coupl = {gsb->( 8\[Pi]^2  Gamma[1-\[Epsilon]])/(4\[Pi])^\[Epsilon] \[Alpha]sB \[Mu]^(2\[Epsilon]), eb->( 8\[Pi]^2  Gamma[1-\[Epsilon]])/(4\[Pi])^\[Epsilon] \[Alpha]B \[Mu]^(2\[Epsilon])};
bracoupl={\[Alpha]sB-> \[Alpha]sontwopi E^(\[Epsilon] EulerGamma)/Gamma[1-\[Epsilon]], \[Alpha]B-> \[Alpha]ontwopi E^(\[Epsilon] EulerGamma)/Gamma[1-\[Epsilon]]};
Emaxchoice = {Emax-> \[Sqrt]s/2};
pcm = {E1-> \[Sqrt]s/2,E2-> \[Sqrt]s/2,E3-> \[Sqrt]s/2, E4-> \[Sqrt]s/2};
simpl= {E1-> EC,E2-> EC,Emax-> EC, eta[1,2]->1};

ListSub = {PolyLog[2, 1 - z_] :>  -PolyLog[2, z] + Pi^2/6 -
         Log[z] Log[1 - z] };

replaceSij = {s12 -> 4 E1 E2 eta[1, 2], s13 -> 4 E1 E3 eta[1, 3],
    s23 -> 4 E2 E3 eta[2, 3], s24 -> 4 E2 E4 eta[2, 4], s14 -> 4 E1 E4 eta[1, 4], s34 -> 4 E3 E4 eta[3, 4]};

LogSubs = {L1 -> 1/2 Log[Emax^2/E1^2], L2 -> 1/2 Log[Emax^2/E2^2],
    L3 -> 1/2 Log[Emax^2/E3^2],  L4 -> 1/2 Log[Emax^2/E4^2]};
 
SpecialLogSubs = {
    FF[Log[Emax^2/c_^2] - Log[Emax^2/(a_ c_)] - Log[Emax^2/(b_ c_)] +
       x___] :> FF[-Log[Emax^2/(a b)] + x],
    FF[-Log[Emax^2/c_^2] + Log[Emax^2/(a_ c_)] + Log[Emax^2/(b_ c_)] +
       x___] :> FF[Log[Emax^2/(a b)] + x],
    FF[- Log[Emax^2/a_^2] - Log[Emax^2/b_^2] + x___] :>
     FF[-2 Log[Emax^2/(a b)] + x],
    FF[ Log[Emax^2/a_^2] + Log[Emax^2/b_^2] + x___] :>
     FF[2 Log[Emax^2/(a b)] + x],
    FF[-s_ Log[Emax^2/a_^2] - s_ Log[Emax^2/b_^2] + x___] :>
     FF[-2 s Log[Emax^2/(a b)] + x],
    FF[s_ Log[Emax^2/a_^2] + s_ Log[Emax^2/b_^2] + x___] :>
     FF[2 s Log[Emax^2/(a b)] + x]
    };
