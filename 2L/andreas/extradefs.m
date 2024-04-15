fromsqrt = {sqrt1[s_, m_] :> Sqrt[s*(s - 4*m^2)], 
   sqrtren1[m1_, m2_] :> Sqrt[4*m1^2 - m2^2], 
   sqrtren2[s_, m1_, m2_] :> 
    Sqrt[m2^4 - 2*m2^2*m1^2 + m1^4 - 2*m2^2*s - 2*m1^2*s + s^2], 
   sqrt2[s_, t_, m_] :> Sqrt[s*t*(s*t - 4*m^2*(t + m^2))], 
   sqrt3[s_, t_, m_] :> 
    Sqrt[s*(m^2*s*(m^2 - 2*t) + (-4*m^2 + s)*t^2)]};
abbrroots = {
   sqrt1[s, mw] -> sqrt[1],
   sqrt1[s, mz] -> sqrt[2],
   sqrt2[s, t, mz] -> sqrt[3],
   sqrt2[s, u, mw] -> sqrt[4],
   sqrt2[s, u, mz] -> sqrt[5],
   sqrt3[s, t, mz] -> sqrt[6],
   sqrt3[s, u, mw] -> sqrt[7],
   sqrt3[s, u, mz] -> sqrt[8],
   sqrtren1[mw, mz] -> sqrt[9]/mz,
   sqrtren1[mw, mh] -> sqrt[10]/mh,
   sqrtren1[mz, mh] -> sqrt[11]/mh,
   sqrtren2[s, mz, mh] -> sqrt[12]};
noexplicittrans = {iPi -> mi[33, 1, 1], zeta[2] -> mi[33, 1, 2], 
   zeta[3] -> mi[33, 1, 3], zeta[4] -> mi[33, 1, 4]};
noexplicitlogs = { L[3430] -> mi[32, 1, 1] + L[3431] };
(* mi[32,1,1]:=Log[t/u]=Log[-x]-Log[1+x]=L[3430]-L[3431] *) 

