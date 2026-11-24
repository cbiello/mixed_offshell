(* usage: NBDIR=<math/subtraction-notebooks> NBOUT=<output dir> wolframscript -file run_W_NNLO_polecheck.wls   -> NBOUT/wpole_finite.m, about 4 min
   Evaluates W_NNLO_polecheck_cleaned.nb headless (gq, photon-quark, photon-gluon, four-quark chapters; the q qb' part, cells 12-60,
   is skipped: it is covered by W.nb and fortran_conventions.wl) and saves, per chapter, the finite remainders gqchri, aqchri,
   agchri, fourqwithqpp, fourqwithqp, the pole series and the real-virtual
   ingredients FVQEDgq, FVQEDgqz1dividedbyz, FLRVQCDγq, FLRVQCDγqz1, PqgRV, PqgNLO. The aux/ files replace the Dropbox paths of the notebook. *)
dir=Environment["NBDIR"]; out=Environment["NBOUT"]; 
Get[dir<>"/aux/replacements.m"]; Get[dir<>"/aux/soft_functions.m"]; Get[dir<>"/aux/splittings.m"]; Get[dir<>"/aux/simpl_functions.m"];
nb=Get[dir<>"/W_NNLO_polecheck_cleaned.nb"];
all=Cases[nb,Cell[c_,sty:("Input"|"Code"|"Title"|"Chapter"|"Section"|"Subsection"|"Subsubsection"|"Text"),___]:>{sty,c},Infinity];
skipQ[b_]:=!FreeQ[b,"Get"|"Export"|"Import"|"SystemOpen"|"HypExp`"|"User"|"Image"|"NumericArray"|"GraphicsBox"];
nlQ[x_]:=StringQ[x]&&(x==="\[IndentingNewLine]"||StringMatchQ[x,Whitespace|""]);
stmts[b_List]:=Flatten[stmts/@DeleteCases[b,_?nlQ],1];
stmts[RowBox[l_List]]/;MemberQ[l,"\[IndentingNewLine]"]:=Flatten[stmts/@(If[Length[#]==1,#[[1]],RowBox[#]]&/@DeleteCases[SplitBy[l,#==="\[IndentingNewLine]"&],{"\[IndentingNewLine]"..}]),1];
stmts[b_]:={b};
prevOut=Null; prevprevOut=Null;
eval[RowBox[{x_,";"}]]:=eval[x];
eval[b_]:=Module[{r}, r=ToExpression[b/.{"%%"->"prevprevOut","%"->"prevOut"},StandardForm]; prevprevOut=prevOut; prevOut=r; r];
saved=<||>; chap=0;
Do[c=all[[i]];
  If[c[[1]]==="Chapter", chap++; Print["== chapter ",chap]];
  If[MemberQ[{"Input","Code"},c[[1]]] && (i<12 || i>60) && !skipQ[c[[2]]],
   t=AbsoluteTiming[r=Last[eval/@stmts[c[[2]]/.BoxData[x_]:>x]]][[1]];
   If[t>5,Print["cell ",i,": ",Round[t]," s, mem ",Round[MaxMemoryUsed[]/10^6]," MB"]]];
  If[i==Length[all] || all[[i+1,1]]==="Chapter",
   Do[If[ToExpression[nm,InputForm,ValueQ], saved[{nm,chap}]=ToExpression[nm]],
     {nm,{"gqchri","aqchri","agchri","fourqwithqpp","fourqwithqp","total","total\[Gamma]qser","total\[Gamma]gser","total4qser"}}];
   (* real-virtual ingredients, for verify_W_polecheck_rv_fortran.wls *)
   Do[If[ToExpression[nm,InputForm,ValueQ], saved[{nm,chap}]=ToExpression[nm]],{nm,{"FVQEDgq","FVQEDgqz1dividedbyz","FLRVQCD\[Gamma]q","FLRVQCD\[Gamma]qz1"}}];
   Do[If[DownValues[Evaluate[ToExpression[nm]]]=!={}, saved[{nm,chap}]=DownValues[Evaluate[ToExpression[nm]]]],{nm,{"PqgRV","PqgNLO"}}];
   Print["end of chapter ",chap,": ",Keys[saved]]; Clear[gqchri,aqchri,agchri,fourqwithqpp,fourqwithqp,total,totalγqser,totalγgser,total4qser]],
 {i,Length[all]}];
Put[saved, out<>"/wpole_finite.m"];
Print["done"];
