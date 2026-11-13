(* usage: NBDIR=<math/subtraction-notebooks> NBOUT=<output dir> wolframscript -file run_W_nb.wls   -> NBOUT/wnb.mx, about 10 min *)
(* Re-evaluate the W.nb derivation (q qb' channel): poles and finite remainder.
   Statements are evaluated one by one and Out[] is filled, so that % works as in the front end. *)
dir=Environment["NBDIR"]; out=Environment["NBOUT"];
Get[dir<>"/aux/replacements.m"]; Get[dir<>"/aux/soft_functions.m"]; Get[dir<>"/aux/splittings.m"]; Get[dir<>"/aux/simpl_functions.m"];
nb=Get[dir<>"/W.nb"];
cells=Cases[nb,Cell[BoxData[b_],"Input",___]:>b,Infinity];
skipQ[b_]:=!FreeQ[b,"Get"|"Export"|"Import"|"mypath"|"HypExp`"|"User"];
last=First@FirstPosition[cells,b_/;!FreeQ[b,"finitecomplete"],{0},{1},Heads->False];
nlQ[x_]:=StringQ[x]&&(x==="\[IndentingNewLine]"||StringMatchQ[x,Whitespace|""]);
stmts[b_List]:=Flatten[stmts/@DeleteCases[b,_?nlQ],1];
stmts[RowBox[l_List]]/;MemberQ[l,"\[IndentingNewLine]"]:=Flatten[stmts/@(If[Length[#]==1,#[[1]],RowBox[#]]&/@DeleteCases[SplitBy[l,#==="\[IndentingNewLine]"&],{"\[IndentingNewLine]"..}]),1];
stmts[b_]:={b};
prevOut=Null; prevprevOut=Null;
eval[RowBox[{x_,";"}]]:=eval[x];
eval[b_]:=Module[{r}, r=Quiet[ToExpression[b/.{"%%"->"prevprevOut","%"->"prevOut"},StandardForm]]; prevprevOut=prevOut; prevOut=r; r];
Print["input cells: ",Length[cells],"  evaluating up to cell ",last];
res=<||>;
Do[
  If[skipQ[cells[[i]]],Null,
   t=AbsoluteTiming[r=Last[eval/@stmts[cells[[i]]]]][[1]];
   res[i]=r;
   If[t>5||i>=41,Print["cell ",i,": ",Round[t]," s, mem ",Round[MaxMemoryUsed[]/10^6]," MB, result ",Short[r,1]]]],
 {i,1,last}];
DumpSave[out<>"/wnb.mx",{TotQQSER,finitecomplete,finitesimpl,TotQQ,doublepole,singlepole,save,res}];
Print["finitesimpl: ",ByteCount[finitesimpl]," bytes"];
