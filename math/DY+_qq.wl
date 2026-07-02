(* Content-type: application/vnd.wolfram.mathematica *)

(*** Wolfram Notebook File ***)
(* http://www.wolfram.com/nb *)

(* CreatedBy='Wolfram 14.3' *)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[       154,          7]
NotebookDataLength[     66555,       1849]
NotebookOptionsPosition[     64476,       1817]
NotebookOutlinePosition[     64870,       1833]
CellTagsIndexPosition[     64827,       1830]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{
Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"CB", 
     RowBox[{"(", 
      RowBox[{"02.07", ".26"}], ")"}], "\[IndentingNewLine]", "Organisation", 
     " ", "of", " ", "my", " ", "finite", " ", "remainder", " ", "of", " ", 
     "the", " ", "subtraction", " ", "for", " ", "q", " ", "qp", " ", 
     "channel", " ", "in", " ", "Neutral", " ", "Drell"}], "-", 
    RowBox[{"Yan", " ", 
     RowBox[{"(", 
      RowBox[{"DY", "+"}], ")"}]}]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"q", " ", "qp"}], " ", "->", " ", 
    RowBox[{"l", " ", 
     RowBox[{"(", "neutrino", ")"}], " ", "lp", " ", 
     RowBox[{"(", 
      RowBox[{"positive", " ", "charge", " ", "electron"}], ")"}]}]}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{"(*", " ", 
   RowBox[{"Convention", ":", " ", 
    RowBox[{
     RowBox[{"func", "[", 
      RowBox[{
       RowBox[{"zi", " ", "q"}], ",", 
       RowBox[{"zj", " ", "qp"}], ",", "...", ",", "zi", ",", "zj"}], "]"}], " ",
      "means", " ", 
     RowBox[{"1", "/", "zi"}], "*", 
     RowBox[{"1", "/", "zj"}], "*", 
     RowBox[{"func", "[", 
      RowBox[{
       RowBox[{"zi", " ", "q"}], ",", 
       RowBox[{"zj", " ", "qp"}], ",", "..."}], "]"}]}]}], " ", 
   "*)"}]}]], "Input",
 CellChangeTimes->{{3.991962890183772*^9, 3.991962959541946*^9}, {
  3.991963010499596*^9, 3.991963013948811*^9}, {3.991963086190227*^9, 
  3.9919631678665743`*^9}},ExpressionUUID->"b30319e1-2140-427c-aa44-\
b44a4fc51180"],

Cell[BoxData[
 RowBox[{"(*", 
  RowBox[{"Useful", " ", "functions"}], "*)"}]], "Input",
 CellChangeTimes->{{3.991962403573283*^9, 
  3.991962407682333*^9}},ExpressionUUID->"f236d87a-10c6-4ab3-8c22-\
437a2885c0f1"],

Cell[BoxData[{
 RowBox[{
  RowBox[{
   RowBox[{"Pqqtilde", "[", 
    RowBox[{"z_", ",", "E1_"}], "]"}], ":=", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"4", " ", 
      RowBox[{"D1", "[", "z", "]"}]}], "-", 
     RowBox[{"2", 
      RowBox[{"(", 
       RowBox[{"1", "+", "z"}], ")"}], 
      RowBox[{"Log", "[", 
       RowBox[{"1", "-", "z"}], "]"}]}], "+", 
     RowBox[{"(", 
      RowBox[{"1", "-", "z"}], ")"}], "+", " ", 
     RowBox[{
      RowBox[{"Log", "[", 
       FractionBox[
        RowBox[{"4", " ", 
         SuperscriptBox["E1", "2"]}], 
        RowBox[{" ", 
         SuperscriptBox["mu", "2"]}]], "]"}], 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"2", " ", 
         RowBox[{"D0", "[", "z", "]"}]}], "-", 
        RowBox[{"(", 
         RowBox[{"1", "+", "z"}], ")"}]}], ")"}]}]}], ")"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"myfinalcalGQqQlp", "=", " ", 
   RowBox[{
    RowBox[{
     RowBox[{"Qq", "^", "2"}], 
     RowBox[{"(", 
      FractionBox[
       RowBox[{"2", " ", 
        SuperscriptBox["\[Pi]", "2"]}], "3"], ")"}]}], " ", "+", 
    RowBox[{
     RowBox[{"Qlp", "^", "2"}], " ", 
     RowBox[{"(", 
      RowBox[{
       FractionBox["13", "2"], "-", 
       SuperscriptBox["\[Pi]", "2"], "+", 
       SuperscriptBox[
        RowBox[{"Log", "[", "E4overEC", "]"}], "2"], "+", 
       FractionBox[
        RowBox[{"3", " ", 
         RowBox[{"Log", "[", "fourEC2overmu2", "]"}]}], "2"]}], ")"}]}], "+", 
    
    RowBox[{"Qlp", " ", "Qq", " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"-", 
        FractionBox[
         RowBox[{"2", " ", 
          SuperscriptBox["\[Pi]", "2"]}], "3"]}], "-", 
       RowBox[{"3", " ", 
        RowBox[{"Log", "[", "fourEC2overmu2", "]"}]}], "+", 
       RowBox[{"3", " ", 
        RowBox[{"Log", "[", 
         RowBox[{"eta", "[", 
          RowBox[{"1", ",", "4"}], "]"}], "]"}]}], "-", 
       RowBox[{"2", " ", 
        RowBox[{"Log", "[", "E4overEC", "]"}], " ", 
        RowBox[{"Log", "[", 
         RowBox[{"eta", "[", 
          RowBox[{"1", ",", "4"}], "]"}], "]"}]}], "+", 
       RowBox[{"2", " ", 
        RowBox[{"PolyLog", "[", 
         RowBox[{"2", ",", 
          RowBox[{"1", "-", 
           RowBox[{"eta", "[", 
            RowBox[{"1", ",", "4"}], "]"}]}]}], "]"}]}]}], ")"}]}], "+", " ", 
    
    RowBox[{"Qlp", " ", "Qqp", " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        RowBox[{"(", 
         RowBox[{"3", "-", 
          RowBox[{"2", " ", 
           RowBox[{"Log", "[", "E4overEC", "]"}]}]}], ")"}], " ", 
        RowBox[{"Log", "[", 
         RowBox[{"eta", "[", 
          RowBox[{"2", ",", "4"}], "]"}], "]"}]}], "+", 
       RowBox[{"2", " ", 
        RowBox[{"PolyLog", "[", 
         RowBox[{"2", ",", 
          RowBox[{"1", "-", 
           RowBox[{"eta", "[", 
            RowBox[{"2", ",", "4"}], "]"}]}]}], "]"}]}]}], ")"}]}]}]}], 
  ";"}]}], "Input",
 CellChangeTimes->{{3.9919624301534967`*^9, 3.991962430277067*^9}, {
   3.99196330975391*^9, 3.991963311119706*^9}, 3.991963342320128*^9, {
   3.991963865762113*^9, 3.991963865895911*^9}, 3.9919648474616756`*^9, {
   3.9919654631407843`*^9, 3.9919654774445953`*^9}, {3.991966043738991*^9, 
   3.9919660439917097`*^9}},
 CellLabel->
  "In[518]:=",ExpressionUUID->"d60e41d1-3362-4c88-bfd0-0fb736a2ac8d"],

Cell[BoxData[
 RowBox[{
  RowBox[{"genGnloQCD", "=", 
   RowBox[{
    RowBox[{
     FractionBox["1", "3"], " ", 
     SuperscriptBox["\[Pi]", "2"], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        RowBox[{"-", "2"}], " ", 
        SuperscriptBox["Qlp", "2"]}], "-", 
       RowBox[{"Qlp", " ", "Qq"}], "-", 
       RowBox[{"Qlp", " ", "Qqp"}], "-", 
       RowBox[{"2", " ", "Qq", " ", "Qqp"}]}], ")"}]}], "+", 
    RowBox[{"Qlp", " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"2", " ", "Qlp"}], "-", "Qq", "-", "Qqp"}], ")"}], " ", 
     SuperscriptBox[
      RowBox[{"Log", "[", "twoE4onMu", "]"}], "2"]}], "-", 
    RowBox[{"Qlp", " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"2", " ", "Qlp"}], "-", 
       RowBox[{"3", " ", "Qq"}], "-", 
       RowBox[{"3", " ", "Qqp"}]}], ")"}], " ", 
     SuperscriptBox[
      RowBox[{"Log", "[", "twoEConMu", "]"}], "2"]}], "+", 
    RowBox[{
     RowBox[{"Log", "[", "twoEConMu", "]"}], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"3", " ", "Qlp", " ", "Qq"}], "+", 
       RowBox[{"3", " ", "Qlp", " ", "Qqp"}], "-", 
       RowBox[{"6", " ", "Qq", " ", "Qqp"}], "+", 
       RowBox[{"2", " ", "Qlp", " ", "Qq", " ", 
        RowBox[{"Log", "[", 
         RowBox[{"eta", "[", 
          RowBox[{"1", ",", "4"}], "]"}], "]"}]}], "+", 
       RowBox[{"2", " ", "Qlp", " ", "Qqp", " ", 
        RowBox[{"Log", "[", 
         RowBox[{"eta", "[", 
          RowBox[{"2", ",", "4"}], "]"}], "]"}]}]}], ")"}]}], "+", 
    RowBox[{
     RowBox[{"Log", "[", "twoE4onMu", "]"}], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"2", " ", "Qlp", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "Qq"}], "-", "Qqp"}], ")"}], " ", 
        RowBox[{"Log", "[", "twoEConMu", "]"}]}], "-", 
       RowBox[{"Qlp", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"3", " ", "Qlp"}], "-", 
          RowBox[{"3", " ", "Qq"}], "-", 
          RowBox[{"3", " ", "Qqp"}], "+", 
          RowBox[{"2", " ", "Qq", " ", 
           RowBox[{"Log", "[", 
            RowBox[{"eta", "[", 
             RowBox[{"1", ",", "4"}], "]"}], "]"}]}], "+", 
          RowBox[{"2", " ", "Qqp", " ", 
           RowBox[{"Log", "[", 
            RowBox[{"eta", "[", 
             RowBox[{"2", ",", "4"}], "]"}], "]"}]}]}], ")"}]}]}], ")"}]}], "+", 
    RowBox[{
     FractionBox["1", "2"], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"13", " ", 
        SuperscriptBox["Qlp", "2"]}], "-", 
       RowBox[{"6", " ", "Qq", " ", "Qqp", " ", 
        RowBox[{"Log", "[", 
         RowBox[{"eta", "[", 
          RowBox[{"1", ",", "2"}], "]"}], "]"}]}], "+", 
       RowBox[{"6", " ", "Qlp", " ", "Qq", " ", 
        RowBox[{"Log", "[", 
         RowBox[{"eta", "[", 
          RowBox[{"1", ",", "4"}], "]"}], "]"}]}], "+", 
       RowBox[{"6", " ", "Qlp", " ", "Qqp", " ", 
        RowBox[{"Log", "[", 
         RowBox[{"eta", "[", 
          RowBox[{"2", ",", "4"}], "]"}], "]"}]}], "-", 
       RowBox[{"4", " ", "Qq", " ", "Qqp", " ", 
        RowBox[{"PolyLog", "[", 
         RowBox[{"2", ",", 
          RowBox[{"1", "-", 
           RowBox[{"eta", "[", 
            RowBox[{"1", ",", "2"}], "]"}]}]}], "]"}]}], "+", 
       RowBox[{"4", " ", "Qlp", " ", "Qq", " ", 
        RowBox[{"PolyLog", "[", 
         RowBox[{"2", ",", 
          RowBox[{"1", "-", 
           RowBox[{"eta", "[", 
            RowBox[{"1", ",", "4"}], "]"}]}]}], "]"}]}], "+", 
       RowBox[{"4", " ", "Qlp", " ", "Qqp", " ", 
        RowBox[{"PolyLog", "[", 
         RowBox[{"2", ",", 
          RowBox[{"1", "-", 
           RowBox[{"eta", "[", 
            RowBox[{"2", ",", "4"}], "]"}]}]}], "]"}]}]}], ")"}]}]}]}], 
  ";"}]], "Input",
 CellChangeTimes->{{3.991965089864251*^9, 3.991965091840464*^9}, {
  3.9919653115072393`*^9, 3.991965317004118*^9}},
 CellLabel->
  "In[510]:=",ExpressionUUID->"712dedb1-0bdc-42cb-bb10-9d500316a9b5"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", 
   RowBox[{
   "equivalent", " ", "to", " ", "myfinalcalGQqQlp", " ", "but", " ", "with", 
    " ", "general", " ", 
    RowBox[{"charges", ":", " ", 
     RowBox[{
     "we", " ", "have", " ", "to", " ", "be", " ", "careful", " ", "with", " ",
       "the", " ", 
      RowBox[{"Qq", "^", "2"}], " ", 
      RowBox[{"shift", ".", " ", "This"}], " ", "allowed", " ", "me", " ", 
      "to", " ", "extend", " ", "also", " ", "to", " ", "leg2", " ", 
      "easily"}]}]}], "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"leg1finalcalG", "=", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{"genGnloQCD", "-", 
        RowBox[{"3", 
         RowBox[{"Qq", "^", "2"}], 
         RowBox[{"Log", "[", "fourEC2overmu2", "]"}]}]}], ")"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"Ql", ":>", "0"}], ",", 
        RowBox[{
         RowBox[{"Log", "[", "twoEConMu", "]"}], ":>", 
         RowBox[{
          RowBox[{"1", "/", "2"}], 
          RowBox[{"Log", "[", "fourEC2overmu2", "]"}]}]}], ",", 
        RowBox[{
         RowBox[{"Log", "[", "twoE4onMu", "]"}], ":>", 
         RowBox[{
          RowBox[{"Log", "[", "E4overEC", "]"}], "+", 
          RowBox[{"Log", "[", "twoEConMu", "]"}]}]}], ",", 
        RowBox[{
         RowBox[{"eta", "[", 
          RowBox[{"1", ",", "2"}], "]"}], ":>", "1"}]}], "}"}]}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"leg2finalcalG", "=", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{"genGnloQCD", "-", 
        RowBox[{"3", 
         RowBox[{"Qqp", "^", "2"}], 
         RowBox[{"Log", "[", "fourEC2overmu2", "]"}]}]}], ")"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"Ql", ":>", "0"}], ",", 
        RowBox[{
         RowBox[{"Log", "[", "twoEConMu", "]"}], ":>", 
         RowBox[{
          RowBox[{"1", "/", "2"}], 
          RowBox[{"Log", "[", "fourEC2overmu2", "]"}]}]}], ",", 
        RowBox[{
         RowBox[{"Log", "[", "twoE4onMu", "]"}], ":>", 
         RowBox[{
          RowBox[{"Log", "[", "E4overEC", "]"}], "+", 
          RowBox[{"Log", "[", "twoEConMu", "]"}]}]}], ",", 
        RowBox[{
         RowBox[{"eta", "[", 
          RowBox[{"1", ",", "2"}], "]"}], ":>", "1"}]}], "}"}]}]}], 
    ";"}]}]}]], "Input",
 CellChangeTimes->{{3.99196607289642*^9, 3.991966188924285*^9}, {
  3.991966264163134*^9, 3.991966277243388*^9}, {3.991966885776874*^9, 
  3.991966917259005*^9}},
 CellLabel->
  "In[603]:=",ExpressionUUID->"e93843a7-a8e6-4d19-b073-2630ae096939"],

Cell[BoxData[
 RowBox[{"(*", "FLVQCDfin", "*)"}]], "Input",
 CellChangeTimes->{{3.9919624470469646`*^9, 3.991962449922947*^9}, {
  3.991962562051062*^9, 
  3.991962562289609*^9}},ExpressionUUID->"275029be-23b0-4453-a5ca-\
ef94fbc4167c"],

Cell[BoxData[
 RowBox[{
  RowBox[{"FINITEDYp", "[", "1", "]"}], "=", 
  RowBox[{
   RowBox[{"(", 
    RowBox[{
     RowBox[{
      SuperscriptBox["Qqp", "2"], " ", 
      RowBox[{"FLVqcdfin", "[", 
       RowBox[{"q", ",", 
        RowBox[{"qp", " ", "z"}], ",", "l", ",", "lp", ",", "z"}], "]"}]}], "+", 
     RowBox[{
      SuperscriptBox["Qq", "2"], " ", 
      RowBox[{"FLVqcdfin", "[", 
       RowBox[{
        RowBox[{"q", " ", "z"}], ",", "qp", ",", "l", ",", "lp", ",", "z"}], 
       "]"}]}]}], ")"}], " ", 
   RowBox[{"(", 
    RowBox[{"1", "-", "z", "+", 
     RowBox[{"4", " ", 
      RowBox[{"D1", "[", "z", "]"}]}], "+", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{"1", "+", "z", "-", 
        RowBox[{"2", " ", 
         RowBox[{"D0", "[", "z", "]"}]}]}], ")"}], " ", 
      RowBox[{"Log", "[", 
       FractionBox["mu2", "s"], "]"}]}], "-", 
     FractionBox[
      RowBox[{"2", " ", 
       RowBox[{"Log", "[", "z", "]"}]}], 
      RowBox[{"1", "-", "z"}]], "+", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{"1", "+", "z"}], ")"}], " ", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{
         RowBox[{"-", "2"}], " ", 
         RowBox[{"Log", "[", 
          RowBox[{"1", "-", "z"}], "]"}]}], "+", 
        RowBox[{"Log", "[", "z", "]"}]}], ")"}]}]}], ")"}]}]}]], "Input",
 CellChangeTimes->{{3.991962470281146*^9, 3.99196254999155*^9}, {
   3.991962591636293*^9, 3.99196259215806*^9}, 3.991963423951234*^9, 
   3.991963516633519*^9, {3.991963727519718*^9, 
   3.991963728100671*^9}},ExpressionUUID->"20fd7f69-dcb4-4ffa-9864-\
3140d26158cd"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"FINITEDYp", "[", "2", "]"}], "=", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{
      FractionBox["1", "6"], " ", 
      RowBox[{"(", 
       RowBox[{"39", "-", 
        RowBox[{"4", " ", 
         SuperscriptBox["\[Pi]", "2"]}]}], ")"}], " ", 
      SuperscriptBox["Qlp", "2"], " ", 
      RowBox[{"FLVqcdfin", "[", 
       RowBox[{"q", ",", "qp", ",", "l", ",", "lp"}], "]"}]}], "+", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{
        RowBox[{
         FractionBox["3", "2"], " ", 
         SuperscriptBox["Qlp", "2"], " ", 
         RowBox[{"FLVqcdfin", "[", 
          RowBox[{"q", ",", "qp", ",", "l", ",", "lp"}], "]"}]}], "-", 
        RowBox[{"3", " ", "QlpQq", " ", 
         RowBox[{"FLVqcdfin", "[", 
          RowBox[{"q", ",", "qp", ",", "l", ",", "lp"}], "]"}]}], "-", 
        RowBox[{"3", " ", "QlpQqp", " ", 
         RowBox[{"FLVqcdfin", "[", 
          RowBox[{"q", ",", "qp", ",", "l", ",", "lp"}], "]"}]}], "+", 
        RowBox[{"3", " ", "QqQqp", " ", 
         RowBox[{"FLVqcdfin", "[", 
          RowBox[{"q", ",", "qp", ",", "l", ",", "lp"}], "]"}]}]}], ")"}], " ", 
      RowBox[{"Log", "[", 
       FractionBox["mu2", "s"], "]"}]}], "+", 
     RowBox[{
      FractionBox["1", "6"], " ", "QqQqp", " ", 
      RowBox[{"FLVqcdfin", "[", 
       RowBox[{"q", ",", "qp", ",", "l", ",", "lp"}], "]"}], " ", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{
         RowBox[{"-", "4"}], " ", 
         SuperscriptBox["\[Pi]", "2"]}], "-", 
        RowBox[{"18", " ", 
         RowBox[{"Log", "[", "eta12", "]"}]}], "-", 
        RowBox[{"12", " ", 
         RowBox[{"PolyLog", "[", 
          RowBox[{"2", ",", 
           RowBox[{"1", "-", "eta12"}]}], "]"}]}]}], ")"}]}], "+", 
     RowBox[{
      FractionBox["1", "6"], " ", "QlpQq", " ", 
      RowBox[{"FLVqcdfin", "[", 
       RowBox[{"q", ",", "qp", ",", "l", ",", "lp"}], "]"}], " ", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{
         RowBox[{"-", "2"}], " ", 
         SuperscriptBox["\[Pi]", "2"]}], "+", 
        RowBox[{"18", " ", 
         RowBox[{"Log", "[", "eta14", "]"}]}], "+", 
        RowBox[{"12", " ", 
         RowBox[{"PolyLog", "[", 
          RowBox[{"2", ",", 
           RowBox[{"1", "-", "eta14"}]}], "]"}]}]}], ")"}]}], "+", 
     RowBox[{
      FractionBox["1", "6"], " ", "QlpQqp", " ", 
      RowBox[{"FLVqcdfin", "[", 
       RowBox[{"q", ",", "qp", ",", "l", ",", "lp"}], "]"}], " ", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{
         RowBox[{"-", "2"}], " ", 
         SuperscriptBox["\[Pi]", "2"]}], "+", 
        RowBox[{"18", " ", 
         RowBox[{"Log", "[", "eta24", "]"}]}], "+", 
        RowBox[{"12", " ", 
         RowBox[{"PolyLog", "[", 
          RowBox[{"2", ",", 
           RowBox[{"1", "-", "eta24"}]}], "]"}]}]}], ")"}]}]}], ")"}]}], 
  ";"}]], "Input",
 CellChangeTimes->{{3.991962514245973*^9, 3.991962548417899*^9}, {
   3.991962589577559*^9, 3.9919625899964848`*^9}, 3.991963541221154*^9, {
   3.9919637255256977`*^9, 
   3.99196372579257*^9}},ExpressionUUID->"c7cac081-5f2c-426e-8998-\
5bddfd52e026"],

Cell[BoxData[
 RowBox[{"(*", "FLVEWfin", "*)"}]], "Input",
 CellChangeTimes->{{3.991962556051277*^9, 
  3.991962567675682*^9}},ExpressionUUID->"3059f3cd-9fd6-4ecd-8b35-\
6d8bbb42698f"],

Cell[BoxData[
 RowBox[{
  RowBox[{"FINITEDY0p", "[", "3", "]"}], "=", 
  RowBox[{"CF", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"FLVewfin", "[", 
      RowBox[{"q", ",", 
       RowBox[{"qp", " ", "z"}], ",", "l", ",", "lp", ",", "z"}], "]"}], "+", 
     
     RowBox[{"FLVewfin", "[", 
      RowBox[{
       RowBox[{"q", " ", "z"}], ",", "qp", ",", "l", ",", "lp", ",", "z"}], 
      "]"}]}], ")"}], " ", 
   RowBox[{"(", 
    RowBox[{"1", "-", "z", "+", 
     RowBox[{"4", " ", 
      RowBox[{"D1", "[", "z", "]"}]}], "+", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{"1", "+", "z", "-", 
        RowBox[{"2", " ", 
         RowBox[{"D0", "[", "z", "]"}]}]}], ")"}], " ", 
      RowBox[{"Log", "[", 
       FractionBox["mu2", "s"], "]"}]}], "-", 
     RowBox[{"2", " ", 
      RowBox[{"(", 
       RowBox[{"1", "+", "z"}], ")"}], " ", 
      RowBox[{"Log", "[", 
       RowBox[{"1", "-", "z"}], "]"}]}], "-", 
     FractionBox[
      RowBox[{"2", " ", 
       RowBox[{"Log", "[", "z", "]"}]}], 
      RowBox[{"1", "-", "z"}]], "+", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{"1", "+", "z"}], ")"}], " ", 
      RowBox[{"Log", "[", "z", "]"}]}]}], ")"}]}]}]], "Input",
 CellChangeTimes->{{3.991962583429284*^9, 3.991962616875204*^9}, {
  3.991962708156433*^9, 3.991962712548923*^9}, {3.991963656946707*^9, 
  3.991963723470315*^9}},ExpressionUUID->"85b59acb-5a88-43b6-b790-\
3d4ef61b48b8"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"FINITEDYp", "[", "4", "]"}], "=", 
   RowBox[{"CF", 
    RowBox[{"(", 
     RowBox[{
      FractionBox[
       RowBox[{"2", " ", 
        SuperscriptBox["\[Pi]", "2"]}], "3"], "-", 
      RowBox[{"3", " ", 
       RowBox[{"Log", "[", 
        FractionBox["mu2", "s"], "]"}]}]}], ")"}], 
    RowBox[{"FLVewfin", "[", 
     RowBox[{"q", ",", "qp", ",", "l", ",", "lp"}], "]"}]}]}], " ", 
  ";"}]], "Input",
 CellChangeTimes->{{3.991962617932856*^9, 3.991962664248831*^9}, {
  3.9919637315418377`*^9, 
  3.991963753061118*^9}},ExpressionUUID->"a30eb47c-ddbe-487e-b668-\
222da2ec68b2"],

Cell[BoxData[
 RowBox[{"(*", 
  RowBox[{"Double", " ", "Boosted"}], "*)"}]], "Input",
 CellChangeTimes->{{3.991962671732389*^9, 
  3.991962680416296*^9}},ExpressionUUID->"894ed87b-f744-4e51-a547-\
41a2bae4dcac"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"FINITEDY0", "[", "5", "]"}], "=", 
   RowBox[{"CF", 
    RowBox[{"(", " ", 
     RowBox[{
      SuperscriptBox["Qq", "2"], " ", "+", 
      SuperscriptBox["Qqp", "2"]}], ")"}], " ", 
    RowBox[{"Identity", "[", 
     RowBox[{"FLM", "[", 
      RowBox[{
       RowBox[{"z1", ".", "p1"}], ",", 
       RowBox[{"z2", ".", "p2"}], ",", "p3", ",", "p4", ",", "z1", ",", 
       "z2"}], "]"}], "]"}], " ", 
    RowBox[{"(", "\[IndentingNewLine]", 
     RowBox[{"(*", "AA", "*)"}], "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{
       RowBox[{"+", 
        RowBox[{"(", 
         RowBox[{"1", "-", "z1", "-", 
          RowBox[{"2", " ", 
           RowBox[{"(", 
            RowBox[{"1", "+", "z1"}], ")"}], " ", 
           RowBox[{"Log", "[", 
            RowBox[{"1", "-", "z1"}], "]"}]}], "+", 
          RowBox[{
           RowBox[{"(", 
            RowBox[{"1", "+", "z1", "-", 
             RowBox[{"2", 
              RowBox[{"D0", "[", "z1", "]"}]}]}], ")"}], " ", 
           RowBox[{"Log", "[", "z1", "]"}]}], "+", 
          RowBox[{"4", " ", 
           RowBox[{"D1", "[", "z1", "]"}]}], "+", 
          RowBox[{
           RowBox[{"(", 
            RowBox[{"1", "+", "z1", "-", 
             RowBox[{"2", 
              RowBox[{"D0", "[", "z1", "]"}]}]}], ")"}], 
           RowBox[{"Log", "[", 
            FractionBox[
             SuperscriptBox["\[Mu]", "2"], "s"], "]"}]}]}], ")"}]}], " ", 
       RowBox[{"(", 
        RowBox[{"1", "-", "z2", "-", 
         RowBox[{"2", 
          RowBox[{"(", 
           RowBox[{"1", "+", "z2"}], ")"}], " ", 
          RowBox[{"Log", "[", 
           RowBox[{"1", "-", "z2"}], "]"}]}], "+", 
         RowBox[{
          RowBox[{"(", 
           RowBox[{"1", "+", "z2", "-", 
            RowBox[{"2", 
             RowBox[{"D0", "[", "z2", "]"}]}]}], ")"}], " ", 
          RowBox[{"Log", "[", "z2", "]"}]}], "+", 
         RowBox[{"4", 
          RowBox[{"D1", "[", "z2", "]"}]}], "+", 
         RowBox[{
          RowBox[{"(", 
           RowBox[{"1", "+", "z2", "-", 
            RowBox[{"2", 
             RowBox[{"D0", "[", "z2", "]"}]}]}], ")"}], 
          RowBox[{"Log", "[", 
           FractionBox[
            SuperscriptBox["\[Mu]", "2"], "s"], "]"}]}]}], ")"}]}], "+", 
      "\[IndentingNewLine]", 
      RowBox[{"(*", "BB", "*)"}], "\[IndentingNewLine]", 
      RowBox[{
       RowBox[{"(", 
        RowBox[{
         RowBox[{
          RowBox[{"(", 
           RowBox[{"1", "-", "z1", "-", " ", 
            RowBox[{"2", 
             RowBox[{"(", 
              RowBox[{"1", "+", "z1"}], ")"}], " ", 
             RowBox[{"Log", "[", 
              RowBox[{"1", "-", "z1"}], "]"}]}], "+", 
            RowBox[{
             RowBox[{"(", 
              RowBox[{"1", "+", "z1", "-", 
               RowBox[{"2", 
                RowBox[{"D0", "[", "z1", "]"}]}]}], ")"}], " ", 
             RowBox[{"Log", "[", "z1", "]"}]}], "+", " ", 
            RowBox[{"4", " ", 
             RowBox[{"D1", "[", "z1", "]"}]}]}], ")"}], 
          RowBox[{"Log", "[", "z1", "]"}]}], "+", 
         RowBox[{
          RowBox[{"(", 
           RowBox[{"1", "+", "z1", "-", 
            RowBox[{"2", 
             RowBox[{"D0", "[", "z1", "]"}]}]}], ")"}], 
          RowBox[{"Log", "[", "z1", "]"}], 
          RowBox[{"Log", "[", 
           FractionBox[
            SuperscriptBox["\[Mu]", "2"], "s"], "]"}]}]}], ")"}], 
       RowBox[{"(", 
        RowBox[{"1", "+", "z2", "-", 
         RowBox[{"2", 
          RowBox[{"D0", "[", "z2", "]"}]}]}], ")"}]}], "+", 
      "\[IndentingNewLine]", 
      RowBox[{"(*", "CC", "*)"}], "\[IndentingNewLine]", 
      RowBox[{
       RowBox[{"(", 
        RowBox[{
         RowBox[{
          RowBox[{"(", 
           RowBox[{"1", "-", "z2", "-", 
            RowBox[{"2", 
             RowBox[{"(", 
              RowBox[{"1", "+", "z2"}], ")"}], " ", 
             RowBox[{"Log", "[", 
              RowBox[{"1", "-", "z2"}], "]"}]}], "+", 
            RowBox[{
             RowBox[{"(", 
              RowBox[{"1", "+", "z2", "-", 
               RowBox[{"2", 
                RowBox[{"D0", "[", "z2", "]"}]}]}], ")"}], " ", 
             RowBox[{"Log", "[", "z2", "]"}]}], "+", " ", 
            RowBox[{"4", " ", 
             RowBox[{"D1", "[", "z2", "]"}]}]}], ")"}], 
          RowBox[{"Log", "[", "z2", "]"}]}], "+", 
         RowBox[{
          RowBox[{"(", 
           RowBox[{"1", "+", "z2", "-", 
            RowBox[{"2", 
             RowBox[{"D0", "[", "z2", "]"}]}]}], ")"}], 
          RowBox[{"Log", "[", "z2", "]"}], 
          RowBox[{"Log", "[", 
           FractionBox[
            SuperscriptBox["\[Mu]", "2"], "s"], "]"}]}]}], ")"}], 
       RowBox[{"(", 
        RowBox[{"1", "+", "z1", "-", 
         RowBox[{"2", 
          RowBox[{"D0", "[", "z1", "]"}]}]}], ")"}]}], " ", "+", 
      "\[IndentingNewLine]", 
      RowBox[{"(*", "DD", "*)"}], "\[IndentingNewLine]", 
      RowBox[{
       RowBox[{"(", 
        RowBox[{"1", "+", "z1", "-", 
         RowBox[{"2", " ", 
          RowBox[{"D0", "[", "z1", "]"}]}]}], ")"}], 
       RowBox[{"Log", "[", "z1", "]"}], " ", " ", 
       RowBox[{"(", 
        RowBox[{"1", "+", "z2", "-", 
         RowBox[{"2", " ", 
          RowBox[{"D0", "[", "z2", "]"}]}]}], ")"}], 
       RowBox[{"Log", "[", "z2", "]"}]}]}], ")"}]}]}], ";"}]], "Input",
 CellChangeTimes->{{3.9919626814994383`*^9, 3.991962705485297*^9}, {
  3.991963789699566*^9, 
  3.991963807604961*^9}},ExpressionUUID->"836b927f-6ea5-4518-859d-\
47c1b6089cad"],

Cell[BoxData[
 RowBox[{"(*", 
  RowBox[{"Single", " ", 
   RowBox[{"boosted", ":", " ", 
    RowBox[{
    "derived", " ", "for", " ", "leg", " ", "1", " ", "but", " ", "similar", " ",
      "for", " ", "leg", " ", "2", " ", "with", " ", "careful", " ", 
     "replacement", " ", "of", " ", "index", " ", "numbers", " ", "and", " ", 
     "charges"}]}]}], "*)"}]], "Input",
 CellChangeTimes->{{3.99196272486541*^9, 3.991962727143035*^9}, {
  3.99196276696769*^9, 3.9919627896524*^9}, {3.991963892413553*^9, 
  3.991963897033119*^9}},ExpressionUUID->"b88ff76d-0dfe-450b-89e9-\
9e20dd8cc655"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"FINITEDYp", "[", "6", "]"}], "=", " ", 
   RowBox[{
    RowBox[{"CF", " ", 
     RowBox[{"FLM", "[", 
      RowBox[{
       RowBox[{"q", " ", "z"}], ",", "qp", ",", "l", ",", "lp", ",", "z"}], 
      "]"}], 
     RowBox[{"(", 
      RowBox[{"1", "-", "z", "-", 
       RowBox[{"2", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "z"}], ")"}], " ", 
        RowBox[{"Log", "[", 
         RowBox[{"1", "-", "z"}], "]"}]}], "+", 
       RowBox[{"4", " ", 
        RowBox[{"D1", "[", "z", "]"}]}], "+", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{"1", "+", "z", "-", 
          RowBox[{"2", " ", 
           RowBox[{"D0", "[", "z", "]"}]}]}], ")"}], " ", 
        RowBox[{"Log", "[", 
         FractionBox[
          SuperscriptBox["\[Mu]", "2"], 
          RowBox[{"4", 
           SuperscriptBox["EC", "2"]}]], "]"}]}]}], " ", ")"}], 
     "leg1finalcalG"}], "+", "\[IndentingNewLine]", 
    RowBox[{"CF", "  ", 
     RowBox[{"FLM", "[", 
      RowBox[{
       RowBox[{"q", " ", "z"}], ",", "qp", ",", "l", ",", "lp", ",", "z"}], 
      "]"}], 
     RowBox[{"(", 
      RowBox[{"1", "-", "z", "-", 
       RowBox[{"2", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "z"}], ")"}], " ", 
        RowBox[{"Log", "[", 
         RowBox[{"1", "-", "z"}], "]"}]}], "+", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{"1", "+", "z", "-", 
          RowBox[{"2", " ", 
           RowBox[{"D0", "[", "z", "]"}]}]}], ")"}], " ", 
        RowBox[{"Log", "[", "z", "]"}]}], "+", 
       RowBox[{"4", " ", 
        RowBox[{"D1", "[", "z", "]"}]}], "+", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{"1", "+", "z", "-", 
          RowBox[{"2", " ", 
           RowBox[{"D0", "[", "z", "]"}]}]}], ")"}], " ", 
        RowBox[{"Log", "[", 
         FractionBox[
          SuperscriptBox["\[Mu]", "2"], "s"], "]"}]}]}], ")"}], 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"-", "2"}], "Qlp", " ", "Qq", " ", 
       RowBox[{"Log", "[", "z", "]"}], 
       RowBox[{"Log", "[", 
        RowBox[{"s14", "/", "s12"}], "]"}]}], ")"}]}], "+", 
    "\[IndentingNewLine]", 
    RowBox[{"CF", "   ", 
     RowBox[{"FLM", "[", 
      RowBox[{
       RowBox[{"q", " ", "z"}], ",", "qp", ",", "l", ",", "lp", ",", "z"}], 
      "]"}], 
     RowBox[{"Qq", "^", "2"}], 
     RowBox[{"(", " ", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{
         FractionBox["1", 
          RowBox[{"1", "-", "z"}]], 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"8", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"4", " ", 
               RowBox[{"Zeta", "[", "2", "]"}]}], "+", 
              RowBox[{"PolyLog", "[", 
               RowBox[{"2", ",", "z"}], "]"}], "+", 
              RowBox[{
               FractionBox["3", "4"], " ", 
               SuperscriptBox[
                RowBox[{"Log", "[", 
                 RowBox[{"1", "-", "z"}], "]"}], "2"]}]}], ")"}], 
            RowBox[{"Log", "[", "z", "]"}]}], "+", 
           RowBox[{"10", " ", 
            SuperscriptBox[
             RowBox[{"Log", "[", "z", "]"}], "3"]}], "+", 
           RowBox[{"12", 
            RowBox[{"Log", "[", 
             RowBox[{"1", "-", "z"}], "]"}], " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"Zeta", "[", "2", "]"}], "-", " ", 
              RowBox[{"PolyLog", "[", 
               RowBox[{"2", ",", "z"}], "]"}], "+", 
              RowBox[{"2", 
               RowBox[{"Log", "[", "z", "]"}]}]}], ")"}]}], "-", 
           RowBox[{"4", " ", 
            RowBox[{"PolyLog", "[", 
             RowBox[{"3", ",", 
              RowBox[{"1", "-", "z"}]}], "]"}]}], "-", 
           RowBox[{"32", " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"PolyLog", "[", 
               RowBox[{"3", ",", "z"}], "]"}], "-", 
              RowBox[{"Zeta", "[", "3", "]"}]}], ")"}]}]}], ")"}]}], "+", 
        "\[IndentingNewLine]", 
        RowBox[{"4", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{"-", "3"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}], "+", 
           RowBox[{"4", " ", 
            RowBox[{"Log", "[", 
             RowBox[{"1", "-", "z"}], "]"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}], "+", 
           RowBox[{"2", " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"PolyLog", "[", 
               RowBox[{"2", ",", "z"}], "]"}], "-", 
              RowBox[{"Zeta", "[", "2", "]"}]}], ")"}]}]}], ")"}], 
         FractionBox[
          RowBox[{"Log", "[", "z", "]"}], 
          RowBox[{"1", "-", "z"}]]}], "-", 
        RowBox[{"4", 
         FractionBox[
          SuperscriptBox[
           RowBox[{"Log", "[", "z", "]"}], "3"], 
          RowBox[{"1", "-", "z"}]]}], "+", 
        RowBox[{"12", 
         RowBox[{"D0", "[", "z", "]"}], " ", 
         RowBox[{"(", 
          RowBox[{"1", "-", 
           RowBox[{
            FractionBox["4", "9"], 
            RowBox[{"Log", "[", "z", "]"}]}]}], ")"}], 
         SuperscriptBox[
          RowBox[{"Log", "[", "z", "]"}], "2"]}], "-", 
        RowBox[{"16", " ", 
         RowBox[{"D2", "[", "z", "]"}], " ", 
         RowBox[{"Log", "[", "z", "]"}]}], "+", "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{"-", 
             FractionBox["2", "3"]}], " ", 
            SuperscriptBox["\[Pi]", "2"], " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}]}], "+", 
           RowBox[{"2", " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"-", "6"}], "+", 
              RowBox[{"5", " ", "z"}]}], ")"}]}], "+", 
           RowBox[{"12", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            SuperscriptBox[
             RowBox[{"Log", "[", 
              RowBox[{"1", "-", "z"}], "]"}], "2"]}], "+", 
           RowBox[{"4", " ", 
            RowBox[{"(", 
             RowBox[{"2", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}], "-", 
           RowBox[{
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            SuperscriptBox[
             RowBox[{"Log", "[", "z", "]"}], "2"]}], "+", 
           RowBox[{
            RowBox[{"Log", "[", 
             RowBox[{"1", "-", "z"}], "]"}], " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"8", " ", 
               RowBox[{"(", 
                RowBox[{"2", "+", "z"}], ")"}]}], "-", 
              RowBox[{"8", " ", 
               RowBox[{"(", 
                RowBox[{"1", "+", "z"}], ")"}], " ", 
               RowBox[{"Log", "[", "z", "]"}]}]}], ")"}]}], "+", 
           RowBox[{
            FractionBox["4", "3"], " ", 
            SuperscriptBox["\[Pi]", "2"], " ", 
            RowBox[{"D0", "[", "z", "]"}]}], "-", 
           RowBox[{"24", 
            RowBox[{"D1", "[", "z", "]"}]}], "-", 
           RowBox[{"24", 
            RowBox[{"D2", "[", "z", "]"}]}]}], ")"}], 
         RowBox[{"Log", "[", "z", "]"}]}], "+", 
        RowBox[{
         SuperscriptBox[
          RowBox[{"Log", "[", "z", "]"}], "2"], 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{"-", "4"}], " ", 
            RowBox[{"(", 
             RowBox[{"2", "+", "z"}], ")"}]}], "-", 
           RowBox[{"4", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", 
             RowBox[{"1", "-", "z"}], "]"}]}], "+", 
           RowBox[{"3", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}]}], ")"}]}], "+", 
        "\[IndentingNewLine]", 
        RowBox[{"2", " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "9"}], "+", 
           RowBox[{"8", " ", "z"}]}], ")"}]}], "-", 
        RowBox[{"8", " ", 
         RowBox[{"(", 
          RowBox[{"1", "+", "z"}], ")"}], " ", 
         SuperscriptBox[
          RowBox[{"Log", "[", 
           RowBox[{"1", "-", "z"}], "]"}], "3"]}], "+", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{"5", "-", 
           RowBox[{"9", " ", "z"}]}], ")"}], " ", 
         RowBox[{"Log", "[", "z", "]"}]}], "+", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{"2", "-", "z"}], ")"}], " ", 
         SuperscriptBox[
          RowBox[{"Log", "[", "z", "]"}], "2"]}], "-", 
        RowBox[{
         FractionBox["13", "6"], " ", 
         RowBox[{"(", 
          RowBox[{"1", "+", "z"}], ")"}], " ", 
         SuperscriptBox[
          RowBox[{"Log", "[", "z", "]"}], "3"]}], "+", 
        RowBox[{
         SuperscriptBox["\[Pi]", "2"], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"3", " ", 
            RowBox[{"(", 
             RowBox[{"1", "-", "z"}], ")"}]}], "+", 
           RowBox[{
            FractionBox["5", "3"], " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", 
             RowBox[{"1", "-", "z"}], "]"}]}], "-", 
           RowBox[{
            FractionBox["8", "3"], " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}]}], ")"}]}], "+", 
        RowBox[{
         SuperscriptBox[
          RowBox[{"Log", "[", 
           RowBox[{"1", "-", "z"}], "]"}], "2"], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"4", " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"-", "1"}], "+", "z"}], ")"}]}], "+", 
           RowBox[{
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}]}], ")"}]}], "+", 
        RowBox[{
         RowBox[{"Log", "[", 
          RowBox[{"1", "-", "z"}], "]"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{"-", "3"}], " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"-", "4"}], "+", 
              RowBox[{"3", " ", "z"}]}], ")"}]}], "-", 
           RowBox[{"8", " ", 
            RowBox[{"(", 
             RowBox[{"2", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}], "+", 
           RowBox[{"4", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            SuperscriptBox[
             RowBox[{"Log", "[", "z", "]"}], "2"]}]}], ")"}]}], "+", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{"6", " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"-", "1"}], "+", "z"}], ")"}]}], "-", 
           RowBox[{"2", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", 
             RowBox[{"1", "-", "z"}], "]"}]}], "-", 
           RowBox[{"4", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}]}], ")"}], " ", 
         RowBox[{"PolyLog", "[", 
          RowBox[{"2", ",", "z"}], "]"}]}], "-", 
        RowBox[{"6", " ", 
         RowBox[{"(", 
          RowBox[{"1", "+", "z"}], ")"}], " ", 
         RowBox[{"PolyLog", "[", 
          RowBox[{"3", ",", 
           RowBox[{"1", "-", "z"}]}], "]"}]}], "+", 
        RowBox[{"12", " ", 
         RowBox[{"(", 
          RowBox[{"1", "+", "z"}], ")"}], " ", 
         RowBox[{"PolyLog", "[", 
          RowBox[{"3", ",", "z"}], "]"}]}], "-", 
        RowBox[{"28", " ", 
         RowBox[{"(", 
          RowBox[{"1", "+", "z"}], ")"}], " ", 
         RowBox[{"Zeta", "[", "3", "]"}]}], "+", "\[IndentingNewLine]", 
        "\[IndentingNewLine]", 
        RowBox[{"16", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"D3", "[", "z", "]"}], "-", 
           RowBox[{
            RowBox[{"Zeta", "[", "2", "]"}], 
            RowBox[{"D1", "[", "z", "]"}]}], " ", "+", 
           RowBox[{"2", 
            RowBox[{"Zeta", "[", "3", "]"}], 
            RowBox[{"D0", "[", "z", "]"}]}]}], ")"}]}], "+", 
        "\[IndentingNewLine]", "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{"+", "4"}], 
            RowBox[{"(", 
             FractionBox[
              RowBox[{
               RowBox[{
                RowBox[{"-", "3"}], " ", 
                RowBox[{"Log", "[", "z", "]"}]}], "+", 
               RowBox[{"4", " ", 
                RowBox[{"Log", "[", 
                 RowBox[{"1", "-", "z"}], "]"}], " ", 
                RowBox[{"Log", "[", "z", "]"}]}], "+", 
               RowBox[{"2", " ", 
                RowBox[{"(", 
                 RowBox[{
                  RowBox[{"PolyLog", "[", 
                   RowBox[{"2", ",", "z"}], "]"}], "-", 
                  RowBox[{"Zeta", "[", "2", "]"}]}], ")"}]}]}], 
              RowBox[{"1", "-", "z"}]], ")"}]}], "-", 
           RowBox[{"8", 
            FractionBox[
             SuperscriptBox[
              RowBox[{"Log", "[", "z", "]"}], "2"], 
             RowBox[{"1", "-", "z"}]]}], "+", 
           RowBox[{"2", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{
               RowBox[{"-", "4"}], " ", 
               RowBox[{"(", 
                RowBox[{"2", "+", "z"}], ")"}]}], "+", 
              RowBox[{"12", " ", 
               RowBox[{"D0", "[", "z", "]"}]}], "+", 
              RowBox[{"8", " ", 
               RowBox[{"D1", "[", "z", "]"}]}], "-", 
              RowBox[{"4", " ", 
               RowBox[{"(", 
                RowBox[{"1", "+", "z"}], ")"}], " ", 
               RowBox[{"Log", "[", 
                RowBox[{"1", "-", "z"}], "]"}]}], "+", 
              RowBox[{"3", " ", 
               RowBox[{"(", 
                RowBox[{"1", "+", "z"}], ")"}], " ", 
               RowBox[{"Log", "[", "z", "]"}]}]}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}], "-", 
           RowBox[{"4", 
            RowBox[{"Zeta", "[", "2", "]"}], " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}]}], "+", 
           RowBox[{"2", " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"-", "6"}], "+", 
              RowBox[{"5", " ", "z"}]}], ")"}]}], "+", 
           RowBox[{"12", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            SuperscriptBox[
             RowBox[{"Log", "[", 
              RowBox[{"1", "-", "z"}], "]"}], "2"]}], "+", 
           RowBox[{"4", " ", 
            RowBox[{"(", 
             RowBox[{"2", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}], "-", 
           RowBox[{
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            SuperscriptBox[
             RowBox[{"Log", "[", "z", "]"}], "2"]}], "+", 
           RowBox[{
            RowBox[{"Log", "[", 
             RowBox[{"1", "-", "z"}], "]"}], " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"8", " ", 
               RowBox[{"(", 
                RowBox[{"2", "+", "z"}], ")"}]}], "-", 
              RowBox[{"8", " ", 
               RowBox[{"(", 
                RowBox[{"1", "+", "z"}], ")"}], " ", 
               RowBox[{"Log", "[", "z", "]"}]}]}], ")"}]}], "+", 
           "\[IndentingNewLine]", 
           RowBox[{"8", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{
               RowBox[{"Zeta", "[", "2", "]"}], " ", 
               RowBox[{"D0", "[", "z", "]"}]}], "-", 
              RowBox[{"3", 
               RowBox[{"D1", "[", "z", "]"}]}], "-", 
              RowBox[{"3", 
               RowBox[{"D2", "[", "z", "]"}]}]}], ")"}]}]}], ")"}], 
         RowBox[{"Log", "[", 
          FractionBox[
           SuperscriptBox["\[Mu]", "2"], "s"], "]"}]}], "+", 
        "\[IndentingNewLine]", "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{"-", "4"}], 
            RowBox[{"(", 
             RowBox[{
              FractionBox[
               RowBox[{"Log", "[", "z", "]"}], 
               RowBox[{"1", "-", "z"}]], "+", 
              RowBox[{"(", 
               RowBox[{"2", "+", "z"}], ")"}], "+", " ", 
              RowBox[{
               RowBox[{"(", 
                RowBox[{"1", "+", "z"}], ")"}], " ", 
               RowBox[{"Log", "[", 
                RowBox[{"1", "-", "z"}], "]"}]}]}], ")"}]}], "+", 
           RowBox[{"3", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}], "+", 
           RowBox[{"12", " ", 
            RowBox[{"D0", "[", "z", "]"}]}], "+", 
           RowBox[{"8", " ", 
            RowBox[{"D1", "[", "z", "]"}]}]}], ")"}], " ", 
         SuperscriptBox[
          RowBox[{"Log", "[", 
           FractionBox[
            SuperscriptBox["\[Mu]", "2"], "s"], "]"}], "2"]}]}], ")"}], 
      ")"}]}]}]}], ";"}]], "Input",
 CellChangeTimes->{{3.991962762485559*^9, 3.9919628157655897`*^9}, {
   3.9919628601422367`*^9, 3.991962860712298*^9}, 3.991963854167644*^9, 
   3.991963888761228*^9, {3.991964016920919*^9, 3.991964030856778*^9}, {
   3.991964859732211*^9, 3.991964922146896*^9}, {3.9919649727899714`*^9, 
   3.991964989724015*^9}, {3.991966298230769*^9, 3.99196629952005*^9}, 
   3.9919669326372423`*^9, 
   3.9919682590776033`*^9},ExpressionUUID->"5d52e1bb-a966-4a6a-8a87-\
f5d106daba3f"],

Cell[BoxData[
 RowBox[{"(*", 
  RowBox[{
   RowBox[{"For", " ", "help"}], ",", " ", 
   RowBox[{
   "i", " ", "implement", " ", "also", " ", "the", " ", "second", " ", "leg", 
    " ", "of", " ", "the", " ", "6", "th", " ", "class", " ", "of", " ", 
    "finite", " ", "terms"}]}], "*)"}]], "Input",
 CellChangeTimes->{{3.991968177704175*^9, 
  3.991968214166353*^9}},ExpressionUUID->"4662e19e-a0c9-405d-956f-\
8d75b4e3c800"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"FINITEDYp", "[", "60", "]"}], "=", " ", 
   RowBox[{
    RowBox[{"CF", " ", 
     RowBox[{"FLM", "[", 
      RowBox[{"q", ",", 
       RowBox[{"qp", " ", "z"}], ",", "l", ",", "lp", ",", "z"}], "]"}], 
     RowBox[{"(", 
      RowBox[{"1", "-", "z", "-", 
       RowBox[{"2", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "z"}], ")"}], " ", 
        RowBox[{"Log", "[", 
         RowBox[{"1", "-", "z"}], "]"}]}], "+", 
       RowBox[{"4", " ", 
        RowBox[{"D1", "[", "z", "]"}]}], "+", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{"1", "+", "z", "-", 
          RowBox[{"2", " ", 
           RowBox[{"D0", "[", "z", "]"}]}]}], ")"}], " ", 
        RowBox[{"Log", "[", 
         FractionBox[
          SuperscriptBox["\[Mu]", "2"], 
          RowBox[{"4", 
           SuperscriptBox["EC", "2"]}]], "]"}]}]}], " ", ")"}], 
     "leg2finalcalG"}], "+", 
    RowBox[{"CF", "  ", 
     RowBox[{"FLM", "[", 
      RowBox[{"q", ",", 
       RowBox[{"qp", " ", "z"}], ",", "l", ",", "lp", ",", "z"}], "]"}], 
     RowBox[{"(", 
      RowBox[{"1", "-", "z", "-", 
       RowBox[{"2", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "z"}], ")"}], " ", 
        RowBox[{"Log", "[", 
         RowBox[{"1", "-", "z"}], "]"}]}], "+", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{"1", "+", "z", "-", 
          RowBox[{"2", " ", 
           RowBox[{"D0", "[", "z", "]"}]}]}], ")"}], " ", 
        RowBox[{"Log", "[", "z", "]"}]}], "+", 
       RowBox[{"4", " ", 
        RowBox[{"D1", "[", "z", "]"}]}], "+", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{"1", "+", "z", "-", 
          RowBox[{"2", " ", 
           RowBox[{"D0", "[", "z", "]"}]}]}], ")"}], " ", 
        RowBox[{"Log", "[", 
         FractionBox[
          SuperscriptBox["\[Mu]", "2"], "s"], "]"}]}]}], ")"}], 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"-", "2"}], "Qlp", " ", "Qqp", " ", 
       RowBox[{"Log", "[", "z", "]"}], 
       RowBox[{"(", 
        RowBox[{"Log", "[", 
         RowBox[{"s24", "/", "s12"}], "]"}], ")"}]}], ")"}]}], "+", 
    RowBox[{"CF", "   ", 
     RowBox[{"FLM", "[", " ", 
      RowBox[{"q", ",", 
       RowBox[{"qp", " ", "z"}], ",", "l", ",", "lp", ",", "z"}], "]"}], 
     RowBox[{"Qqp", "^", "2"}], 
     RowBox[{"(", " ", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{
         FractionBox["1", 
          RowBox[{"1", "-", "z"}]], 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"8", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"4", " ", 
               RowBox[{"Zeta", "[", "2", "]"}]}], "+", 
              RowBox[{"PolyLog", "[", 
               RowBox[{"2", ",", "z"}], "]"}], "+", 
              RowBox[{
               FractionBox["3", "4"], " ", 
               SuperscriptBox[
                RowBox[{"Log", "[", 
                 RowBox[{"1", "-", "z"}], "]"}], "2"]}]}], ")"}], 
            RowBox[{"Log", "[", "z", "]"}]}], "+", 
           RowBox[{"10", " ", 
            SuperscriptBox[
             RowBox[{"Log", "[", "z", "]"}], "3"]}], "+", 
           RowBox[{"12", 
            RowBox[{"Log", "[", 
             RowBox[{"1", "-", "z"}], "]"}], " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"Zeta", "[", "2", "]"}], "-", " ", 
              RowBox[{"PolyLog", "[", 
               RowBox[{"2", ",", "z"}], "]"}], "+", 
              RowBox[{"2", 
               RowBox[{"Log", "[", "z", "]"}]}]}], ")"}]}], "-", 
           RowBox[{"4", " ", 
            RowBox[{"PolyLog", "[", 
             RowBox[{"3", ",", 
              RowBox[{"1", "-", "z"}]}], "]"}]}], "-", 
           RowBox[{"32", " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"PolyLog", "[", 
               RowBox[{"3", ",", "z"}], "]"}], "-", 
              RowBox[{"Zeta", "[", "3", "]"}]}], ")"}]}]}], ")"}]}], "+", 
        "\[IndentingNewLine]", 
        RowBox[{"4", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{"-", "3"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}], "+", 
           RowBox[{"4", " ", 
            RowBox[{"Log", "[", 
             RowBox[{"1", "-", "z"}], "]"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}], "+", 
           RowBox[{"2", " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"PolyLog", "[", 
               RowBox[{"2", ",", "z"}], "]"}], "-", 
              RowBox[{"Zeta", "[", "2", "]"}]}], ")"}]}]}], ")"}], 
         FractionBox[
          RowBox[{"Log", "[", "z", "]"}], 
          RowBox[{"1", "-", "z"}]]}], "-", 
        RowBox[{"4", 
         FractionBox[
          SuperscriptBox[
           RowBox[{"Log", "[", "z", "]"}], "3"], 
          RowBox[{"1", "-", "z"}]]}], "+", 
        RowBox[{"12", 
         RowBox[{"D0", "[", "z", "]"}], " ", 
         RowBox[{"(", 
          RowBox[{"1", "-", 
           RowBox[{
            FractionBox["4", "9"], 
            RowBox[{"Log", "[", "z", "]"}]}]}], ")"}], 
         SuperscriptBox[
          RowBox[{"Log", "[", "z", "]"}], "2"]}], "-", 
        RowBox[{"16", " ", 
         RowBox[{"D2", "[", "z", "]"}], " ", 
         RowBox[{"Log", "[", "z", "]"}]}], "+", "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{"-", 
             FractionBox["2", "3"]}], " ", 
            SuperscriptBox["\[Pi]", "2"], " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}]}], "+", 
           RowBox[{"2", " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"-", "6"}], "+", 
              RowBox[{"5", " ", "z"}]}], ")"}]}], "+", 
           RowBox[{"12", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            SuperscriptBox[
             RowBox[{"Log", "[", 
              RowBox[{"1", "-", "z"}], "]"}], "2"]}], "+", 
           RowBox[{"4", " ", 
            RowBox[{"(", 
             RowBox[{"2", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}], "-", 
           RowBox[{
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            SuperscriptBox[
             RowBox[{"Log", "[", "z", "]"}], "2"]}], "+", 
           RowBox[{
            RowBox[{"Log", "[", 
             RowBox[{"1", "-", "z"}], "]"}], " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"8", " ", 
               RowBox[{"(", 
                RowBox[{"2", "+", "z"}], ")"}]}], "-", 
              RowBox[{"8", " ", 
               RowBox[{"(", 
                RowBox[{"1", "+", "z"}], ")"}], " ", 
               RowBox[{"Log", "[", "z", "]"}]}]}], ")"}]}], "+", 
           RowBox[{
            FractionBox["4", "3"], " ", 
            SuperscriptBox["\[Pi]", "2"], " ", 
            RowBox[{"D0", "[", "z", "]"}]}], "-", 
           RowBox[{"24", 
            RowBox[{"D1", "[", "z", "]"}]}], "-", 
           RowBox[{"24", 
            RowBox[{"D2", "[", "z", "]"}]}]}], ")"}], 
         RowBox[{"Log", "[", "z", "]"}]}], "+", 
        RowBox[{
         SuperscriptBox[
          RowBox[{"Log", "[", "z", "]"}], "2"], 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{"-", "4"}], " ", 
            RowBox[{"(", 
             RowBox[{"2", "+", "z"}], ")"}]}], "-", 
           RowBox[{"4", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", 
             RowBox[{"1", "-", "z"}], "]"}]}], "+", 
           RowBox[{"3", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}]}], ")"}]}], "+", 
        "\[IndentingNewLine]", 
        RowBox[{"2", " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "9"}], "+", 
           RowBox[{"8", " ", "z"}]}], ")"}]}], "-", 
        RowBox[{"8", " ", 
         RowBox[{"(", 
          RowBox[{"1", "+", "z"}], ")"}], " ", 
         SuperscriptBox[
          RowBox[{"Log", "[", 
           RowBox[{"1", "-", "z"}], "]"}], "3"]}], "+", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{"5", "-", 
           RowBox[{"9", " ", "z"}]}], ")"}], " ", 
         RowBox[{"Log", "[", "z", "]"}]}], "+", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{"2", "-", "z"}], ")"}], " ", 
         SuperscriptBox[
          RowBox[{"Log", "[", "z", "]"}], "2"]}], "-", 
        RowBox[{
         FractionBox["13", "6"], " ", 
         RowBox[{"(", 
          RowBox[{"1", "+", "z"}], ")"}], " ", 
         SuperscriptBox[
          RowBox[{"Log", "[", "z", "]"}], "3"]}], "+", 
        RowBox[{
         SuperscriptBox["\[Pi]", "2"], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"3", " ", 
            RowBox[{"(", 
             RowBox[{"1", "-", "z"}], ")"}]}], "+", 
           RowBox[{
            FractionBox["5", "3"], " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", 
             RowBox[{"1", "-", "z"}], "]"}]}], "-", 
           RowBox[{
            FractionBox["8", "3"], " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}]}], ")"}]}], "+", 
        RowBox[{
         SuperscriptBox[
          RowBox[{"Log", "[", 
           RowBox[{"1", "-", "z"}], "]"}], "2"], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"4", " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"-", "1"}], "+", "z"}], ")"}]}], "+", 
           RowBox[{
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}]}], ")"}]}], "+", 
        RowBox[{
         RowBox[{"Log", "[", 
          RowBox[{"1", "-", "z"}], "]"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{"-", "3"}], " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"-", "4"}], "+", 
              RowBox[{"3", " ", "z"}]}], ")"}]}], "-", 
           RowBox[{"8", " ", 
            RowBox[{"(", 
             RowBox[{"2", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}], "+", 
           RowBox[{"4", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            SuperscriptBox[
             RowBox[{"Log", "[", "z", "]"}], "2"]}]}], ")"}]}], "+", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{"6", " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"-", "1"}], "+", "z"}], ")"}]}], "-", 
           RowBox[{"2", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", 
             RowBox[{"1", "-", "z"}], "]"}]}], "-", 
           RowBox[{"4", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}]}], ")"}], " ", 
         RowBox[{"PolyLog", "[", 
          RowBox[{"2", ",", "z"}], "]"}]}], "-", 
        RowBox[{"6", " ", 
         RowBox[{"(", 
          RowBox[{"1", "+", "z"}], ")"}], " ", 
         RowBox[{"PolyLog", "[", 
          RowBox[{"3", ",", 
           RowBox[{"1", "-", "z"}]}], "]"}]}], "+", 
        RowBox[{"12", " ", 
         RowBox[{"(", 
          RowBox[{"1", "+", "z"}], ")"}], " ", 
         RowBox[{"PolyLog", "[", 
          RowBox[{"3", ",", "z"}], "]"}]}], "-", 
        RowBox[{"28", " ", 
         RowBox[{"(", 
          RowBox[{"1", "+", "z"}], ")"}], " ", 
         RowBox[{"Zeta", "[", "3", "]"}]}], "+", "\[IndentingNewLine]", 
        "\[IndentingNewLine]", 
        RowBox[{"16", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"D3", "[", "z", "]"}], "-", 
           RowBox[{
            RowBox[{"Zeta", "[", "2", "]"}], 
            RowBox[{"D1", "[", "z", "]"}]}], " ", "+", 
           RowBox[{"2", 
            RowBox[{"Zeta", "[", "3", "]"}], 
            RowBox[{"D0", "[", "z", "]"}]}]}], ")"}]}], "+", 
        "\[IndentingNewLine]", "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{"+", "4"}], 
            RowBox[{"(", 
             FractionBox[
              RowBox[{
               RowBox[{
                RowBox[{"-", "3"}], " ", 
                RowBox[{"Log", "[", "z", "]"}]}], "+", 
               RowBox[{"4", " ", 
                RowBox[{"Log", "[", 
                 RowBox[{"1", "-", "z"}], "]"}], " ", 
                RowBox[{"Log", "[", "z", "]"}]}], "+", 
               RowBox[{"2", " ", 
                RowBox[{"(", 
                 RowBox[{
                  RowBox[{"PolyLog", "[", 
                   RowBox[{"2", ",", "z"}], "]"}], "-", 
                  RowBox[{"Zeta", "[", "2", "]"}]}], ")"}]}]}], 
              RowBox[{"1", "-", "z"}]], ")"}]}], "-", 
           RowBox[{"8", 
            FractionBox[
             SuperscriptBox[
              RowBox[{"Log", "[", "z", "]"}], "2"], 
             RowBox[{"1", "-", "z"}]]}], "+", 
           RowBox[{"2", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{
               RowBox[{"-", "4"}], " ", 
               RowBox[{"(", 
                RowBox[{"2", "+", "z"}], ")"}]}], "+", 
              RowBox[{"12", " ", 
               RowBox[{"D0", "[", "z", "]"}]}], "+", 
              RowBox[{"8", " ", 
               RowBox[{"D1", "[", "z", "]"}]}], "-", 
              RowBox[{"4", " ", 
               RowBox[{"(", 
                RowBox[{"1", "+", "z"}], ")"}], " ", 
               RowBox[{"Log", "[", 
                RowBox[{"1", "-", "z"}], "]"}]}], "+", 
              RowBox[{"3", " ", 
               RowBox[{"(", 
                RowBox[{"1", "+", "z"}], ")"}], " ", 
               RowBox[{"Log", "[", "z", "]"}]}]}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}], "-", 
           RowBox[{"4", 
            RowBox[{"Zeta", "[", "2", "]"}], " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}]}], "+", 
           RowBox[{"2", " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"-", "6"}], "+", 
              RowBox[{"5", " ", "z"}]}], ")"}]}], "+", 
           RowBox[{"12", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            SuperscriptBox[
             RowBox[{"Log", "[", 
              RowBox[{"1", "-", "z"}], "]"}], "2"]}], "+", 
           RowBox[{"4", " ", 
            RowBox[{"(", 
             RowBox[{"2", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}], "-", 
           RowBox[{
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            SuperscriptBox[
             RowBox[{"Log", "[", "z", "]"}], "2"]}], "+", 
           RowBox[{
            RowBox[{"Log", "[", 
             RowBox[{"1", "-", "z"}], "]"}], " ", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{"8", " ", 
               RowBox[{"(", 
                RowBox[{"2", "+", "z"}], ")"}]}], "-", 
              RowBox[{"8", " ", 
               RowBox[{"(", 
                RowBox[{"1", "+", "z"}], ")"}], " ", 
               RowBox[{"Log", "[", "z", "]"}]}]}], ")"}]}], "+", 
           "\[IndentingNewLine]", 
           RowBox[{"8", 
            RowBox[{"(", 
             RowBox[{
              RowBox[{
               RowBox[{"Zeta", "[", "2", "]"}], " ", 
               RowBox[{"D0", "[", "z", "]"}]}], "-", 
              RowBox[{"3", 
               RowBox[{"D1", "[", "z", "]"}]}], "-", 
              RowBox[{"3", 
               RowBox[{"D2", "[", "z", "]"}]}]}], ")"}]}]}], ")"}], 
         RowBox[{"Log", "[", 
          FractionBox[
           SuperscriptBox["\[Mu]", "2"], "s"], "]"}]}], "+", 
        "\[IndentingNewLine]", "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{"-", "4"}], 
            RowBox[{"(", 
             RowBox[{
              FractionBox[
               RowBox[{"Log", "[", "z", "]"}], 
               RowBox[{"1", "-", "z"}]], "+", 
              RowBox[{"(", 
               RowBox[{"2", "+", "z"}], ")"}], "+", " ", 
              RowBox[{
               RowBox[{"(", 
                RowBox[{"1", "+", "z"}], ")"}], " ", 
               RowBox[{"Log", "[", 
                RowBox[{"1", "-", "z"}], "]"}]}]}], ")"}]}], "+", 
           RowBox[{"3", " ", 
            RowBox[{"(", 
             RowBox[{"1", "+", "z"}], ")"}], " ", 
            RowBox[{"Log", "[", "z", "]"}]}], "+", 
           RowBox[{"12", " ", 
            RowBox[{"D0", "[", "z", "]"}]}], "+", 
           RowBox[{"8", " ", 
            RowBox[{"D1", "[", "z", "]"}]}]}], ")"}], " ", 
         SuperscriptBox[
          RowBox[{"Log", "[", 
           FractionBox[
            SuperscriptBox["\[Mu]", "2"], "s"], "]"}], "2"]}]}], ")"}], 
      ")"}]}]}]}], ";"}]], "Input",
 CellChangeTimes->{{3.991968207790625*^9, 3.991968312839275*^9}},
 CellLabel->
  "In[748]:=",ExpressionUUID->"a7a40bee-fba5-4706-a8fb-44bc1e5c622b"],

Cell[BoxData[
 RowBox[{"(*", 
  RowBox[{"Elastic", " ", "bit"}], "*)"}]], "Input",
 CellChangeTimes->{{3.9919628299425*^9, 
  3.991962833106448*^9}},ExpressionUUID->"124777c3-157f-406b-87f4-\
d89758f9ab4f"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"FINITEDYp", "[", "7", "]"}], "=", " ", 
   RowBox[{
    RowBox[{"CF", " ", 
     RowBox[{"FLM", "[", 
      RowBox[{"q", ",", "qbp", ",", "nu", ",", "l"}], "]"}], 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"4", 
        RowBox[{"Zeta", "[", "2", "]"}]}], "-", 
       RowBox[{"3", 
        RowBox[{"Log", "[", 
         FractionBox[
          SuperscriptBox["\[Mu]", "2"], 
          RowBox[{"4", " ", 
           SuperscriptBox["EC", "2"]}]], "]"}]}]}], ")"}], "genGnloQCD"}], "+",
     "\[IndentingNewLine]", " ", 
    RowBox[{"CF", " ", 
     RowBox[{"FLM", "[", 
      RowBox[{"q", ",", "qp", ",", "l", ",", "lp"}], "]"}], 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        RowBox[{"-", 
         FractionBox["2", "45"]}], 
        RowBox[{"(", 
         RowBox[{
          SuperscriptBox["Qq", "2"], "+", 
          SuperscriptBox["Qqp", "2"]}], ")"}], 
        SuperscriptBox["\[Pi]", "4"]}], "+", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           RowBox[{"-", "6"}], 
           SuperscriptBox[
            RowBox[{"(", 
             RowBox[{"Qlp", "-", "Qq", "-", "Qqp"}], ")"}], "2"]}], "-", 
          RowBox[{"16", 
           RowBox[{"(", 
            RowBox[{
             SuperscriptBox["Qq", "2"], "+", 
             SuperscriptBox["Qqp", "2"]}], ")"}], 
           RowBox[{"Log", "[", 
            FractionBox[
             SuperscriptBox["\[Mu]", "2"], 
             RowBox[{"4", " ", 
              SuperscriptBox["EC", "2"]}]], "]"}]}]}], ")"}], 
        RowBox[{"Zeta", "[", "3", "]"}]}], "+", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           RowBox[{"-", "9"}], 
           SuperscriptBox["Qlp", "2"]}], "+", 
          RowBox[{
           RowBox[{"(", 
            RowBox[{
             RowBox[{"-", "9"}], "-", 
             FractionBox[
              RowBox[{"8", " ", 
               SuperscriptBox["\[Pi]", "2"]}], "3"]}], ")"}], " ", 
           SuperscriptBox["Qq", "2"]}], "-", 
          RowBox[{"18", "Qq", " ", "Qqp"}], "-", 
          RowBox[{
           FractionBox["1", "3"], " ", 
           RowBox[{"(", 
            RowBox[{"27", "+", 
             RowBox[{"8", " ", 
              SuperscriptBox["\[Pi]", "2"]}]}], ")"}], " ", 
           SuperscriptBox["Qqp", "2"]}], "+", 
          RowBox[{"18", "Qlp", " ", 
           RowBox[{"(", 
            RowBox[{"Qq", "+", "Qqp"}], ")"}]}]}], ")"}], 
        FractionBox["1", "4"], 
        SuperscriptBox[
         RowBox[{"(", 
          RowBox[{"Log", "[", 
           FractionBox[
            SuperscriptBox["\[Mu]", "2"], 
            RowBox[{"4", " ", 
             SuperscriptBox["EC", "2"]}]], "]"}], ")"}], "2"]}], "-", 
       RowBox[{
        SuperscriptBox["\[Pi]", "2"], " ", 
        SuperscriptBox[
         RowBox[{"(", 
          RowBox[{"Qlp", "-", "Qq", "-", "Qqp"}], ")"}], "2"], " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", 
           FractionBox["21", "8"]}], "-", " ", 
          RowBox[{"Log", "[", 
           FractionBox[
            SuperscriptBox["\[Mu]", "2"], 
            RowBox[{"4", " ", 
             SuperscriptBox["EC", "2"]}]], "]"}]}], ")"}]}]}], ")"}]}]}]}], 
  ";"}]], "Input",
 CellChangeTimes->{{3.991962850580179*^9, 3.991962873307253*^9}, 
   3.991965129100658*^9, 3.991965171831678*^9, {3.991966322965094*^9, 
   3.991966338641944*^9}, {3.991966380284135*^9, 3.991966380475478*^9}, 
   3.991966422600358*^9},
 CellLabel->
  "In[747]:=",ExpressionUUID->"721d599e-8ca3-4293-b7b4-8bace57c4d23"]
},
WindowSize->{809, 949},
WindowMargins->{{Automatic, -850}, {Automatic, 16}},
FrontEndVersion->"14.3 for Mac OS X ARM (64-bit) (July 8, 2025)",
StyleDefinitions->"Default.nb",
ExpressionUUID->"8d19c329-992d-4bc5-924a-bc072f23f581"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[554, 20, 1516, 39, 111, "Input",ExpressionUUID->"b30319e1-2140-427c-aa44-b44a4fc51180"],
Cell[2073, 61, 213, 5, 29, "Input",ExpressionUUID->"f236d87a-10c6-4ab3-8c22-437a2885c0f1"],
Cell[2289, 68, 3313, 99, 203, "Input",ExpressionUUID->"d60e41d1-3362-4c88-bfd0-0fb736a2ac8d"],
Cell[5605, 169, 3919, 107, 247, "Input",ExpressionUUID->"712dedb1-0bdc-42cb-bb10-9d500316a9b5"],
Cell[9527, 278, 2530, 68, 192, "Input",ExpressionUUID->"e93843a7-a8e6-4d19-b073-2630ae096939"],
Cell[12060, 348, 236, 5, 29, "Input",ExpressionUUID->"275029be-23b0-4453-a5ca-ef94fbc4167c"],
Cell[12299, 355, 1576, 46, 67, "Input",ExpressionUUID->"20fd7f69-dcb4-4ffa-9864-3140d26158cd"],
Cell[13878, 403, 3111, 84, 252, "Input",ExpressionUUID->"c7cac081-5f2c-426e-8998-5bddfd52e026"],
Cell[16992, 489, 184, 4, 29, "Input",ExpressionUUID->"3059f3cd-9fd6-4ecd-8b35-6d8bbb42698f"],
Cell[17179, 495, 1410, 41, 65, "Input",ExpressionUUID->"85b59acb-5a88-43b6-b790-3d4ef61b48b8"],
Cell[18592, 538, 628, 19, 48, "Input",ExpressionUUID->"a30eb47c-ddbe-487e-b668-222da2ec68b2"],
Cell[19223, 559, 211, 5, 29, "Input",ExpressionUUID->"894ed87b-f744-4e51-a547-41a2bae4dcac"],
Cell[19437, 566, 5538, 152, 520, "Input",ExpressionUUID->"836b927f-6ea5-4518-859d-47c1b6089cad"],
Cell[24978, 720, 589, 12, 70, "Input",ExpressionUUID->"b88ff76d-0dfe-450b-89e9-9e20dd8cc655"],
Cell[25570, 734, 17518, 481, 889, "Input",ExpressionUUID->"5d52e1bb-a966-4a6a-8a87-f5d106daba3f"],
Cell[43091, 1217, 425, 10, 29, "Input",ExpressionUUID->"4662e19e-a0c9-405d-956f-8d75b4e3c800"],
Cell[43519, 1229, 17159, 473, 919, "Input",ExpressionUUID->"a7a40bee-fba5-4706-a8fb-44bc1e5c622b"],
Cell[60681, 1704, 206, 5, 29, "Input",ExpressionUUID->"124777c3-157f-406b-87f4-d89758f9ab4f"],
Cell[60890, 1711, 3582, 104, 203, "Input",ExpressionUUID->"721d599e-8ca3-4293-b7b4-8bace57c4d23"]
}
]
*)

