(* Content-type: application/vnd.wolfram.mathematica *)

(*** Wolfram Notebook File ***)
(* http://www.wolfram.com/nb *)

(* CreatedBy='Wolfram 14.3' *)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[       154,          7]
NotebookDataLength[     39564,       1112]
NotebookOptionsPosition[     37879,       1084]
NotebookOutlinePosition[     38274,       1100]
CellTagsIndexPosition[     38231,       1097]
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
     " ", "of", " ", "the", " ", "finite", " ", "remainder", " ", "of", " ", 
     "the", " ", "subtraction", " ", "for", " ", "q", " ", "qx", " ", 
     "channel", " ", "in", " ", "Neutral", " ", "Drell"}], "-", 
    RowBox[{"Yan", " ", 
     RowBox[{"(", "DY0", ")"}], " ", "from", " ", "Federico", " ", 
     RowBox[{"Buccioni", "."}]}]}], " ", "*)"}], 
  RowBox[{"(*", " ", 
   RowBox[{"Convention", ":", " ", 
    RowBox[{
     RowBox[{"func", "[", 
      RowBox[{
       RowBox[{"p1", ".", "zi"}], ",", 
       RowBox[{"p2", ".", "zj"}], ",", "...", ",", "zi", ",", "zj"}], "]"}], " ",
      "means", " ", 
     RowBox[{"1", "/", "zi"}], "*", 
     RowBox[{"1", "/", "zj"}], "*", 
     RowBox[{"func", "[", 
      RowBox[{
       RowBox[{"p1", ".", "zi"}], ",", 
       RowBox[{"p2", ".", "zj"}], ",", "..."}], "]"}]}]}], " ", 
   "*)"}]}]], "Input",
 CellChangeTimes->{{3.991962890183772*^9, 3.991962959541946*^9}, {
  3.991963018509359*^9, 3.991963077879485*^9}, {3.991963144482834*^9, 
  3.991963153990328*^9}},ExpressionUUID->"b30319e1-2140-427c-aa44-\
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
   RowBox[{"Ptildeqq", "[", 
    RowBox[{"z_", ",", "E1_"}], "]"}], ":=", 
   RowBox[{
    RowBox[{"4", " ", 
     RowBox[{"D1", "[", "z", "]"}]}], "-", 
    RowBox[{"2", 
     RowBox[{"(", 
      RowBox[{"1", "+", "z"}], ")"}], 
     RowBox[{"Log", "[", 
      RowBox[{"1", "-", "z"}], "]"}]}], "+", 
    RowBox[{"(", 
     RowBox[{"1", "-", "z"}], ")"}], "+", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{
       RowBox[{"2", 
        RowBox[{"D0", "[", "z", "]"}]}], "-", 
       RowBox[{"(", 
        RowBox[{"1", "+", "z"}], ")"}]}], ")"}], "2", 
     RowBox[{"Log", "[", 
      FractionBox[
       RowBox[{"2", " ", "E1"}], "\[Mu]"], "]"}]}]}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"PqgAP0", "[", "z_", "]"}], ":=", 
   RowBox[{
    SuperscriptBox[
     RowBox[{"(", 
      RowBox[{"1", "-", "z"}], ")"}], "2"], "+", 
    SuperscriptBox["z", "2"]}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Gq", "[", 
   RowBox[{"Qq_", ",", "Qe_"}], "]"}], ":=", 
  RowBox[{"(", 
   RowBox[{
    RowBox[{
     SuperscriptBox["Qe", "2"], 
     RowBox[{"(", 
      RowBox[{"13", "-", 
       RowBox[{"4", 
        RowBox[{"Zeta", "[", "2", "]"}]}], "+", 
       SuperscriptBox[
        RowBox[{"Log", "[", 
         FractionBox["E3", "E4"], "]"}], "2"], "+", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{"3", "-", 
          RowBox[{"2", 
           RowBox[{"Log", "[", 
            FractionBox[
             RowBox[{"E3", " ", "E4"}], 
             SuperscriptBox["EC", "2"]], "]"}]}]}], ")"}], 
        RowBox[{"Log", "[", 
         RowBox[{"eta", "[", 
          RowBox[{"3", ",", "4"}], "]"}], "]"}]}], " ", "+", 
       RowBox[{"2", " ", 
        RowBox[{"PolyLog", "[", 
         RowBox[{"2", ",", 
          RowBox[{"1", "-", 
           RowBox[{"eta", "[", 
            RowBox[{"3", ",", "4"}], "]"}]}]}], "]"}]}]}], ")"}]}], "+", 
    "\[IndentingNewLine]", 
    RowBox[{"2", " ", "Qe", " ", "Qq", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        RowBox[{"(", 
         RowBox[{
          FractionBox["3", "2"], "+", 
          RowBox[{"Log", "[", 
           FractionBox["EC", "E4"], "]"}]}], ")"}], " ", 
        RowBox[{"Log", "[", 
         FractionBox[
          RowBox[{"eta", "[", 
           RowBox[{"2", ",", "4"}], "]"}], 
          RowBox[{"eta", "[", 
           RowBox[{"1", ",", "4"}], "]"}]], "]"}]}], "+", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{
          FractionBox["3", "2"], "+", 
          RowBox[{"Log", "[", 
           FractionBox["EC", "E3"], "]"}]}], ")"}], " ", 
        RowBox[{"Log", "[", 
         FractionBox[
          RowBox[{"eta", "[", 
           RowBox[{"1", ",", "3"}], "]"}], 
          RowBox[{"eta", "[", 
           RowBox[{"2", ",", "3"}], "]"}]], "]"}]}], "+", 
       RowBox[{"PolyLog", "[", 
        RowBox[{"2", ",", 
         RowBox[{"1", "-", 
          RowBox[{"eta", "[", 
           RowBox[{"1", ",", "3"}], "]"}]}]}], "]"}], "-", 
       RowBox[{"PolyLog", "[", 
        RowBox[{"2", ",", 
         RowBox[{"1", "-", 
          RowBox[{"eta", "[", 
           RowBox[{"1", ",", "4"}], "]"}]}]}], "]"}], "-", 
       RowBox[{"PolyLog", "[", 
        RowBox[{"2", ",", 
         RowBox[{"1", "-", 
          RowBox[{"eta", "[", 
           RowBox[{"2", ",", "3"}], "]"}]}]}], "]"}], "+", 
       RowBox[{"PolyLog", "[", 
        RowBox[{"2", ",", 
         RowBox[{"1", "-", 
          RowBox[{"eta", "[", 
           RowBox[{"2", ",", "4"}], "]"}]}]}], "]"}]}], ")"}]}], "+", 
    RowBox[{
     SuperscriptBox["Qq", "2"], "4", 
     RowBox[{"Zeta", "[", "2", "]"}]}]}], ")"}]}]}], "Input",
 CellChangeTimes->{{3.9919624301534967`*^9, 3.991962430277067*^9}, {
  3.9919632915726843`*^9, 
  3.991963302612094*^9}},ExpressionUUID->"d60e41d1-3362-4c88-bfd0-\
0fb736a2ac8d"],

Cell[BoxData[
 RowBox[{"(*", "FLVQCDfin", "*)"}]], "Input",
 CellChangeTimes->{{3.9919624470469646`*^9, 3.991962449922947*^9}, {
  3.991962562051062*^9, 
  3.991962562289609*^9}},ExpressionUUID->"275029be-23b0-4453-a5ca-\
ef94fbc4167c"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"FINITEDY0", "[", "1", "]"}], "=", 
   RowBox[{
    SuperscriptBox["Qq", "2"], "  ", 
    RowBox[{"Identity", "[", 
     RowBox[{
      RowBox[{"FLVQCDfin", "[", 
       RowBox[{"p1", ",", 
        RowBox[{"z", ".", "p2"}], ",", "p3", ",", "p4", ",", "z"}], "]"}], "+", 
      RowBox[{"FLVQCDfin", "[", 
       RowBox[{
        RowBox[{"z", ".", "p1"}], ",", "p2", ",", "p3", ",", "p4", ",", "z"}],
        "]"}]}], "]"}], " ", 
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
         RowBox[{"Log", "[", "z", "]"}]}], ")"}]}]}], ")"}]}]}], 
  ";"}]], "Input",
 CellChangeTimes->{{3.991962470281146*^9, 3.99196254999155*^9}, {
  3.991962591636293*^9, 3.99196259215806*^9}, {3.991963384951714*^9, 
  3.9919633896285954`*^9}},ExpressionUUID->"20fd7f69-dcb4-4ffa-9864-\
3140d26158cd"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"FINITEDY0", "[", "2", "]"}], "=", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{
      RowBox[{"(", 
       RowBox[{"13", "-", 
        FractionBox[
         RowBox[{"2", " ", 
          SuperscriptBox["\[Pi]", "2"]}], "3"]}], ")"}], " ", 
      SuperscriptBox["Qe", "2"], " ", 
      RowBox[{"FLVQCDfin", "[", 
       RowBox[{"p1", ",", "p2", ",", "p3", ",", "p4"}], "]"}]}], "+", 
     RowBox[{
      FractionBox["2", "3"], " ", 
      SuperscriptBox["\[Pi]", "2"], " ", 
      SuperscriptBox["Qq", "2"], " ", 
      RowBox[{"FLVQCDfin", "[", 
       RowBox[{"p1", ",", "p2", ",", "p3", ",", "p4"}], "]"}]}], "-", 
     RowBox[{"3", " ", 
      SuperscriptBox["Qq", "2"], " ", 
      RowBox[{"FLVQCDfin", "[", 
       RowBox[{"p1", ",", "p2", ",", "p3", ",", "p4"}], "]"}], " ", 
      RowBox[{"Log", "[", 
       FractionBox["mu2", "s"], "]"}]}], "+", 
     RowBox[{"Qe", " ", "Qq", " ", 
      RowBox[{"FLVQCDfin", "[", 
       RowBox[{"p1", ",", "p2", ",", "p3", ",", "p4"}], "]"}], " ", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"3", " ", 
         RowBox[{"Log", "[", 
          FractionBox[
           RowBox[{"s13", " ", "s24"}], 
           RowBox[{"s14", " ", "s23"}]], "]"}]}], "+", 
        RowBox[{"4", " ", 
         RowBox[{"PolyLog", "[", 
          RowBox[{"2", ",", 
           RowBox[{"1", "-", "eta13"}]}], "]"}]}], "-", 
        RowBox[{"4", " ", 
         RowBox[{"PolyLog", "[", 
          RowBox[{"2", ",", 
           RowBox[{"1", "-", "eta14"}]}], "]"}]}]}], ")"}]}]}], ")"}]}], 
  ";"}]], "Input",
 CellChangeTimes->{{3.991962514245973*^9, 3.991962548417899*^9}, {
  3.991962589577559*^9, 
  3.9919625899964848`*^9}},ExpressionUUID->"c7cac081-5f2c-426e-8998-\
5bddfd52e026"],

Cell[BoxData[
 RowBox[{"(*", "FLVEWfin", "*)"}]], "Input",
 CellChangeTimes->{{3.991962556051277*^9, 
  3.991962567675682*^9}},ExpressionUUID->"3059f3cd-9fd6-4ecd-8b35-\
6d8bbb42698f"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"FINITEDY0", "[", "3", "]"}], "=", 
   RowBox[{"CF", "  ", 
    RowBox[{"Identity", "[", 
     RowBox[{
      RowBox[{"FLVQEDfin", "[", 
       RowBox[{"p1", ",", 
        RowBox[{"z", ".", "p2"}], ",", "p3", ",", "p4", ",", "z"}], "]"}], "+", 
      RowBox[{"FLVQEDfin", "[", 
       RowBox[{
        RowBox[{"z", ".", "p1"}], ",", "p2", ",", "p3", ",", "p4", ",", "z"}],
        "]"}]}], "]"}], " ", 
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
       RowBox[{"Log", "[", "z", "]"}]}]}], ")"}]}]}], ";"}]], "Input",
 CellChangeTimes->{{3.991962583429284*^9, 3.991962616875204*^9}, {
   3.991962708156433*^9, 3.991962712548923*^9}, 
   3.991963624842144*^9},ExpressionUUID->"85b59acb-5a88-43b6-b790-\
3d4ef61b48b8"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"FINITEDY0", "[", "4", "]"}], "=", 
   RowBox[{"CF", 
    RowBox[{"(", 
     RowBox[{
      FractionBox[
       RowBox[{"2", " ", 
        SuperscriptBox["\[Pi]", "2"]}], "3"], "-", 
      RowBox[{"3", " ", 
       RowBox[{"Log", "[", 
        FractionBox["mu2", "s"], "]"}]}]}], ")"}], 
    RowBox[{"FLVQEDfin", "[", 
     RowBox[{"p1", ",", "p2", ",", "p3", ",", "p4"}], "]"}]}]}], 
  ";"}]], "Input",
 CellChangeTimes->{{3.991962617932856*^9, 3.991962664248831*^9}, {
  3.991963762051001*^9, 
  3.991963766082056*^9}},ExpressionUUID->"a30eb47c-ddbe-487e-b668-\
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
   RowBox[{"2", "CF", " ", 
    SuperscriptBox["Qq", "2"], " ", 
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
  3.99196362838428*^9, 
  3.991963629717799*^9}},ExpressionUUID->"836b927f-6ea5-4518-859d-\
47c1b6089cad"],

Cell[BoxData[
 RowBox[{"(*", 
  RowBox[{"Single", " ", 
   RowBox[{"boosted", ":", " ", 
    RowBox[{
    "derived", " ", "for", " ", "leg", " ", "1", " ", "but", " ", "similar", " ",
      "for", " ", "leg", " ", "2", " ", "with", " ", "careful", " ", 
     "replacement", " ", "of", " ", "index", " ", "numbers"}]}]}], 
  "*)"}]], "Input",
 CellChangeTimes->{{3.99196272486541*^9, 3.991962727143035*^9}, {
  3.99196276696769*^9, 
  3.9919627896524*^9}},ExpressionUUID->"b88ff76d-0dfe-450b-89e9-9e20dd8cc655"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"FINITEDY0", "[", "6", "]"}], "=", 
   RowBox[{
    RowBox[{"CF", " ", 
     RowBox[{"Identity", "[", 
      RowBox[{"FLM", "[", 
       RowBox[{
        RowBox[{"z", ".", "p1"}], ",", "p2", ",", "p3", ",", "p4", ",", "z"}],
        "]"}], "]"}], 
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
     RowBox[{"Gq", "[", 
      RowBox[{"Qq", ",", "Qe"}], "]"}]}], "+", "\[IndentingNewLine]", 
    RowBox[{"CF", " ", "\[Alpha]", " ", "\[Alpha]s", " ", 
     RowBox[{"IdentityMatrix", "[", 
      RowBox[{"FLM", "[", 
       RowBox[{
        RowBox[{"z", ".", "p1"}], ",", "p2", ",", "p3", ",", "p4", ",", "z"}],
        "]"}], "]"}], " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"-", "2"}], " ", "Qe", " ", "Qq", " ", 
       RowBox[{"(", 
        RowBox[{
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
         RowBox[{"Log", "[", "z", "]"}], 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"0", " ", 
            RowBox[{"Log", "[", 
             FractionBox["s13", "s14"], "]"}]}], "+", 
           RowBox[{"1", " ", 
            RowBox[{"Log", "[", 
             FractionBox["s24", "s23"], "]"}]}]}], ")"}]}], ")"}]}], ")"}]}], 
    "+", "\[IndentingNewLine]", 
    RowBox[{"CF", " ", 
     RowBox[{"Identity", "[", 
      RowBox[{"FLM", "[", 
       RowBox[{
        RowBox[{"z", ".", "p1"}], ",", "p2", ",", "p3", ",", "p4", ",", "z"}],
        "]"}], "]"}], " ", 
     SuperscriptBox["Qq", "2"], 
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
  3.9919628601422367`*^9, 3.991962860712298*^9}, {3.991963634023703*^9, 
  3.991963638391954*^9}, {3.9919645400186577`*^9, 
  3.99196454141101*^9}},ExpressionUUID->"5d52e1bb-a966-4a6a-8a87-\
f5d106daba3f"],

Cell[BoxData[
 RowBox[{"(*", 
  RowBox[{"Elastic", " ", "bit"}], "*)"}]], "Input",
 CellChangeTimes->{{3.9919628299425*^9, 
  3.991962833106448*^9}},ExpressionUUID->"124777c3-157f-406b-87f4-\
d89758f9ab4f"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"FINITEDY0", "[", "7", "]"}], "=", " ", 
   RowBox[{
    RowBox[{"CF", " ", 
     RowBox[{"FLM", "[", 
      RowBox[{"p1", ",", "p2", ",", "p3", ",", "p4"}], "]"}], 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"4", 
        RowBox[{"Zeta", "[", "2", "]"}]}], "-", 
       RowBox[{"3", 
        RowBox[{"Log", "[", 
         FractionBox[
          SuperscriptBox["\[Mu]", "2"], 
          RowBox[{"4", " ", 
           SuperscriptBox["EC", "2"]}]], "]"}]}]}], ")"}], 
     RowBox[{"Gq", "[", 
      RowBox[{"Qq", ",", "Qe"}], "]"}]}], "+", "\[IndentingNewLine]", 
    "\[IndentingNewLine]", " ", 
    RowBox[{"CF", " ", 
     RowBox[{"FLM", "[", 
      RowBox[{"p1", ",", "p2", ",", "p3", ",", "p4"}], "]"}], " ", 
     SuperscriptBox["Qq", "2"], 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        RowBox[{"-", 
         FractionBox["4", "45"]}], " ", 
        SuperscriptBox["\[Pi]", "4"]}], "-", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{
          RowBox[{"2", 
           SuperscriptBox["\[Pi]", "2"]}], "+", 
          RowBox[{"32", " ", 
           RowBox[{"Zeta", "[", "3", "]"}]}]}], ")"}], " ", 
        RowBox[{"Log", "[", 
         FractionBox[
          SuperscriptBox["\[Mu]", "2"], 
          RowBox[{"4", " ", 
           SuperscriptBox["EC", "2"]}]], "]"}]}], " ", "+", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{"9", "-", 
          RowBox[{
           FractionBox["4", "3"], " ", 
           SuperscriptBox["\[Pi]", "2"]}]}], ")"}], " ", 
        SuperscriptBox[
         RowBox[{"Log", "[", 
          FractionBox[
           SuperscriptBox["\[Mu]", "2"], 
           RowBox[{"4", " ", 
            SuperscriptBox["EC", "2"]}]], "]"}], "2"]}]}], ")"}]}]}]}], 
  ";"}]], "Input",
 CellChangeTimes->{{3.991962850580179*^9, 3.991962873307253*^9}},
 CellLabel->
  "In[284]:=",ExpressionUUID->"721d599e-8ca3-4293-b7b4-8bace57c4d23"]
},
WindowSize->{808, 762},
WindowMargins->{{Automatic, -1681}, {Automatic, 16}},
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
Cell[554, 20, 1252, 31, 90, "Input",ExpressionUUID->"b30319e1-2140-427c-aa44-b44a4fc51180"],
Cell[1809, 53, 213, 5, 29, "Input",ExpressionUUID->"f236d87a-10c6-4ab3-8c22-437a2885c0f1"],
Cell[2025, 60, 3834, 118, 258, "Input",ExpressionUUID->"d60e41d1-3362-4c88-bfd0-0fb736a2ac8d"],
Cell[5862, 180, 236, 5, 29, "Input",ExpressionUUID->"275029be-23b0-4453-a5ca-ef94fbc4167c"],
Cell[6101, 187, 1521, 44, 65, "Input",ExpressionUUID->"20fd7f69-dcb4-4ffa-9864-3140d26158cd"],
Cell[7625, 233, 1758, 49, 172, "Input",ExpressionUUID->"c7cac081-5f2c-426e-8998-5bddfd52e026"],
Cell[9386, 284, 184, 4, 29, "Input",ExpressionUUID->"3059f3cd-9fd6-4ecd-8b35-6d8bbb42698f"],
Cell[9573, 290, 1453, 41, 65, "Input",ExpressionUUID->"85b59acb-5a88-43b6-b790-3d4ef61b48b8"],
Cell[11029, 333, 624, 19, 48, "Input",ExpressionUUID->"a30eb47c-ddbe-487e-b668-222da2ec68b2"],
Cell[11656, 354, 211, 5, 29, "Input",ExpressionUUID->"894ed87b-f744-4e51-a547-41a2bae4dcac"],
Cell[11870, 361, 5454, 149, 520, "Input",ExpressionUUID->"836b927f-6ea5-4518-859d-47c1b6089cad"],
Cell[17327, 512, 510, 11, 70, "Input",ExpressionUUID->"b88ff76d-0dfe-450b-89e9-9e20dd8cc655"],
Cell[17840, 525, 17892, 490, 991, "Input",ExpressionUUID->"5d52e1bb-a966-4a6a-8a87-f5d106daba3f"],
Cell[35735, 1017, 206, 5, 29, "Input",ExpressionUUID->"124777c3-157f-406b-87f4-d89758f9ab4f"],
Cell[35944, 1024, 1931, 58, 113, "Input",ExpressionUUID->"721d599e-8ca3-4293-b7b4-8bace57c4d23"]
}
]
*)

