#!/usr/bin/env wolframscript
(* Headless driver for the DY+ finite-remainder cross-checks.

   It evaluates the input cells of the DY+_qq.wl notebook (so that genGnloQCD,
   FINITEDYp[7], leg1finalcalG, leg2finalcalG are defined), then runs the
   symbolic comparisons in validate_DYp7.wl.

   Usage (from anywhere):
       wolframscript -file /path/to/math/run_validate.wls
*)

dir = DirectoryName[$InputFileName];
If[dir === "" || dir === $Failed, dir = Directory[]];

nb = Get[FileNameJoin[{dir, "DY+_qq.wl"}]];
inputs = Cases[nb, Cell[BoxData[b_], "Input", ___] :> b, Infinity];
Scan[Function[b, Quiet[ToExpression[b, StandardForm]]], inputs];

Get[FileNameJoin[{dir, "validate_DYp7.wl"}]];
