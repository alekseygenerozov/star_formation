(* ::Package:: *)

#!/usr/local/bin/MathematicaScript/ -script


AppendTo[$Path,"/Users/aleksey/code/mathematica/star_formation_code"];
Needs["vwEffBumpy`"]


gamma=0.1
mhalos=10.^Range[10.8,14.,0.2] MS;
mstarTots=mstarTot/@mhalos;
mbhs=MbhMbulge/@mstarTots;
Timing[tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], gamma], {i, 1, Length[mhalos]}]]
vweffTotsStandard= tmp[[1 ;; All,5]]; 


(* tScale = 10^8, \[Epsilon]=0.001 *)
vwEffBumpy`\[Epsilon]Floor=0.001;
vwEffBumpy`tScale=10.^8year;

mhalos=10.^Range[10.8,14.,0.2] MS;
mstarTots=mstarTot/@mhalos;
mbhs=MbhMbulge/@mstarTots;
\[CapitalGamma]s=\[CapitalGamma]fitM/@mbhs;
Timing[tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], gamma], {i, 1, Length[mhalos]}]]

vweffTotsA = tmp[[1 ;; All,5]]; 


(* tScale = 10^7, \[Epsilon]=0.001 *)
vwEffBumpy`\[Epsilon]Floor=0.001;
vwEffBumpy`tScale=10.^7year;

mhalos=10.^Range[10.8,14.,0.2] MS;
mstarTots=mstarTot/@mhalos;
mbhs=MbhMbulge/@mstarTots;
\[CapitalGamma]s=\[CapitalGamma]fitM/@mbhs;
Timing[tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], gamma], {i, 1, Length[mhalos]}]]

vweffTotsB = tmp[[1 ;; All,5]]; 


(* tScale = 10^6, \[Epsilon]=0.001 *)
vwEffBumpy`\[Epsilon]Floor=0.001;
vwEffBumpy`tScale=10.^6year;

mhalos=10.^Range[10.8,14.,0.2] MS;
mstarTots=mstarTot/@mhalos;
mbhs=MbhMbulge/@mstarTots;
\[CapitalGamma]s=\[CapitalGamma]fitM/@mbhs;

Timing[tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], gamma], {i, 1, Length[mhalos]}]]
vweffTotsC= tmp[[1 ;; All,5]]; 


(* tScale = 10^7, \[Epsilon]=0.01 *)
vwEffBumpy`\[Epsilon]Floor=0.01;
vwEffBumpy`tScale=10.^7year;

mhalos=10.^Range[10.8,14.,0.2] MS;
mstarTots=mstarTot/@mhalos;
mbhs=MbhMbulge/@mstarTots;
\[CapitalGamma]s=\[CapitalGamma]fitM/@mbhs;

Timing[tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], gamma], {i, 1, Length[mhalos]}]]
vweffTotsD= tmp[[1 ;; All,5]]; 


(* tScale = 10^7, \[Epsilon]=0.1 *)
vwEffBumpy`\[Epsilon]Floor=0.1;
vwEffBumpy`tScale=10.^7year;

mhalos=10.^Range[10.8,14.,0.2] MS;
mstarTots=mstarTot/@mhalos;
mbhs=MbhMbulge/@mstarTots;
\[CapitalGamma]s=\[CapitalGamma]fitM/@mbhs;

Timing[tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], gamma], {i, 1, Length[mhalos]}]]
vweffTotsE= tmp[[1 ;; All,5]]; 


(* tScale = 10^7, \[Epsilon]=0.3 *)
vwEffBumpy`\[Epsilon]Floor=0.3;
vwEffBumpy`tScale=10.^7year;

mhalos=10.^Range[10.8,14.,0.2] MS;
mstarTots=mstarTot/@mhalos;
mbhs=MbhMbulge/@mstarTots;
\[CapitalGamma]s=\[CapitalGamma]fitM/@mbhs;

Timing[tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], gamma], {i, 1, Length[mhalos]}]]
vweffTotsF= tmp[[1 ;; All,5]]; 


(* tScale = 10^8, \[Epsilon]=0.1 *)
vwEffBumpy`\[Epsilon]Floor=0.1;
vwEffBumpy`tScale=10.^8year;

mhalos=10.^Range[10.8,14.,0.2] MS;
mstarTots=mstarTot/@mhalos;
mbhs=MbhMbulge/@mstarTots;
\[CapitalGamma]s=\[CapitalGamma]fitM/@mbhs;

Timing[tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], gamma], {i, 1, Length[mhalos]}]]
vweffTotsG= tmp[[1 ;; All,5]]; 


(*ListLogLogPlot[{Transpose[{mbhs/MS, vweffTotsStandard}], Transpose[{mbhs/MS, vweffTotsA}], Transpose[{mbhs/MS, vweffTotsB}], Transpose[{mbhs/MS, vweffTotsC}], Transpose[{mbhs/MS, vweffTotsE}], Transpose[{mbhs/MS, vweffTotsF}], Transpose[{mbhs/MS, vweffTotsG}]}, \
Joined\[Rule]True,PlotStyle\[Rule]{Black, {Red, Dashed}, {Blue, Dashed}, {Purple, Dashed}, {Blue, Dotted}, {Blue, DotDashed}, {Red, Dotted}}, Frame\[Rule]True, FrameStyle\[Rule]Directive[Thick, 16], AxesStyle\[Rule]Directive[Thick, 16], FrameLabel\[Rule]{Style["Subscript[M, BH] [Subscript[M, Sun]]", 20], Style["Subscript[v, w] [cm/s]", 20]}]*)


dat=Transpose[{mbhs/MS, vweffTotsStandard, vweffTotsA, vweffTotsB, vweffTotsC, vweffTotsE, vweffTotsF, vweffTotsG}]
Export["/Users/aleksey/Second_Year_Project/star_formation/bumpy.csv",dat, TableHeadings->{"M", "Stand", "A", "B", "C", "E", "F", "G"}];
