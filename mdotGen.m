(* ::Package:: *)

#!/usr/local/bin/MathematicaScript/ -script


Needs["vwEff`"]


SetDirectory["/Users/aleksey/Second_Year_Project/star_formation"];
(*Reading saved up quantities from file*)
etaf=Import["eta.csv", "HeaderLines"->1];
(*Black hole masses*)
mbhs=etaf[[;;,1]]*MS;
(*etas*)
\[Eta]s=etaf[[;;,2]];
(*Total vwEff*)
vweffTots=Import["vwSources.csv", "HeaderLines"->1][[;;,-1]]*10.^5;
(*Ia radius*)
radiiIa=Import["Ia.csv", "HeaderLines"->1][[;;,2]];


mdots = Table[mdotsol[mbhs[[i]], vweffTots[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}]; 
(*mdotCores = Table[mdotsol[mbhs[[i]], vweffTotsCore[[i]], 0.1, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}]; *)

mdotEdds = mdotEdd /@ mbhs; 
LEdds=LEdd/@mbhs;


mdotIasCore = Table[mdotIA[mbhs[[i]], radiiIa[[i]], 0.1, \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}];
mdotIas = Table[mdotIA[mbhs[[i]], radiiIa[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}];
mdotInfs = (\[Eta]s*mbhs)/th; 

mdotMaxCools = Table[eddrMaxCool[mbhs[[i]], 0.1, \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}]*mdotEdds;
mdotMaxCoolsCore = Table[eddrMaxCool[mbhs[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}]*mdotEdds;

(*rhoCusps=Table[rhoRs[mbhs[[i]], vweffTots[[i]], 0.8, etas[[i]]],{i, 1, Length[mbhs]}]
  rhoCores=Table[rhoRs[mbhs[[i]], vweffTots[[i]], 0.1, etas[[i]]],{i, 1, Length[mbhs]}]*)


hcs = Table[hc[mbhs[[i]], vweffTots[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}];
(*hcsCore = Table[hc[mbhs[[i]], vweffTotsCore[[i]], 0.1, \[Eta]s[[i]]], 
  {i, 1, Length[mbhs]}]; *)


Export["hc.csv", Transpose[{mbhs/MS, hcs}], "TableHeadings" -> {"Mbh", "Cusp"}]

Export["mdots.csv", Transpose[{mbhs/MS, mdots, (*mdotCores,*) mdotIas, mdotIasCore, mdotInfs, mdotMaxCools,\
mdotMaxCoolsCore, LEdds}],TableHeadings->{"Mbh","Cusp","Core","Ia","IaCore","Inf",\
"MaxCools","MaxCoolsCore","LEdds"}]




