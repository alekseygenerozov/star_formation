(* ::Package:: *)

#!/usr/local/bin/MathematicaScript/ -script


AppendTo[$Path,"/Users/aleksey/code/mathematica/star_formation_code"];
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
vweffTotsCore=Import["vwSourcesCore.csv", "HeaderLines"->1][[;;,-1]]*10.^5;
vweffTotsGamma=Import["vwSourcesGamma.csv", "HeaderLines"->1][[;;,-1]]*10.^5;
(*Ia radius*)
radiiIa=Import["Ia.csv", "HeaderLines"->1][[;;,2]];

\[CapitalGamma]s=\[CapitalGamma]fitM/@mbhs;


mdots = Table[mdotsol[mbhs[[i]], vweffTots[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}]; 
mdotsCore = Table[mdotsol[mbhs[[i]], vweffTotsCore[[i]], 0.1, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}]; 
mdotsGamma = Table[mdotsol[mbhs[[i]], vweffTotsGamma[[i]], \[CapitalGamma]s[[i]], \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}]; 

mdotEdds = mdotEdd /@ mbhs; 
LEdds=LEdd/@mbhs;


mdotIasCore = Table[mdotIA[mbhs[[i]], radiiIa[[i]], 0.1, \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}];
mdotIas = Table[mdotIA[mbhs[[i]], radiiIa[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}];
mdotIasGamma = Table[mdotIA[mbhs[[i]], radiiIa[[i]], \[CapitalGamma]s[[i]], \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}];

mdotInfs = (\[Eta]s*mbhs)/th; 

mdotMaxCoolsCore = Table[eddrMaxCool[mbhs[[i]], 0.1, \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}]*mdotEdds;
mdotMaxCools = Table[eddrMaxCool[mbhs[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}]*mdotEdds;
mdotMaxCoolsGamma = Table[eddrMaxCool[mbhs[[i]], \[CapitalGamma]s[[i]], \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}]*mdotEdds;

(*rhoCusps=Table[rhoRs[mbhs[[i]], vweffTots[[i]], 0.8, etas[[i]]],{i, 1, Length[mbhs]}]
  rhoCores=Table[rhoRs[mbhs[[i]], vweffTots[[i]], 0.1, etas[[i]]],{i, 1, Length[mbhs]}]*)


hcs = Table[hc[mbhs[[i]], vweffTots[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}];
(*hcsCore = Table[hc[mbhs[[i]], vweffTotsCore[[i]], 0.1, \[Eta]s[[i]]], 
  {i, 1, Length[mbhs]}]; *)


Needs["SigFig`"]

Export["hc.csv", Transpose[{mbhs/MS, hcs}], "TableHeadings" -> {"Mbh", "Cusp"}]

mdotsAll=Map[OutputForm[SigForm[#,3, scientific->True]]&\
,Transpose[{mbhs/MS, mdots, mdotsCore, mdotsGamma, mdotIas, mdotIasCore,\
 mdotIasGamma, mdotInfs, mdotMaxCools, mdotMaxCoolsCore, mdotMaxCoolsGamma, LEdds}],2]

Export["mdots.csv", mdotsAll, TableHeadings->
{"Mbh","Cusp","Core","Gamma",\
 "Ia","IaCore","IaGamma","Inf","MaxCools","MaxCoolsCore","MaxCoolsGamma","LEdds"}]




