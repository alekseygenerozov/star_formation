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


mdotMaxCoolsCore = Table[mdotMaxCool[mbhs[[i]], 0.1, \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}];
mdotMaxCools = Table[mdotMaxCool[mbhs[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}];
mdotMaxCoolsGamma = Table[mdotMaxCool[mbhs[[i]], \[CapitalGamma]s[[i]], \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}];


TempC=10.^9;

mdotComptonsCore = Table[mdotCompton[mbhs[[i]], 0.1, \[Eta]s[[i]], TempC], 
    {i, 1, Length[\[Eta]s]}];
mdotComptons = Table[mdotCompton[mbhs[[i]], 0.8, \[Eta]s[[i]], TempC], 
    {i, 1, Length[\[Eta]s]}];
mdotComptonsGamma = Table[mdotCompton[mbhs[[i]], \[CapitalGamma]s[[i]], \[Eta]s[[i]],  TempC], 
    {i, 1, Length[\[Eta]s]}];


hcs = Table[hc[mbhs[[i]], vweffTots[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}];
hcsCore = Table[hc[mbhs[[i]], vweffTotsCore[[i]], 0.1, \[Eta]s[[i]]], 
  {i, 1, Length[mbhs]}]; 


Needs["SigFig`"]

Export["hc.csv", Transpose[{mbhs/MS, hcs, hcsCore}], "TableHeadings" -> {"Mbh", "Cusp","Core"}]

mdotsAll=Map[OutputForm[SigForm[#,3, scientific->True]]&\
,Transpose[{mbhs/MS, mdots, mdotsCore, mdotsGamma, mdotIas, mdotIasCore,\
	    mdotIasGamma, mdotMaxCools, mdotMaxCoolsCore, mdotMaxCoolsGamma, mdotComptons, mdotComptonsCore, mdotComptonsGamma, LEdds}],{2}];

Export["mdots.csv", mdotsAll, TableHeadings->
{"Mbh","Cusp","Core","Gamma",\
 "Ia","IaCore","IaGamma","MaxCools","MaxCoolsCore","MaxCoolsGamma", "mdotComptons", "mdotComptonsCore", "mdotComptonsGamma", "LEdds"}];

