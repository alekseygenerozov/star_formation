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
vweffIasCorrected=Import["vwSources.csv", "HeaderLines"->1][[;;,3]]*10.^5;
vweffOther=(vweffTots^2.-vweffIasCorrected^2.)^0.5;

vweffTotsCore=Import["vwSourcesCore.csv", "HeaderLines"->1][[;;,-1]]*10.^5;
vweffIasCoreCorrected=Import["vwSourcesCore.csv", "HeaderLines"->1][[;;,3]]*10.^5;
vweffOtherCore=(vweffTotsCore^2.-vweffIasCoreCorrected^2.)^0.5;

vweffTotsGamma=Import["vwSourcesGamma.csv", "HeaderLines"->1][[;;,-1]]*10.^5;
vweffIasGammaCorrected=Import["vwSourcesGamma.csv", "HeaderLines"->1][[;;,3]]*10.^5;
(*Ia radius*)
radiiIa=Import["Ia.csv", "HeaderLines"->1][[;;,2]];


rsCusp = Table[rs[mbhs[[i]], Sqrt[vweffTots[[i]]^2-vweffIasCorrected[[i]]^2], 0.8], {i, 1, Length[mbhs]}];
rsCore = Table[rs[mbhs[[i]], Sqrt[vweffTotsCore[[i]]^2-vweffIasCoreCorrected[[i]]^2], 0.1], {i, 1, Length[mbhs]}];
rsGamma = Table[rs[mbhs[[i]], Sqrt[vweffTotsGamma[[i]]^2-vweffIasGammaCorrected[[i]]^2], 0.1], {i, 1, Length[mbhs]}];

rsCusp2 = Table[rs[mbhs[[i]], vweffTots[[i]], 0.8], {i, 1, Length[mbhs]}];
rsCore2 = Table[rs[mbhs[[i]], vweffTotsCore[[i]], 0.1], {i, 1, Length[mbhs]}];
rsGamma2 = Table[rs[mbhs[[i]], vweffTotsGamma[[i]], 0.1], {i, 1, Length[mbhs]}];


hcs = Table[hc[mbhs[[i]], vweffTots[[i]], 0.8, \[Eta]s[[i]]], {i, 1, Length[mbhs]}];
hcsCore = Table[hc[mbhs[[i]], vweffTotsCore[[i]], 0.1, \[Eta]s[[i]]], 
  {i, 1, Length[mbhs]}]; 

(*Calculating heating/cooling ratio using Ia formalism. Using total heating--
implicitly assuming heating is dominated by Ia's. Only take this when rs=rIa*)
hcIas=Table[hcIa[mbhs[[i]], vweffTots[[i]], radiiIa[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}];
hcIasCore=Table[hcIa[mbhs[[i]], vweffTots[[i]], radiiIa[[i]], 0.1, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}];
hcOverall=Table[If[(vweffIasCorrected[[i]]>vweffOther[[i]])&&(rsCusp2[[i]]<radiiIa[[i]]), hcIas[[i]], hcs[[i]]], {i, 1, Length[hcs]}];
hcOverallCore=Table[If[(vweffIasCoreCorrected[[i]]>vweffOtherCore[[i]])&&(rsCore2[[i]]<radiiIa[[i]]), hcIasCore[[i]], hcs[[i]]], {i, 1, Length[hcs]}];


rsOverall=Table[If[rsCusp[[i]]>radiiIa[[i]], radiiIa[[i]], rsCusp[[i]]], {i, 1, Length[rsCusp]}];
rsOverallCore=Table[If[rsCore[[i]]>radiiIa[[i]], radiiIa[[i]], rsCore[[i]]], {i, 1, Length[rsCusp]}];
rsOverallGamma=Table[If[rsCusp[[i]]>radiiIa[[i]], radiiIa[[i]], rsGamma[[i]]], {i, 1, Length[rsCusp]}];


\[CapitalGamma]s=\[CapitalGamma]fitM/@mbhs;


mdotEdds = mdotEdd /@ mbhs; 
LEdds=LEdd/@mbhs;


mdotIasCore = Table[mdotIA[mbhs[[i]], radiiIa[[i]], 0.1, \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}];
mdotIas = Table[mdotIA[mbhs[[i]], radiiIa[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}];
mdotIasGamma = Table[mdotIA[mbhs[[i]], radiiIa[[i]], \[CapitalGamma]s[[i]], \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}];


mdotMaxCoolsCore = Table[mdotMaxCool[mbhs[[i]], 0.1, \[Eta]s[[i]]], {i, 1, Length[\[Eta]s]}];
mdotMaxCools = Table[mdotMaxCool[mbhs[[i]], 0.8, \[Eta]s[[i]]], {i, 1, Length[\[Eta]s]}];
mdotMaxCoolsGamma = Table[mdotMaxCool[mbhs[[i]], \[CapitalGamma]s[[i]], \[Eta]s[[i]]],{i, 1, Length[\[Eta]s]}];


TempC=10.^9;

mdotComptonsCore = Table[mdotCompton[mbhs[[i]], 0.1, \[Eta]s[[i]], TempC], 
    {i, 1, Length[\[Eta]s]}];
mdotComptons = Table[mdotCompton[mbhs[[i]], 0.8, \[Eta]s[[i]], TempC], 
    {i, 1, Length[\[Eta]s]}];
mdotComptonsGamma = Table[mdotCompton[mbhs[[i]], \[CapitalGamma]s[[i]], \[Eta]s[[i]],  TempC], 
    {i, 1, Length[\[Eta]s]}];


mdotsCore = Table[mdotsol[mbhs[[i]], vweffTotsCore[[i]], 0.1, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}]; 
mdots = Table[mdotsol[mbhs[[i]], vweffTots[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}]; 
mdotsGamma = Table[mdotsol[mbhs[[i]], vweffTotsGamma[[i]], \[CapitalGamma]s[[i]], \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}]; 


mdotsOverallCore =  Table[If[rsCore[[i]]<radiiIa[[i]], mdotsCore[[i]], mdotIasCore[[i]]], 
    {i, 1, Length[mbhs]}];
mdotsOverall =  Table[If[rsCusp[[i]]<radiiIa[[i]], mdots[[i]], mdotIas[[i]]], 
    {i, 1, Length[mbhs]}];  
mdotsOverallGamma = Table[If[rsGamma[[i]]<radiiIa[[i]], mdotsGamma[[i]], mdotIasGamma[[i]]], 
    {i, 1, Length[mbhs]}]; 


Needs["SigFig`"]
Export["hc.csv", Transpose[{mbhs/MS, hcs, hcsCore, hcIas, hcIasCore, hcOverall, hcOverallCore}], "TableHeadings" -> {"Mbh", "Cusp","Core", "Ia", "IaCore", "Overall", "OverallCore"}]

mdotsAll=Map[OutputForm[SigForm[#,3, scientific->True]]&\
,Transpose[{mbhs/MS, mdots, mdotsCore, mdotsGamma, mdotIas, mdotIasCore,\
mdotIasGamma, mdotMaxCools, mdotMaxCoolsCore, mdotMaxCoolsGamma, mdotComptons, mdotComptonsCore, mdotComptonsGamma,\
mdotsOverall, mdotsOverallCore, mdotsOverallGamma,  LEdds}],{2}];

Export["mdots.csv", mdotsAll, TableHeadings->
{"Mbh","Cusp","Core","Gamma",\
 "Ia","IaCore","IaGamma","MaxCools","MaxCoolsCore","MaxCoolsGamma", "mdotComptons", "mdotComptonsCore", "mdotComptonsGamma", "overall", "overallCore", "overallGamma", "LEdds"}];

Export["rs.csv", Transpose[{mbhs/MS, rsOverall, rsOverallCore,rsOverallGamma,radiiIa}], "TableHeadings"->
{"Mbh", "rsCusp", "rsCore", "rsGamma", "rIa"}]



