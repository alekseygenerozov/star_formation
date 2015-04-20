#!/usr/local/bin/MathematicaScript/ -script


AppendTo[$Path,"/Users/aleksey/code/mathematica/star_formation_code"];
Needs["vwEff`"]


mhalos=10.^Range[10.8, 14.,0.2] MS;
mstarTots=mstarTot/@mhalos;
mbhs=MbhMbulge/@mstarTots;
\[CapitalGamma]s=\[CapitalGamma]fitM/@mbhs;


(*heating sources*)
tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], 0.8], {i, 1, Length[mhalos]}];
vweffStars=tmp[[1 ;; All,1]]; 
vweffMSPs=tmp[[1 ;; All,2]]; 
vwc = tmp[[1 ;; All,3]]; 
vweffIas = tmp[[1 ;; All,4]]; 
vweffTots = tmp[[1 ;; All,5]]; 
vweffIasCorrected= tmp[[1;;All,6]];

tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], 0.1], {i, 1, Length[mhalos]}];
vweffStars=tmp[[1 ;; All,1]]; 
vweffMSPs=tmp[[1 ;; All,2]]; 
vwcCore = tmp[[1 ;; All,3]]; 
vweffIas = tmp[[1 ;; All,4]]; 
vweffTotsCore = tmp[[1 ;; All,5]];  
vweffIasCoreCorrected= tmp[[1;;All,6]];

tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], \[CapitalGamma]s[[i]]], {i, 1, Length[mhalos]}];
vweffStars=tmp[[1 ;; All,1]]; 
vweffMSPs=tmp[[1 ;; All,2]]; 
vwcGamma = tmp[[1 ;; All,3]]; 
vweffIas = tmp[[1 ;; All,4]]; 
vweffTotsGamma = tmp[[1 ;; All,5]]; 
vweffIasGammaCorrected= tmp[[1;;All,6]];


vweffTots-vweffIasCorrected


(*Ia radii*)
radiiIa = Table[radiusIa[mbhs[[i]], mhalos[[i]]], {i, 1, Length[mbhs]}]; 


(*etas...*)
\[Eta]s = \[Eta]/@mhalos;


SetDirectory["/Users/aleksey/Second_Year_Project/star_formation"]
Export["vwSources.csv", Transpose[{mbhs*(10.^5/MS), vweffIas, vweffIasCorrected, vweffStars, vweffMSPs, vwc, vweffTots}/10.^5], 
  "TableHeadings" -> {"Mbh", "Ias", "IasCorrected", "Stars", "MSPs", "Compton", "Total"}]
Export["vwSourcesCore.csv", Transpose[{mbhs*(10.^5/MS), vwcCore, vweffIasCoreCorrected, vweffTotsCore}/10.^5], 
  "TableHeadings" -> {"mbh", "Compton", "IasCorrected", "Total"}]
Export["vwSourcesGamma.csv", Transpose[{mbhs*(10.^5/MS), vwcGamma, vweffIasGammaCorrected, vweffTotsGamma}/10.^5], 
  "TableHeadings" -> {"mbh", "Compton", "IasCorrected", "Total"}]

Export["Ia.csv", Transpose[{mbhs/MS, radiiIa}], "TableHeadings" -> {"Mbh", "rIa", "rsCusp", "rsCore", "rsGamma"}]
Export["eta.csv", Transpose[{mbhs/MS, \[Eta]s}], "TableHeadings" -> {"Mbh", "eta"}]

