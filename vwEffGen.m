(* ::Package:: *)

#!/usr/local/bin/MathematicaScript/ -script


AppendTo[$Path,"/Users/aleksey/code/mathematica/star_formation_code"];
Needs["vwEff`"]


mhalos=10.^Range[10.8, 14.,0.2] MS;
mstarTots=mstarTot/@mhalos;
mbhs=MbhMbulge/@mstarTots;
\[CapitalGamma]s=\[CapitalGamma]fitM/@mbhs;


tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], 0.8], {i, 1, Length[mhalos]}];
vweffStars=tmp[[1 ;; All,1]]; 
vweffMSPs=tmp[[1 ;; All,2]]; 
vwc = tmp[[1 ;; All,3]]; 
vweffIas = tmp[[1 ;; All,4]]; 
vweffTots = tmp[[1 ;; All,5]]; 
vweffIasCorrected= tmp[[1;;All,6]];

tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], 0.1], {i, 1, Length[mhalos]}];
vweffStars2=tmp[[1 ;; All,1]]; 
vweffMSPs=tmp[[1 ;; All,2]]; 
vwc = tmp[[1 ;; All,3]]; 
vweffIas = tmp[[1 ;; All,4]]; 
vweffTots = tmp[[1 ;; All,5]];  
vweffIasCorrected= tmp[[1;;All,6]];

tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], \[CapitalGamma]s[[i]]], {i, 1, Length[mhalos]}];
vweffStars=tmp[[1 ;; All,1]]; 
vweffMSPs=tmp[[1 ;; All,2]]; 
vwc = tmp[[1 ;; All,3]]; 
vweffIas = tmp[[1 ;; All,4]]; 
vweffTots = tmp[[1 ;; All,5]]; 
vweffIasCorrected= tmp[[1;;All,6]];


(*(*Impulsive limit mbh=10^6 msun...*)
mbh=10.^6 MS;
timeImps=10.^Range[6.6,10.,0.2]*year;
\[Eta]Imps=\[Eta]Imp/@timeImps;
tmp=vweffTotImp[mbh, #,0.8]&/@timeImps;
(*heating from different sources in the impulsive limit*)
vwcImps=tmp[[;;,1]];
vweffIaImps=tmp[[;;, 2]];
vweffTotImps=tmp[[;;, 3]];
vweffIaImps2=tmp[[;;, 4]];

vweffStarImps=vweffStarImp/@timeImps;
vweffMSPImps=vweffMSPImp/@timeImps;
vweffNoBHImps=Sqrt[vweffTotImps^2-vwcImps^2];
(*Heating over cooling ratio*)
hcImps=Table[hc[mbh,vweffNoBHImps[[i]], 0.8, \[Eta]Imps[[i]]], {i, 1, Length[timeImps]}];
hcIaImps=Table[hcIa[mbh,vweffIaImps2[[i]], 0.8, \[Eta]Imps[[i]]], {i, 1, Length[timeImps]}];

Transpose[{timeImps/year*10.^5, vweffIaImps, vweffIaImps2, vweffStarImps, vweffMSPImps, vwcImps, vweffNoBHImps}/10.^5]//Export["/Users/aleksey/Second_Year_Project/star_formation/vwSourcesImps6.csv", #, TableHeadings->{"Time","Ias", "IasCorrected", "Stars", "MSPs", "Compton", "Total"}]&;
Transpose[{timeImps/year, hcImps, hcIaImps}]//Export["/Users/aleksey/Second_Year_Project/star_formation/hcImps6.csv", #, TableHeadings->{"Time", "Cusp", "Ia"}]&;*)


(*(*Impulsive limit mbh=10^8 msun...*)
mbh=10.^8 MS;
timeImps=10.^Range[6.6,10.,0.2]*year;
\[Eta]Imps=\[Eta]Imp/@timeImps;
tmp=vweffTotImp[mbh, #,0.8]&/@timeImps;
(*heating from different sources in the impulsive limit*)
vwcImps=tmp[[;;,1]];
vweffIaImps=tmp[[;;, 2]];
vweffTotImps=tmp[[;;, 3]];
vweffIaImps2=tmp[[;;, 4]];

vweffStarImps=vweffStarImp/@timeImps;
vweffMSPImps=vweffMSPImp/@timeImps;
vweffNoBHImps=Sqrt[vweffTotImps^2-vwcImps^2];
(*Heating over cooling ratio*)
hcImps=Table[hc[mbh,vweffNoBHImps[[i]], 0.8, \[Eta]Imps[[i]]], {i, 1, Length[timeImps]}];
hcIaImps=Table[hcIa[mbh,vweffIaImps2[[i]], 0.8, \[Eta]Imps[[i]]], {i, 1, Length[timeImps]}];

Transpose[{timeImps*10.^5/year, vweffIaImps, vweffIaImps2,  vweffStarImps, vweffMSPImps, vwcImps, vweffNoBHImps}/10.^5]//Export["/Users/aleksey/Second_Year_Project/star_formation/vwSourcesImps8.csv", #, TableHeadings->{"Time","Ias", "IasCorrected", "Stars", "MSPs", "Compton", "Total"}]&;
Transpose[{timeImps/year, hcImps, hcIaImps}]//Export["/Users/aleksey/Second_Year_Project/star_formation/hcImps8.csv", #, TableHeadings->{"Time", "Cusp", "Ia"}]&;*)


SetDirectory["/Users/aleksey/Second_Year_Project/star_formation"]
Export["vwSources.csv", Transpose[{mbhs*(10.^5/MS), vweffIas, vweffStars, vweffMSPs, vwc, vweffTots}/10.^5], 
  "TableHeadings" -> {"Mbh", "Ias", "IasCorrected", "Stars", "MSPs", "Compton", "Total"}]
Export["vwSourcesCore.csv", Transpose[{mbhs*(10.^5/MS), vwcCore, vweffIasCoreCorrected, vweffTotsCore}/10.^5], 
  "TableHeadings" -> {"mbh", "Compton", "IasCorrected", "Total"}]
Export["vwSourcesGamma.csv", Transpose[{mbhs*(10.^5/MS), vwcGamma, vweffIasGammaCorrected, vweffTotsGamma}/10.^5], 
  "TableHeadings" -> {"mbh", "Compton", "IasCorrected", "Total"}]

Export["Ia.csv", Transpose[{mbhs/MS, radiiIa, rsCusp, rsCore, rsGamma}], "TableHeadings" -> {"Mbh", "rIa", "rsCusp", "rsCore", "rsGamma"}]
Export["eta.csv", Transpose[{mbhs/MS, \[Eta]s}], "TableHeadings" -> {"Mbh", "eta"}]
Export["etaImp.csv", Transpose[{timeImps/year, \[Eta]Imps}], "TableHeadings" -> {"time", "eta"}]
Export["vweffStarImps.csv", Transpose[{timeImps/year, vweffStarImps}], "TableHeadings" -> {"time", "vw"}]
