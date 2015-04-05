#!/usr/local/bin/MathematicaScript/ -script


AppendTo[$Path,"/Users/aleksey/code/mathematica/star_formation_code"];
Needs["vwEff`"]


mhalos=10.^Range[10.8, 14.,0.2] MS;
mstarTots=mstarTot/@mhalos;
mbhs=MbhMbulge/@mstarTots;
\[CapitalGamma]s=\[CapitalGamma]fitM/@mbhs;


Transpose[{mbhs/MS,mhalos/MS}]


vweffStars = vweffStar /@ mhalos; 
vweffMSPs = vweffMSP /@ mhalos; 
vweffIas = vweffIa /@ mhalos;


tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], 0.8], {i, 1, Length[mhalos]}]
vwc = tmp[[1 ;; All,1]]; 
vweffTots = tmp[[1 ;; All,3]]; 

tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], 0.1], {i, 1, Length[mhalos]}]
vwcCore = tmp[[1 ;; All,1]]; 
vweffTotsCore = tmp[[1 ;; All,3]];

tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], \[CapitalGamma]s[[i]]], {i, 1, Length[mhalos]}]
vwcGamma = tmp[[1 ;; All,1]]; 
vweffTotsGamma = tmp[[1 ;; All,3]];


\[Eta]s = \[Eta] /@ mhalos; 


radiiIa = Table[radiusIa[mbhs[[i]], mhalos[[i]]], {i, 1, Length[mbhs]}]; 


rsCore = Table[rs[mbhs[[i]], vweffTots[[i]], 0.1], {i, 1, Length[mbhs]}];
rsCusp = Table[rs[mbhs[[i]], vweffTots[[i]], 0.8], {i, 1, Length[mbhs]}];


(*Impulsive limit mbh=10^6 msun...*)
mbh=10.^6 MS;
timeImps=10.^Range[6.6,10.,0.2]*year;
\[Eta]Imps=\[Eta]Imp/@timeImps;
tmp=vweffTotImp[mbh, #,0.8]&/@timeImps;
(*heating from different sources in the impulsive limit*)
vwcImps=tmp[[;;,1]];
vweffIaImps=tmp[[;;, 2]];
vweffTotImps=tmp[[;;, 3]];
vweffStarImps=vweffStarImp/@timeImps;
vweffMSPImps=vweffMSPImp/@timeImps;
vweffNoBHImps=Sqrt[vweffTotImps^2-vwcImps^2];
(*Heating over cooling ratio*)
hcImps=Table[hc[mbh,vweffNoBHImps[[i]], 0.8, \[Eta]Imps[[i]]], {i, 1, Length[timeImps]}];

Transpose[{timeImps/year*10.^5, vweffIaImps, vweffStarImps, vweffMSPImps, vwcImps, vweffNoBHImps}/10.^5]//Export["/Users/aleksey/Second_Year_Project/star_formation/vwSourcesImps6.csv", #, TableHeadings->{"Time","Ias", "Stars", "MSPs", "Compton", "Total"}]&;
Transpose[{timeImps/year, hcImps}]//Export["/Users/aleksey/Second_Year_Project/star_formation/hcImps6.csv", #, TableHeadings->{"Time", "Cusp"}]&;


(*Impulsive limit mbh=10^8 msun...*)
mbh=10.^8 MS;
timeImps=10.^Range[6.6,10.,0.2]*year;
\[Eta]Imps=\[Eta]Imp/@timeImps;
tmp=vweffTotImp[mbh, #,0.8]&/@timeImps;
(*heating from different sources in the impulsive limit*)
vwcImps=tmp[[;;,1]];
vweffIaImps=tmp[[;;, 2]];
vweffTotImps=tmp[[;;, 3]];
vweffStarImps=vweffStarImp/@timeImps;
vweffMSPImps=vweffMSPImp/@timeImps;
vweffNoBHImps=Sqrt[vweffTotImps^2-vwcImps^2];
(*Heating over cooling ratio*)
hcImps=Table[hc[mbh,vweffNoBHImps[[i]], 0.8, \[Eta]Imps[[i]]], {i, 1, Length[timeImps]}]

Transpose[{timeImps*10.^5/year, vweffIaImps, vweffStarImps, vweffMSPImps, vwcImps, vweffNoBHImps}/10.^5]//Export["/Users/aleksey/Second_Year_Project/star_formation/vwSourcesImps8.csv", #, TableHeadings->{"Time","Ias", "Stars", "MSPs", "Compton", "Total"}]&;
Transpose[{timeImps/year, hcImps}]//Export["/Users/aleksey/Second_Year_Project/star_formation/hcImps8.csv", #, TableHeadings->{"Time", "Cusp"}]&;


SetDirectory["/Users/aleksey/Second_Year_Project/star_formation"]
Export["vwSources.csv", Transpose[{mbhs*(10.^5/MS), vweffIas, vweffStars, vweffMSPs, vwc, vweffTots}/10.^5], 
  "TableHeadings" -> {"Mbh", "Ias", "Stars", "MSPs", "Compton", "Total"}]
Export["vwSourcesCore.csv", Transpose[{mbhs*(10.^5/MS), vwcCore, vweffTotsCore}/10.^5], 
  "TableHeadings" -> {"IaCorrected", "Compton", "Total"}]
Export["vwSourcesGamma.csv", Transpose[{mbhs*(10.^5/MS), vwcGamma, vweffTotsGamma}/10.^5], 
  "TableHeadings" -> {"Compton", "Total"}]

Export["Ia.csv", Transpose[{mbhs/MS, radiiIa, rsCusp, rsCore}], "TableHeadings" -> {"Mbh", "rIa", "rsCusp", "rsCore"}]
Export["eta.csv", Transpose[{mbhs/MS, \[Eta]s}], "TableHeadings" -> {"Mbh", "eta"}]
Export["etaImp.csv", Transpose[{timeImps/year, \[Eta]Imps}], "TableHeadings" -> {"time", "eta"}]
Export["vweffStarImps.csv", Transpose[{timeImps/year, vweffStarImps}], "TableHeadings" -> {"time", "vw"}]
