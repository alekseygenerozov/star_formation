(* ::Package:: *)

#!/usr/local/bin/MathematicaScript/ -script


AppendTo[$Path,"/Users/aleksey/code/mathematica/star_formation_code"];
Needs["vwEff`"]


(*Impulsive limit mbh=10^6 msun...*)
mbh=10.^6*MS;
timeImps=10.^Range[6.6,10.,0.1]*year;
\[Eta]Imps=\[Eta]Imp/@timeImps;
tmp=vweffTotImp[mbh, #,0.8]&/@timeImps;
(*heating from different sources in the impulsive limit*)
vwcImps=tmp[[;;,1]];
vweffIaImps=tmp[[;;, 2]];
vweffTotImps=tmp[[;;, 3]];
vweffIaImps2=tmp[[;;, 4]];
rslist = Table[rs[mbh, vweffTotImps[[i]], 0.8], {i, 1, Length[timeImps]}];

vweffStarImps=vweffStarImp/@timeImps;
vweffMSPImps=vweffMSPImp/@timeImps;
vweffNoBHImps=Sqrt[vweffTotImps^2-vwcImps^2];

radiiIa=Table[radiusIaImp[mbh, timeImps[[i]]], {i, 1, Length[timeImps]}];
(*Heating over cooling ratio*)
hcImps=Table[hc[mbh,vweffNoBHImps[[i]], 0.8, \[Eta]Imps[[i]]], {i, 1, Length[timeImps]}];
hcIaImps=Table[hcIa[mbh,vweffIaImps2[[i]], radiiIa[[i]], 0.8, \[Eta]Imps[[i]]], {i, 1, Length[timeImps]}];
hcOverall=Table[If[(vweffIaImps2[[i]]>vweffStarImps[[i]])&&(vweffIaImps2[[i]]>vwcImps[[i]])&& (radiiIa[[i]]>rslist[[i]]), hcIaImps[[i]], hcImps[[i]]], {i, 1, Length[timeImps]}];

Transpose[{timeImps/year, \[Eta]Imps}] //Export[ "/Users/aleksey/Second_Year_Project/star_formation/etaImp.csv", #,\
TableHeadings->{"time", "eta"}]&;
Transpose[{timeImps/year*10.^5, vweffIaImps, vweffIaImps2, vweffStarImps, vweffMSPImps, vwcImps, vweffNoBHImps}/10.^5]//Export["/Users/aleksey/Second_Year_Project/star_formation/vwSourcesImps6.csv", #, TableHeadings->{"Time","Ias", "IasCorrected", "Stars", "MSPs", "Compton", "Total"}]&;
Transpose[{timeImps/year, hcImps, hcIaImps, hcOverall}]//Export["/Users/aleksey/Second_Year_Project/star_formation/hcImps6.csv", #, TableHeadings->{"Time", "Cusp", "Ia", "Overall"}]&;


(*Impulsive limit mbh=10^8 msun...*)
mbh=10.^8 MS;
gamma=0.8
timeImps=10.^Range[6.6,10.,0.1]*year;
\[Eta]Imps=\[Eta]Imp/@timeImps;
tmp=vweffTotImp[mbh, #, gamma]&/@timeImps;
(*heating from different sources in the impulsive limit*)
vwcImps=tmp[[;;,1]];
vweffIaImps=tmp[[;;, 2]];
vweffTotImps=tmp[[;;, 3]];
vweffIaImps2=tmp[[;;, 4]];
rslist = Table[rs[mbh, vweffTotImps[[i]], gamma], {i, 1, Length[timeImps]}];

vweffStarImps=vweffStarImp/@timeImps;
vweffMSPImps=vweffMSPImp/@timeImps;
vweffNoBHImps=Sqrt[vweffTotImps^2-vwcImps^2];

radiiIa=Table[radiusIaImp[mbh, timeImps[[i]]], {i, 1, Length[timeImps]}];
(*Heating over cooling ratio*)
hcImps=Table[hc[mbh,vweffNoBHImps[[i]], gamma, \[Eta]Imps[[i]]], {i, 1, Length[timeImps]}];
hcIaImps=Table[hcIa[mbh,vweffIaImps2[[i]], radiiIa[[i]], gamma, \[Eta]Imps[[i]]], {i, 1, Length[timeImps]}];
hcOverall=Table[If[(vweffIaImps2[[i]]>vweffStarImps[[i]])&&(vweffIaImps2[[i]]>vwcImps[[i]])&& (radiiIa[[i]]>rslist[[i]]), hcIaImps[[i]], hcImps[[i]]], {i, 1, Length[timeImps]}];

 
Transpose[{timeImps/year*10.^5, vweffIaImps, vweffIaImps2, vweffStarImps, vweffMSPImps, vwcImps, vweffNoBHImps}/10.^5]//Export["/Users/aleksey/Second_Year_Project/star_formation/vwSourcesImps8.csv", #, TableHeadings->{"Time","Ias", "IasCorrected", "Stars", "MSPs", "Compton", "Total"}]&;
Transpose[{timeImps/year, hcImps, hcIaImps, hcOverall}]//Export["/Users/aleksey/Second_Year_Project/star_formation/hcImps8.csv", #, TableHeadings->{"Time", "Cusp", "Ia", "Overall"}]&;
