(* ::Package:: *)

(*Impulsive limit mbh=10^6 msun...*)
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
Transpose[{timeImps/year, hcImps, hcIaImps}]//Export["/Users/aleksey/Second_Year_Project/star_formation/hcImps6.csv", #, TableHeadings->{"Time", "Cusp", "Ia"}]&;


(*Impulsive limit mbh=10^8 msun...*)
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
Transpose[{timeImps/year, hcImps, hcIaImps}]//Export["/Users/aleksey/Second_Year_Project/star_formation/hcImps8.csv", #, TableHeadings->{"Time", "Cusp", "Ia"}]&;


(*Export["etaImp.csv", Transpose[{timeImps/year, \[Eta]Imps}], "TableHeadings" -> {"time", "eta"}]
Export["vweffStarImps.csv", Transpose[{timeImps/year, vweffStarImps}], "TableHeadings" -> {"time", "vw"}]*)
