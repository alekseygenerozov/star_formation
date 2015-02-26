(* ::Package:: *)

#!/usr/local/bin/MathematicaScript/ -script


Needs["vwEff`"]


mhalos=10.^Range[10.8, 14.,0.2] MS;
mstarTots=mstarTot/@mhalos;
mbhs=MbhMbulge/@mstarTots


vweffStars = vweffStar /@ mhalos; 
vweffMSPs = vweffMSP /@ mhalos; 
vweffIas = vweffIa /@ mhalos;


tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], 0.8], {i, 1, Length[mhalos]}]
vwc = tmp[[1 ;; All,1]]; 
vweffIasCorrected = tmp[[1 ;; All,2]]; 
vweffTots = tmp[[1 ;; All,3]]; 

tmp = Table[vweffTot[mbhs[[i]], mhalos[[i]], 0.1], {i, 1, Length[mhalos]}]
vwcCore = tmp[[1 ;; All,1]]; 
vweffIasCoreCorrected = tmp[[1 ;; All,2]]; 
vweffTotsCore = tmp[[1 ;; All,3]];


\[Eta]s = \[Eta] /@ mhalos; 


radiiIa = Table[radiusIa[mbhs[[i]], mhalos[[i]]], {i, 1, Length[mbhs]}]; 


stags = Table[rs[mbhs[[i]], vweffTots[[i]]], {i, 1, Length[mbhs]}];


timeImps = 10.^Range[6.6, 10., 0.2]*year; 
vweffStarImps = vweffStarImp /@ timeImps; 
\[Eta]Imps = \[Eta]Imp /@ timeImps; 
mbhsSparse = 10.^Range[6., 9., 1.]*MS; 
eddrImps = Table[mdotsol[mbhsSparse[[j]], vweffStarImps[[i]], 0.8, \[Eta]Imps[[i]]]/
     mdotEdd[mbhsSparse[[j]]], {i, 1, Length[timeImps]}, 
    {j, 1, Length[mbhsSparse]}]; 


SetDirectory["/Users/aleksey/Second_Year_Project/star_formation"]
Export["vwSources.csv", Transpose[{mbhs*(10.^5/MS), vweffIas, vweffStars, vweffMSPs, vwc, vweffTots}/10.^5], 
  "TableHeadings" -> {"Mbh", "Ias", "Stars", "MSPs", "Compton", "Total"}]
Export["vwSourcesCore.csv", Transpose[{mbhs*(10.^5/MS), vweffIas, vweffStars, vweffMSPs, vwcCore, vweffTotsCore}/10.^5], 
  "TableHeadings" -> {"Mbh", "Ias", "Stars", "MSPs", "Compton", "Total"}]
Export["Ia.csv", Transpose[{mbhs/MS, radiiIa}], "TableHeadings" -> {"Mbh", "rIa"}]

Export["eta.csv", Transpose[{mbhs/MS, \[Eta]s}], "TableHeadings" -> {"Mbh", "eta"}]
Export["etaImp.csv", Transpose[{timeImps/year, \[Eta]Imps}], "TableHeadings" -> {"time", "eta"}]
Export["vweffStarImps.csv", Transpose[{timeImps/year, vweffStarImps}], "TableHeadings" -> {"time", "vw"}]
