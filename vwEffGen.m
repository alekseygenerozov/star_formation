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
mdots = Table[mdotsol[mbhs[[i]], vweffTots[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}]; 
mdotsCore = Table[mdotsol[mbhs[[i]], vweffTotsCore[[i]], 0.1, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}]; 
mdotEdds = mdotEdd /@ mbhs; 
eddrs = mdots/mdotEdds; 
eddrsCore = mdotsCore/mdotEdds; 


radiiIa = Table[radiusIa[mbhs[[i]], mhalos[[i]]], {i, 1, Length[mbhs]}]; 
mdotIas = Table[mdotIA[mbhs[[i]], radiiIa[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}]; 
eddrIas = mdotIas/mdotEdds; 
mdotInfs = (\[Eta]s*mbhs)/th; 
eddrInfs = mdotInfs/mdotEdds; 
eddrMaxCools = Table[eddrMaxCool[mbhs[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[\[Eta]s]}]; 


stags = Table[rs[mbhs[[i]], vweffTots[[i]]], {i, 1, Length[mbhs]}]; 


hcs = Table[hc[mbhs[[i]], vweffTots[[i]], 0.8, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}]; 
hcsCore = Table[hc[mbhs[[i]], vweffTots[[i]], 0.1, \[Eta]s[[i]]], 
    {i, 1, Length[mbhs]}]; 


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
Export["hc.csv", Transpose[{mbhs/MS, hcs, hcsCore}], "TableHeadings" -> {"Mbh", "Cusp", "Core"}]
Export["Ia.csv", Transpose[{mbhs/MS, radiiIa, stags}], "TableHeadings" -> {"Mbh", "rIa", "rs"}]
Export["eta.csv", Transpose[{mbhs/MS, \[Eta]s}], "TableHeadings" -> {"Mbh", "eta"}]
Export["etaImp.csv", Transpose[{timeImps/year, \[Eta]Imps}], "TableHeadings" -> {"time", "eta"}]
Export["vweffStarImps.csv", Transpose[{timeImps/year, vweffStarImps}], "TableHeadings" -> {"time", "vw"}]
Export["eddrs.csv", Transpose[{mbhs/MS, eddrs, eddrsCore}], "TableHeadings" -> {"Mbh", "Cusp", "Core"}]
Export["eddrIas.csv", Transpose[{mbhs/MS, eddrIas}], "TableHeadings" -> {"Mbh", "Cusp"}]
Export["eddrInfs.csv", Transpose[{mbhs/MS, eddrInfs}], "TableHeadings" -> {"Mbh", "Cusp"}]
Export["eddrMaxCools.csv", Transpose[{mbhs/MS, eddrMaxCools}], "TableHeadings" -> {"Mbh", "Cusp"}]
