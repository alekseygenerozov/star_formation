(* ::Package:: *)

BeginPackage["vwEff`"];


MS=1.989 10^33;
th=4.35 10.^17;
H0=1./(th);
\[CapitalOmega]m=0.3;
\[CapitalOmega]\[Lambda]=0.7;
year=3.15 10^7;
ad=5./3.;
G=6.674*10^-8;
Rsun=6.955 10^10;
Lsun=3.846 10^33;
c=3. 10^10;
pc=3.08568*10^18;
mp=1.672621777 10^\[Minus]24;
kb=1.3806488 10^\[Minus]16;
\[Mu]=0.62;
\[Sigma]sb=5.67 10^-5;
me=9.11 10^-31;
\[Sigma]Thomson=0.4 me;



(*Thermal propeties of gas*)
lambdaC[T_]:=1.1 10.^-22(T/10.^6)^-0.7
cooling[\[Rho]_, T_]:=(\[Rho]/(\[Mu] mp))^2 lambdaC[T]

\[Kappa][T_]:=2. 10^-6 T^(5./2.)
cs[T_]:=Sqrt[ kb (ad T)/(\[Mu] mp)]



(*Cosmology*)
zu=10.;
age=th NIntegrate[1/((1+zp)Sqrt[\[CapitalOmega]m (1+zp)^3+\[CapitalOmega]\[Lambda]]), {zp, 0, \[Infinity]}];
(*Fitting formula for look-back time from Peebles; age is returned in years--eq. 13.20*)
tL[z_]:=(age-th 2./(3.(1-\[CapitalOmega]m)^(1./2.)) ArcSinh[((1-\[CapitalOmega]m)/\[CapitalOmega]m)^(1./2.) (1./(1.+z))^(3./2.)])
(*Inverse of lookback time formula above--accepts time in years.*)
z[t_]:=((\[CapitalOmega]m/(1-\[CapitalOmega]m))^(1/2) Sinh[((age-t )/th)((3(1-\[CapitalOmega]m)^(1/2))/2)])^(-2./3.)-1


(*Black hole scaling relations*)
Mhalo[mbh_]:=10.^13 (mbh/(10.^8.18 MS)) ^(1./1.55) MS
Mbh[mhalo_]:=(mhalo/(10.^13 MS))^1.55 10.^8.18 MS
vcirc[mhalo_]:=2. 10^7 10.0^((Log10[mhalo/(10.^12 MS)]-0.15)/3.32)
(*m-sigma relationship--from Gultekin...*)
\[Sigma][mbh_]:=(1./2.)^(1./5.1) (mbh/(10^8 MS))^(1./5.1) 2 10^7
MbhMbulge[mbulge_]:=10.^8.46 (mbulge/(10.^11 MS))^1.05 MS
(*Influence radius*)
rinf[mbh_]:=14*(mbh/(10.^8 MS))^0.6 pc
(*Scaling relation for the break radius for cores*)
rbCore[mbh_]:=106*(mbh/(10.^8 MS))^0.39*pc
(*ratio of break radius to influence radius for core galaxies*)
rbCorerinf[mbh_]:=rbCore[mbh]/rinf[mbh]


(*IMF*)
m0=0.1MS;
\[Alpha]=2.35;
\[Mu]sal[M_]:=(1./16.58) (M/MS)^-\[Alpha] 1/MS;
(*Mean zams mass for the above IMF*)
mavg=0.35 MS;



(*Fits from Ciotti & Ostriker 2007*)
(*Minimum timescale to evolve off of the main sequence*)
tmin=3. 10^6;
(*Fit for the turnoff mass from Ciotti & Ostriker 2007*)
(*Mt0FitCO[t_]:=MS 10.`^(0.0558` Log10[t/year]^2-1.338` Log10[t/year]+7.764`)*)

(*Fit for turnoff mass, asymptotes to result in Ciotti & Ostriker 2007 after a few hundred Myrs.*)
Mt0FitGen=0.24+4.77 10^-6  E^(-4.58 x)-0.34 x+0.068 x^2;
Mt0Fit1[t_]:=MS 10.^(Mt0FitGen/.x-> Log10[t/(10.^9  year)])
(*Inverse of above function*)
tt0Fit1[Mt0_]:=t/.FindRoot[Mt0Fit1[t]==Mt0,{t, 10.^7 year}];
(*Have a floor on the stallar lifetime--this should be around 3 Myr. In the generalized fitting formula for Mt0 above this is approximately (though not exactly) the lifetime of a 100 solar mass star. I find the exact value (tmin) for the 100 MS lifetime below, 
assuming our fiting formula*)
tmin=t/.FindRoot[Mt0Fit1[t]/MS==100., {t, 3. 10^6 year}];
(*If time is lass than tmin above then return 100 MS for the turnoff mass, since the derivative of the turnoff mass is 0, this ensures that the mass return rate from evolved stars for the first 3 Myr is 0.*)

Mt0Fit[t_]:=Piecewise[{{100MS,t<tmin}}, Mt0Fit1[t]];
(*tmin=3. 10^6 year;
Mt0Fit[t_]:=Piecewise[{{100. MS, t<tmin}}, Mt0FitCO[t]];*)
(*Lifetime of the most massive stars*)

dMt0dt=Mt0Fit';
(*Note that the above formula for the turn-off mass does not account for the maximum mass within our IMF*) 
(*dMt0dt[t_]:=If[t<tmin, 0., dMt0dt1[t]]*)
(*Mass loss as a function of time of time for a coevolving population from Ciotti & Ostriker 2007.*)
\[CapitalDelta]M1[t_]:=(0.945 Mt0Fit[t]-0.503 MS );
ttrans=tt0Fit1[9. MS];
\[CapitalDelta]M[t_]:=Piecewise[{{Mt0Fit[t]-1.4 MS,t<=ttrans}},\[CapitalDelta]M1[t] ];


(*Fits from Moster et al. 2012*)
(*Fitting function for the star formation history from Moster et al. 2012*)
f10=2.658 MS/year;
f20=5.507;
f30=-0.915;
f11=11.937 ;
f21=2.437;
f31=0.696;
f12=0.424;
f32=-0.159;
f1[mhalo_]:= f10 Exp[-(Log10[mhalo/MS]-f11)^2/(2 f12^2)]
f2[mhalo_]:=f20 +f21 Log10[(mhalo/MS)/10.^12]
f3[mhalo_]:=10.^(f30+f31 ((mhalo/MS)/10.^12)^f32)
(*Note that we divide by the mean stellar mass (in solar masses) to obtain the rate of star formation--in g s^-1*)
dSdtForm[z_, mhalo_]:=f1[mhalo] (1+z)^-f2[mhalo] Exp[f2[mhalo]/f3[mhalo] z/(1+z)] 
dNdtForm[z_, mhalo_]:=1/mavg dSdtForm[z,mhalo]



(*Stellar mass accreted*)
g10=1.633 MS/year;
g20=0.855;
g11=-2.017;
g21=0.098;
g12=-0.806;
g22=-0.797;

mhaloAcc=10.^13.*MS
g1[mhalo_]:=g10*((mhalo/(10.^12.5  MS))^g11+(mhalo/(10.^12.5 MS))^g12)^-1
g2[mhalo_]:=g20+g21 (mhalo/(10.^12 MS))^g22
dSdtAcc[z_, mhalo_]:=If[mhalo>mhaloAcc, g1[mhalo] Exp[-z/g2[mhalo]],0.]
dNdtAcc[z_, mhalo_]:=dSdtAcc[z, mhalo]/mavg


dSdt[z_, mhalo_]:=dSdtAcc[z, mhalo]+dSdtForm[z, mhalo]
dNdt[z_, mhalo_]:=dNdtAcc[z, mhalo]+dNdtForm[z, mhalo]
(*Total stellar mass accumulated within galaxy including stars which are accreted in mergers
To find the total stellar mass at present one would have to account for stars dying off. 
However, to compute this requires a more computationally intensive double integral 
and to a good approximation the total stellar mass at present is the total stellar mass
 formed.*)
mstarTotAcc[M_]:=mavg If[M>mhaloAcc, NIntegrate[Abs[dNdtAcc[z[t], M]], {t, 0., tL[zu]}, Method->"DoubleExponential"], 0.]
mstarTotForm[M_]:=mavg NIntegrate[Abs[dNdtForm[z[t], M]], {t, 0., tL[zu]}, Method->"DoubleExponential"]
mstarTot[M_]:=mstarTotAcc[M]+mstarTotForm[M]
(*mstarTot[M_]:=mavg NIntegrate[Abs[dNdt[z[t], M]], {t, 0., tL[zu]}, Method->"DoubleExponential"]*)


(*Stellar properties*)
(*Luminosity and radius in units of solar from KW scaling relations*)
Lstar[Mstar_]:=If[Mstar>MS, (Mstar/MS)^3.2,(Mstar/MS)^2.5]Lsun
Rstar[Mstar_]:=If[Mstar>MS,  (Mstar/MS)^0.57, (Mstar/MS)^0.8] Rsun
Teff[Mstar_]:=(Lstar[Mstar]/(4. \[Pi] \[Sigma]sb Rstar[Mstar]^2))^(1/4)
gstar[Mstar_]:=G Mstar/Rstar[Mstar]^2
gsun=G MS/Rsun^2;
vesc[Mstar_]:=Sqrt[2.0  G Mstar/(Rstar[Mstar])]
(*Wind velocities for main sequence winds (km/s). Currently just use value for the Sun. On the lower main sequence (<1.2 Msun), the wind velocity is a weak function of mass and rotation rate. (Although SMR and FMR stars) have different significant differences. See Holzwarth and Jardine 2008*)
vwMS[Mstar_]:=4.3 10^7 Sqrt[(Mstar/MS)/(Rstar[Mstar]/Rsun)]
(*Number of stars in stellar population from Voss--100 stars with stellar masses between 8 and 120 Msun with a Salpeter IMF. Should resolve slight discrepancy between upper mass limit in Voss and the one used here.*)
nVoss=100./NIntegrate[\[Mu]sal[Mstar],{Mstar, 8. MS, 100. MS}];
(*Wind power Voss et al. 2009--energy budget is actually dominated by Wolf-Rayet stars.*)
twrcut=10.^8*year;
edotWR[t_]:=Piecewise[{{ 10.^36.1/(nVoss/100.),t<4.0 10^6 year}, {0.,t>twrcut}}, 10.^36.1/(nVoss/100.)*(t/(4. 10^6 year))^-3.73]

(*Reimers' prescription for stellar mass loss return rate. Notice that this would significantly overestimate the mass return rate for the sun, although from observation there is a very wide scatter about this relation. Even an individual star could show orders of magnitude variability in mdot.*)
mdotStarReimers[Mstar_]:=4 10^-13 (Mstar/MS)^-1 (Lstar[Mstar]/Lsun)(Rstar[Mstar]/Rsun) MS/year
(*Improved prescription from Schroder & Cuntz 2005*)
mdotStar[Mstar_]:=mdotStarReimers[Mstar]
mdotStar[Mstar_]:=8 10^-14 (Mstar/MS)^-1 (Lstar[Mstar]/Lsun)(Rstar[Mstar]/Rsun)(Teff[Mstar]/4000.)^3.5 (1+gsun/(4300. gstar[MS])) MS/year
enStar[Mstar_]:=0.5*mdotStar[Mstar]*\[Mu]sal[Mstar]*vwMS[Mstar]^2
(*Fitting formula to mass-loss rate per star from Voss 2009*)
mdotVoss[t_]:=If[t> 4. 10^6 year, 1./(nVoss/100.) 10.^-6.1 (t/(4. 10^6 year))^-1.8 MS/year, 1./(nVoss/100.) 10.^-6.1 MS/year]



(*Supernova properties*)

(*Ia properties*)
DTD[t_]:=0.03((t/(10.^8 year))^(-1.12 ))*(1./(10.^10 MS))*(1/year)
(*Ia rate per stellar mass--from star formed within galaxty*)
rateIaSpecific[t_?NumericQ, mhalo_]:=NIntegrate[DTD[t1]*dSdtForm[z[t1],mhalo], {t1,  Max[t,tmin], tL[zu]}]\
/NIntegrate[dSdtForm[z[t1],mhalo], {t1,  Max[t,tmin], tL[zu]}]
(*Total Ia rate from in situ star formation and then accreted stars*)
rateIaFormTot[mhalo_]:=rateIaSpecific[0., mhalo]*mstarTotForm[mhalo]
rateIaAccTot[mhalo_]:=NIntegrate[rateIaSpecific[t1, mhalo]*dSdtAcc[z[t1], mhalo], {t1, 0., tL[3.]}]
rateIa[mhalo_]:=(rateIaFormTot[mhalo]+rateIaAccTot[mhalo])/mstarTot[mhalo]

(*Ia rate within the impulsive limit*)
rateIaImp[t_]:=DTD[t]
fracII=NIntegrate[\[Mu]sal[mstar], {mstar, 8.MS, 100.MS}];
(*Type II properties--so far just include the stars which were formed in situ.*)
RII[mhalo_]:=dNdtForm[0, mhalo]fracII
RIIsp[mhalo_]:=RII[mhalo]/mstarTot[mhalo]
(*Radius for which time between successive Ias is the same as the dynamical time*)
radiusIa[mbh_, mhalo_]:=
Sqrt[G /(\[Sigma][mbh] rateIa[mhalo])]
radiusIaImp[mbh_, t_]:=
Sqrt[G/(\[Sigma][mbh]DTD[t])]


(*Consistenncy of convention Mbh vs. M...*)
radiusII[mbh_, rateII_, \[CapitalGamma]_:1]:=(rinf[mbh]^(2-\[CapitalGamma])/(rateII*mbh) Sqrt[G mbh])^(1./(3.5-\[CapitalGamma]))
radiusII[mbh_, \[CapitalGamma]_:1]:=radiusII[mbh, RIIsp[Mhalo[mbh]], \[CapitalGamma]]


(*Mass and energy injectiuon as a function of Halo mass*)
mdotImp[t_]:=Abs[Mt0Fit'[t]] \[CapitalDelta]M[t] \[Mu]sal[Mt0Fit[t]]
(*Mass shed by turn-off stars*)
(*Mass injection per star for Moster star formation history truncated at lookback time t.*) 
(*mdotSpecific4[t_?NumericQ, mhalo_]:=NIntegrate[dNdtForm[z[t1],mhalo] Abs[dMt0dt[t1]] \[CapitalDelta]M[t1] \[Mu]sal[Mt0Fit[t1]],{t1,t,tmin,ttrans,tL[zu]}]\
/NIntegrate[dNdtForm[z[tl], mhalo], {tl, t, tL[zu]}]
(*mdotSpecific[t_?NumericQ, mhalo_]:=NIntegrate[dNdtForm[z[t1],mhalo] Abs[dMt0dt[t1]] \[CapitalDelta]M[t1] \[Mu]sal[Mt0Fit[t1]],{t1,t,tmin,ttrans,tL[zu]}]\
/NIntegrate[dNdtForm[z[tl], mhalo], {tl, t, tL[zu]}]*)
mdotSpecific2[t_?NumericQ, mhalo_]:=NIntegrate[dNdtForm[z[t1],mhalo] Abs[dMt0dt[t1]] \[CapitalDelta]M[t1] \[Mu]sal[Mt0Fit[t1]],{t1,t,tL[zu]}]\
/NIntegrate[dNdtForm[z[tl], mhalo], {tl, t, tL[zu]}]
mdotSpecific3[t_?NumericQ, mhalo_]:=(NIntegrate[dNdtForm[z[t1],mhalo] Abs[dMt0dt[t1]] \[CapitalDelta]M[t1] \[Mu]sal[Mt0Fit[t1]],{t1,t,ttrans}(*, Method->"AdaptiveMonteCarlo"*)]+\
NIntegrate[dNdtForm[z[t1],mhalo] Abs[dMt0dt[t1]] \[CapitalDelta]M[t1] \[Mu]sal[Mt0Fit[t1]],{t1,ttrans, tL[zu]}(*, Method->"AdaptiveMonteCarlo"*)])\
/NIntegrate[dNdtForm[z[tl], mhalo], {tl, t, tL[zu]}]
*)
mdotSpecific[t_?NumericQ, mhalo_]:=Piecewise[{{NIntegrate[dNdtForm[z[t1],mhalo] Abs[dMt0dt[t1]] \[CapitalDelta]M[t1] \[Mu]sal[Mt0Fit[t1]],{t1,t,tmin,ttrans,tL[zu]}]\
, t<=tmin}, {NIntegrate[dNdtForm[z[t1],mhalo] Abs[dMt0dt[t1]] \[CapitalDelta]M[t1] \[Mu]sal[Mt0Fit[t1]],{t1,t,ttrans,tL[zu]}], t>tmin&&t<=ttrans}},\
NIntegrate[dNdtForm[z[t1],mhalo] Abs[dMt0dt[t1]] \[CapitalDelta]M[t1] \[Mu]sal[Mt0Fit[t1]],{t1,t,ttrans,tL[zu]}]]/NIntegrate[dNdtForm[z[tl], mhalo], {tl, t, tL[zu]}]


mdotForm[mhalo_]:=mdotSpecific[0, mhalo]*NIntegrate[dNdtForm[z[t1], mhalo], {t1, 0, tL[zu]}]
mdotAcc[mhalo_]:=If[mhalo>mhaloAcc, NIntegrate[mdotSpecific[t1, mhalo]*dNdtAcc[z[t1], mhalo], {t1, 0., tL[3.]}],0]
(*mdotAcc2[mhalo_]:=If[mhalo>mhaloAcc, NIntegrate[mdotSpecific2[t1, mhalo]*dNdtAcc[z[t1], mhalo], {t1, 0., tL[3.]}],0]*)
mdot[mhalo_]:=mdotAcc[mhalo]+mdotForm[mhalo]

(*lookup table for computing heating from main sequence stellar winds*)
lookupTableOrds=10.^Range[-0.8,2, 0.2]MS;
lookupTable=NIntegrate[enStar[ms], {ms, 0.1*MS, #}]&/@lookupTableOrds;
enStarIntInterp=Transpose[{Log10[lookupTableOrds], Log10[lookupTable]}]//Interpolation;
enStarInt[mt0_]:=10.^enStarIntInterp[Log10[mt0]]


I1=NIntegrate[enStar[Mstar], {Mstar,0.1*MS,MS}];
(*Turnoff and main sequence energy injection per star for Moster star formation histories truncated at time t.*)
edotTOSpecific[t_?NumericQ, mhalo_]:=NIntegrate[dNdtForm[z[t1], mhalo]*edotWR[t1], {t1, t, tL[zu]} ]/NIntegrate[dNdtForm[z[t1], mhalo], {t1, t, tL[zu]} ]
edotMSSpecific[t_?NumericQ, mhalo_]:=(NIntegrate[dNdtForm[z[t1], mhalo]*enStarInt[Mt0Fit[t1]], {t1, t, tmin, tL[zu]},\
Method->{"InterpolationPointsSubdivision","MaxSubregions"->1+Length[First@enStarInt["Coordinates"]]}])\
/NIntegrate[dNdtForm[z[t1], mhalo], {t1, t, tL[zu]}]

(*Turnoff and main sequence contributions to energy injection.*)
edotTOForm[mhalo_]:= edotTOSpecific[0, mhalo]*NIntegrate[dNdtForm[z[t1], mhalo], {t1, 0., tL[zu]} ]
edotMSForm[mhalo_]:= edotMSSpecific[0, mhalo]*NIntegrate[dNdtForm[z[t1], mhalo], {t1, 0., tL[zu]} ]
(*Energy injection per accreted star for stars accreted at look-back time t*)
edotTOAcc[mhalo_]:= If[mhalo>mhaloAcc, NIntegrate[dNdtAcc[z[t], mhalo] edotTOSpecific[t, mhalo], {t, 0., 0.99*twrcut} ],0.]
edotMSAcc[mhalo_]:= If[mhalo>mhaloAcc, NIntegrate[dNdtAcc[z[t], mhalo] edotMSSpecific[t, mhalo], {t, 0., tL[3.]} ],0.]
(*edotMSAcc[mhalo_]:=0.5 NIntegrate[dNdtForm[z[t], mhalo]mdotStar[Mstar]\[Mu]sal[Mstar] vwMS[Mstar]^2, {t,0., tL[zu]},{Mstar, 0.1 MS ,MS}]\
+0.5 NIntegrate[dNdtForm[z[t], mhalo]mdotStar[Mstar]\[Mu]sal[Mstar] vwMS[Mstar]^2, {t,0., tL[zu]},{Mstar ,MS, Mt0Fit[t]}];*)

(*Contribution of main sequence stars to energy injected*)
(*Overall effective wind velocity.*)
vweffStar[mhalo_]:=Sqrt[2 (edotMSForm[mhalo]+edotTOForm[mhalo]+edotMSAcc[mhalo]+edotTOAcc[mhalo])/(mdotAcc[mhalo]+mdotForm[mhalo])]

\[Eta][mhalo_]:=(mdotForm[mhalo]+mdotAcc[mhalo])/mstarTot[mhalo] th
(*Energy injected by MS stars in impulsive limit--note that unlike the continuous star formation limit here we have the mass and energy injected per star. Maybe make the impulsive and the continuous limits more consistent.*)
edotMSImp[t_?NumericQ]:=0.5 NIntegrate[enStar[Mstar],{Mstar, 0.1 , Mt0Fit[t]}]
edotTOImp[t_]:=edotWR[t]

mdotMSImp[t_?NumericQ]:= NIntegrate[mdotStar[Mstar]\[Mu]sal[Mstar],{Mstar, 0.1 , Mt0Fit[t]}]

vweffStarImp[t_]:=Sqrt[2 (edotMSImp[t]+edotTOImp[t])/mdotImp[t]]
vweffStarImp2[t_]:=Sqrt[2 (edotMSImp[t]+edotTOImp[t])/(mdotImp[t]+mdotMSImp[t])]
\[Eta]Imp[t_]:=mdotImp[t]/mavg th

(*Effective vws for different heating sources for different star formation histories*)
vweffIa[mhalo_, \[Epsilon]Ia_:0.4]:=Sqrt[(2.th rateIa[mhalo] \[Epsilon]Ia 10.^51)/\[Eta][mhalo]]

(*Effective wind velocity for Ias in the impulsive star formation limit*)
vweffIaImp[t_, \[Epsilon]Ia_:0.4]:=Sqrt[(2. th rateIaImp[t] \[Epsilon]Ia 10.^51)/\[Eta]Imp[t]]

vweffMSP[mhalo_, \[Epsilon]msp_:0.1,Lsd_:10.^34]:=3. 10^6 (\[Epsilon]msp/0.1)^0.5 (Lsd/10.^34)^(1/2) (\[Eta][mhalo]/0.02)^(-1/2)
vweffMSPImp[t_, \[Epsilon]msp_:0.1,Lsd_:10.^34]:=3. 10^6 (\[Epsilon]msp/0.1)^0.5 (Lsd/10.^34)^(1/2) (\[Eta]Imp[t]/0.02)^(-1/2)

edotCompton[mbh_, vw_, \[CapitalGamma]_:1, \[Eta]_:1, Tc_:10.^9]:=4.1*10^-35*nRs[mbh, vw, \[CapitalGamma], \[Eta]]^2 (10.*(mdotsol[mbh, vw, \[CapitalGamma], \[Eta]]^2/mdotEdd[mbh]) c^2)/(nRs[mbh, vw, \[CapitalGamma], \[Eta]] rs[mbh, vw, \[CapitalGamma]]^2)*(Tc)
vwComptonGen[mbh_,vw_, \[CapitalGamma]_:1., \[Eta]_:1, Tc_:10.^9]:=Sqrt[ ((2. th edotCompton[mbh, vw, \[CapitalGamma], \[Eta], Tc])/(\[Eta] rhoStarRs[mbh, vw, \[CapitalGamma]]))]

vwComptonDom[mbh_, \[CapitalGamma]_:1., \[Eta]_:1., Tc_:10.^9]:=vw1/.FindRoot[vwComptonGen[mbh, vw1, \[CapitalGamma], \[Eta], Tc]==vw1, {vw1, 5.*10^7}]


vweffTot[mbh_, mhalo_, \[CapitalGamma]_:1, \[Epsilon]msp_:0.1, Lsd_:10.^34, \[Epsilon]Ia_:0.4, Tc_:10.^9]:=Module[{vw0,  rs1, \[Eta]1, vwIa0, vwc, vwIa},
vw0=(vweffStar[mhalo]^2+vweffMSP[mhalo,\[Epsilon]msp, Lsd]^2)^(1/2);
vwIa0=vweffIa[mhalo, \[Epsilon]Ia];

{vwc, vwIa, Sqrt[vw0^2+vwIa^2+vwc^2]}/.FindRoot[{vwc==vwComptonGen[mbh, Sqrt[vw0^2+vwc^2+vwIa^2], \[CapitalGamma], \[Eta][mhalo],  Tc], vwIa==vwIa0 E^(-radiusIa[mbh, mhalo]/rs[mbh,Sqrt[vw0^2+vwc^2+vwIa^2],\[CapitalGamma]])}, {vwc,vw0}, {vwIa, vwIa0}]

 ]
vweffTotImp[mbh_, t_,\[CapitalGamma]_:1,  \[Epsilon]msp_:0.1,Lsd_:10.^34, \[Epsilon]Ia_:0.4, Tc_:10.^9]:=Module[{vw0,  rs1, \[Eta]1,vwIa0, vwc, vwIa},
vw0=(vweffStarImp[t]^2+vweffMSPImp[t,\[Epsilon]msp, Lsd]^2)^(1/2);
vwIa0=vweffIaImp[t, \[Epsilon]Ia];

{vwc, vwIa, Sqrt[vw0^2+vwIa^2+vwc^2], (vwIa-vwIa0 E^(-radiusIaImp[mbh, t]/rs[mbh,Sqrt[vw0^2+vwc^2+vwIa^2],\[CapitalGamma]]))/vwIa}/.FindRoot[{vwc==vwComptonGen[mbh, Sqrt[vw0^2+vwc^2+vwIa^2], \[CapitalGamma], \[Eta]Imp[t], Tc], vwIa==vwIa0 E^(-radiusIaImp[mbh, t]/rs[mbh,Sqrt[vw0^2+vwc^2+vwIa^2]])}, {vwc,vw0}, {vwIa, vwIa0}]

 ]


(*estimate for the gas density slope at rs*)
densSlope[\[CapitalGamma]_]:=-(1./6.*(1.-4.*(1+\[CapitalGamma])))
\[CapitalGamma]fitM[mbh_]:=0.3*(mbh/10.^8/MS)^-0.24
(*normalized heating rate and critical value of this parameter beyond which rs approaches rb*)
zeta[mbh_,vw_]:=(vw^2.+3*sigma[mbh]^2.)^0.5/(3.**0.5*sigma[mbh])
zetac[rbrinf_, \[CapitalGamma]_]:=(rbrinf)^(0.5*(1.-\[CapitalGamma]))
(**)
(*result for stagnation radius--follows generozov's law unless we are below critical
heating rate.*)
rs[mbh_,vw_,\[CapitalGamma]_:1, rbrinf_:1]:=Module[{zc, z},
	zc=zetc[rbrinf, \[CapitalGamma]];
	z=zeta[mbh, vw];

	If[z>=zc, ((13.+8.\[CapitalGamma])/(4.+2.\[CapitalGamma])-densSlope[\[CapitalGamma]]*(3./(2.+\[CapitalGamma])))G mbh/vw^2/densSlope[\[CapitalGamma]], 
	(rbrinf)*rinf[mbh]
	]
]
(*temperature at the stagnation radius accounting for the velocity dispersion of the black 
hole.*)
tempRs[vw_, \[CapitalGamma]_:1]:=(ad-1)/ad*\[Mu]*mp*((3.+8.*\[CapitalGamma])/(3.+8.*\[CapitalGamma]-6.*densSlope[\[CapitalGamma]]))*vw^2/(2.*kb)

rhoStarRs[mbh_, vw_, \[CapitalGamma]_:1.]:=mbh/((4.*\[Pi]) rinf[mbh]^3)*(2.-\[CapitalGamma])*(rs[mbh,vw, \[CapitalGamma]]/rinf[mbh])^(-1.-\[CapitalGamma])
mencRs[mbh_,vw_, \[CapitalGamma]_:1.]:=mbh*(rs[mbh,vw, \[CapitalGamma]]/rinf[mbh])^(2.-\[CapitalGamma])

(*accretion rate onto BH for a given solution, assuming eta=1*)
LEdd[mbh_]:=(4 \[Pi] G mbh me c)/(\[Sigma]Thomson);
mdotEdd[mbh_]:=(4 \[Pi] G mbh me )/(\[Sigma]Thomson 0.1 c);
mdotsol[mbh_, vw_, \[CapitalGamma]_:1., \[Eta]_:1.]:=\[Eta] mencRs[mbh,vw,\[CapitalGamma]]/th
mdotIA[mbh_, rIa_, \[CapitalGamma]_:1., \[Eta]_:1.]:=\[Eta] mbh/th (rIa/rinf[mbh])^(2.-\[CapitalGamma])

(*mdotsolHalo[mbh_, mhalo_,\[CapitalGamma]_:1]:=mdotsol[mbh, vweffTot[mhalo], \[CapitalGamma], \[Eta][mhalo]]*)

qRs[mbh_,vw_,  \[CapitalGamma]_:1, \[Eta]_:1]:=\[Eta] rhoStarRs[mbh,vw,\[CapitalGamma]]/th
tff[r_, mbh_]:=r^1.5/(G*mbh)^0.5

rhoRs[mbh_,vw_, \[CapitalGamma]_:1, \[Eta]_:1]:=mdotsol[mbh,vw, \[CapitalGamma], \[Eta]]*tff[rs[mbh,vw,\[CapitalGamma]], mbh]/(4.\[Pi]/3.*rs[mbh,vw,\[CapitalGamma]]^3.)
nRs[mbh_,vw_, \[CapitalGamma]_:1, \[Eta]_:1]:=rhoRs[mbh,vw, \[CapitalGamma], \[Eta]]/(\[Mu]*mp)


heatingRs[mbh_,vw_, \[CapitalGamma]_:1, \[Eta]_:1]:=0.5 qRs[mbh,vw, \[CapitalGamma],\[Eta]]*((3.5)/(3.5-densSlope[\[CapitalGamma]]))*vw^2 
coolingRs[mbh_, vw_, \[CapitalGamma]_:1., \[Eta]_:1.]:=(rhoRs[mbh, vw, \[CapitalGamma], \[Eta]]/(\[Mu]*mp))^2*lambdaC[tempRs[vw, \[CapitalGamma]]]
hc[mbh_,vw_, \[CapitalGamma]_:1., \[Eta]_:1.]:=heatingRs[mbh,vw, \[CapitalGamma], \[Eta]]/coolingRs[mbh,vw,\[CapitalGamma], \[Eta]]

(*Maximum Mdot befor thermal instability sets in*)
vwMaxCool[mbh_, \[CapitalGamma]_:1, \[Eta]_:1, hcCrit_:1]:=vw1/.FindRoot[hc[mbh, vw1, \[CapitalGamma], \[Eta]]==hcCrit, {vw1,3.*10^7.}]
mdotMaxCool[mbh_, \[CapitalGamma]_:1, \[Eta]_:1, hcCrit_:1]:=mdotsol[mbh, vwMaxCool[mbh, \[CapitalGamma], \[Eta], hcCrit], \[CapitalGamma], \[Eta]]
mdotCompton[mbh_, \[CapitalGamma]_:1., \[Eta]_:1., Tc_:10.^9]:=mdotsol[mbh, vwComptonDom[mbh, \[CapitalGamma], \[Eta], Tc], \[CapitalGamma], \[Eta]]


EndPackage[]
