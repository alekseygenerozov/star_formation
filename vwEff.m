(* ::Package:: *)

BeginPackage["vwEff`"]


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
zu=10;
age=th NIntegrate[1/((1+zp)Sqrt[\[CapitalOmega]m (1+zp)^3+\[CapitalOmega]\[Lambda]]), {zp, 0, \[Infinity]}];
(*Fitting formula for look-back time from Peebles; age is returned in years--eq. 13.20*)
tL[z_]:=(age-th 2./(3.(1-\[CapitalOmega]m)^(1./2.)) ArcSinh[((1-\[CapitalOmega]m)/\[CapitalOmega]m)^(1./2.) (1./(1.+z))^(3./2.)])
(*Inverse of lookback time formula above--accepts time in years.*)
z[t_]:=((\[CapitalOmega]m/(1-\[CapitalOmega]m))^(1/2) Sinh[((age-t )/th)((3(1-\[CapitalOmega]m)^(1/2))/2)])^(-2./3.)-1


(*Black hole scaling relations*)
Mhalo[Mbh_]:=10.^13 (Mbh/(10.^8.18 MS)) ^(1./1.55) MS
Mbh[Mhalo_]:=(Mhalo/(10.^13 MS))^1.55 10.^8.18 MS
vcirc[Mhalo_]:=2. 10^7 10.0^((Log10[Mhalo/(10.^12 MS)]-0.15)/3.32)
(*m-sigma relationship--from Gultekin...*)
\[Sigma][Mbh_]:=(1./2.)^(1./5.1) (Mbh/(10^8 MS))^(1./5.1) 2 10^7
(*Influence radius*)
rinf[Mbh_]:=14*(Mbh/(10.^8 MS))^0.6 pc




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
tmin=t/.FindRoot[Mt0Fit1[t]==100.MS, {t, 3. 10^6 year}];
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
f1[M_]:= f10 Exp[-(Log10[M/MS]-f11)^2/(2 f12^2)]
f2[M_]:=f20 +f21 Log10[(M/MS)/10.^12]
f3[M_]:=10.^(f30+f31 ((M/MS)/10.^12)^f32)
(*Note that we divide by the mean stellar mass (in solar masses) to obtain the rate of star formation--in g s^-1*)
dSdt[z_, M_]:=f1[M] (1+z)^-f2[M] Exp[f2[M]/f3[M] z/(1+z)] 
dNdt[z_, M_]:=1/mavg dSdt[z,M]
(*To find the total stellar mass at present one would have to account for stars dying off. However, to compute this requires a more computationally intensive double integral and to a good approximation the total stellar mass at present is the total stellar mass formed.--this is the commented out function mstarTotComp*)
mstarTot[M_]:=mavg NIntegrate[Abs[dNdt[z[t], M]], {t, 0., tL[zu]}, Method->"DoubleExponential"]



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
edotWR[t_]:=Piecewise[{{ 10.^36.1/(nVoss/100.),t<4.0 10^6 year}}, 10.^36.1/(nVoss/100.)*(t/(4. 10^6 year))^-3.73]

(*Reimers' prescription for stellar mass loss return rate. Notice that this would significantly overestimate the mass return rate for the sun, although from observation there is a very wide scatter about this relation. Even an individual star could show orders of magnitude variability in mdot.*)
mdotStarReimers[Mstar_]:=4 10^-13 (Mstar/MS)^-1 (Lstar[Mstar]/Lsun)(Rstar[Mstar]/Rsun) MS/year
(*Improved prescription from Schroder & Cuntz 2005*)
mdotStar[Mstar_]:=mdotStarReimers[Mstar]
mdotStar[Mstar_]:=8 10^-14 (Mstar/MS)^-1 (Lstar[Mstar]/Lsun)(Rstar[Mstar]/Rsun)(Teff[Mstar]/4000.)^3.5 (1+gsun/(4300. gstar[MS])) MS/year
(*Fitting formula to mass-loss rate per star from Voss 2009*)
mdotVoss[t_]:=If[t> 4. 10^6 year, 1./(nVoss/100.) 10.^-6.1 (t/(4. 10^6 year))^-1.8 MS/year, 1./(nVoss/100.) 10.^-6.1 MS/year]



(*Supernova properties*)

(*Ia properties*)
DTD[t_]:=0.03((t/(10.^8 year))^(-1.12 ))*(1./(10.^10 MS))*(1/year)
(*Ia rate per stellar mass--continuous star formation limit*)
rateIa[M_]:=NIntegrate[DTD[t]*dSdt[z[t],M], {t,  tmin, tL[zu]}]/NIntegrate[dSdt[z[t],M], {t,  tmin, tL[zu]}]
(*Ia rate within the impulsive limit*)
rateIaImp[t_]:=DTD[t]
fracII=NIntegrate[\[Mu]sal[mstar], {mstar, 8.MS, 100.MS}];
(*Type II properties*)
RII[M_]:=dNdt[0, M]fracII
RIIsp[M_]:=RII[M]/mstarTot[M]
(*Radius for which time between successive Ias is the same as the dynamical time*)
radiusIa[M_]:=Module[{mbh},
mbh=Mbh[M];
Sqrt[G /(\[Sigma][mbh] rateIa[M])]
]
radiusIaImp[Mbh_, t_]:=
Sqrt[G/(\[Sigma][Mbh]DTD[t])]


(*Consistenncy of convention Mbh vs. M...*)
radiusII[Mbh_, rateII_, \[CapitalGamma]_:1]:=(rinf[Mbh]^(2-\[CapitalGamma])/(rateII*Mbh) Sqrt[G Mbh])^(1./(3.5-\[CapitalGamma]))
radiusII[Mbh_, \[CapitalGamma]_:1]:=radiusII[Mbh, RIIsp[Mhalo[Mbh]], \[CapitalGamma]]





(*Mass and energy injectiuon as a function of Halo mass*)
mdotImp[t_]:=Abs[Mt0Fit'[t]] \[CapitalDelta]M[t] \[Mu]sal[Mt0Fit[t]]
(*Mass shed by turn-off stars*)
mdot[M_]:=NIntegrate[dNdt[z[t], M]Abs[dMt0dt[t]] \[CapitalDelta]M[t] \[Mu]sal[Mt0Fit[t]], {t,  0., tL[zu]}, Method->"DoubleExponential"]
(*Turnoff contribution to energy*)
edotTO[M_]:= NIntegrate[dNdt[z[t], M] edotWR[t], {t, 0., tL[zu]} ]
I1=0.5*NIntegrate[mdotStar[Mstar]\[Mu]sal[Mstar] vwMS[Mstar]^2, {Mstar, 0.1 MS ,100. MS}];
edotMS[M_]:=(*I1*NIntegrate[dNdt[z[t], M], {t,0, tmin}]+*)0.5 NIntegrate[dNdt[z[t], M]mdotStar[Mstar]\[Mu]sal[Mstar] vwMS[Mstar]^2, {t,0., tL[zu]},{Mstar, 0.1 MS ,MS}]+0.5 NIntegrate[dNdt[z[t], M]mdotStar[Mstar]\[Mu]sal[Mstar] vwMS[Mstar]^2, {t,0., tL[zu]},{Mstar ,MS, Mt0Fit[t]}];
(*Contribution of main sequence stars to energy injected*)
(*Overall effective wind velocity.*)

vweffStar[M_]:=Sqrt[2 (edotMS[M]+edotTO[M])/mdot[M]]
vweffMS[M_]:=Sqrt[2 (edotMS[M])/mdot[M]]
\[Eta][M_]:=mdot[M]/mstarTot[M] th
(*Energy injected by MS stars in impulsive limit--note that unlike the continuous star formation limit here we have the mass and energy injected per star. Maybe make the impulsive and the continuous limits more consistent.*)
edotMSImp[t_?NumericQ]:=0.5 NIntegrate[mdotStar[Mstar]\[Mu]sal[Mstar] vwMS[Mstar]^2,{Mstar, 0.1 , Mt0Fit[t]}]
edotTOImp[t_]:=edotWR[t]

mdotMSImp[t_?NumericQ]:= NIntegrate[mdotStar[Mstar]\[Mu]sal[Mstar],{Mstar, 0.1 , Mt0Fit[t]}]

vweffStarImp[t_]:=Sqrt[2 (edotMSImp[t]+edotTOImp[t])/mdotImp[t]]
vweffStarImp2[t_]:=Sqrt[2 (edotMSImp[t]+edotTOImp[t])/(mdotImp[t]+mdotMSImp[t])]
\[Eta]Imp[t_]:=mdotImp[t]/mavg th

(*Effective vws for different heating sources for different star formation histories*)
vweffIa[M_, \[Epsilon]Ia_:0.4]:=Sqrt[(2.th rateIa[M] \[Epsilon]Ia 10.^51)/\[Eta][M]]

(*Effective wind velocity for Ias in the impulsive star formation limit*)
vweffIaImp[t_, \[Epsilon]Ia_:0.4]:=Sqrt[(2. th rateIaImp[t] \[Epsilon]Ia 10.^51)/\[Eta]Imp[t]]

vweffMSP[M_, \[Epsilon]msp_:0.1,Lsd_:10.^34]:=3. 10^6 (\[Epsilon]msp/0.1)^0.5 (Lsd/10.^34)^(1/2) (\[Eta][M]/0.02)^(-1/2)
vweffMSPImp[t_, \[Epsilon]msp_:0.1,Lsd_:10.^34]:=3. 10^6 (\[Epsilon]msp/0.1)^0.5 (Lsd/10.^34)^(1/2) (\[Eta]Imp[t]/0.02)^(-1/2)
(*Total vweff including Compton heating--note that for now I have adopted a fixed efficiency of 10^-2. Note that this Ia heating is not included, even without the Compton heating the stagation radius is generally inside of the Ia radius across our parameter space. Note that we need a more consistent convention for Mbh vs. Mhalo as function arguments.*)
vwCompton[M_,vw_, \[CapitalGamma]_:1., \[Epsilon]_:0.01, Tc_:10.^9]:=Module[{mbh, Tc1},
mbh=Mbh[M];
3.14 10^7*0.43^(0.5*(5.-\[CapitalGamma])) (mbh/10.^8/MS)^(-0.60+0.2*(5. -\[CapitalGamma])) Sqrt[Tc/10.^9] (vw*\[Chi][mbh,vw]/(5. 10^7))^(-1.5+\[CapitalGamma]) Sqrt[\[Epsilon]/10.^-2] Sqrt[\[Eta][M]/0.02]/Sqrt[2.-\[CapitalGamma]]
]
vwComptonImp[Mbh_,t_,vw_, \[CapitalGamma]_:1., \[Epsilon]_:0.01, Tc_:10.^9]:=Module[{ Tc1},
3.14 10^7*0.43^(0.5*(5.-\[CapitalGamma])) (Mbh/10.^8/MS)^(-0.60+0.2*(5. -\[CapitalGamma])) Sqrt[Tc/10.^9] (vw*\[Chi][Mbh,vw]/(5. 10^7))^(-1.5+\[CapitalGamma]) Sqrt[\[Epsilon]/10.^-2] Sqrt[\[Eta]Imp[t]/0.02]/Sqrt[2.-\[CapitalGamma]]
]


vweffTot[M_,\[CapitalGamma]_:1,  \[Epsilon]msp_:0.1,Lsd_:10.^34, \[Epsilon]Ia_:0.4, \[Epsilon]_:0.01 , Tc_:10.^9]:=Module[{vw0, mbh, rs1, \[Eta]1,vwIa0, vwc, vwIa},
vw0=(vweffStar[M]^2+vweffMSP[M,\[Epsilon]msp, Lsd]^2)^(1/2);
vwIa0=vweffIa[M, \[Epsilon]Ia];
mbh=Mbh[M];
{vwc, vwIa, Sqrt[vw0^2+vwIa^2+vwc^2]}/.FindRoot[{vwc==vwCompton[M, Sqrt[vw0^2+vwc^2+vwIa^2], \[CapitalGamma], \[Epsilon], Tc], vwIa==vwIa0 E^(-radiusIa[M]/rs[mbh,Sqrt[vw0^2+vwc^2+vwIa^2]])}, {vwc,vw0}, {vwIa, vwIa0}]

 ]
vweffTotImp[Mbh_, t_,\[CapitalGamma]_:1,  \[Epsilon]msp_:0.1,Lsd_:10.^34, \[Epsilon]Ia_:0.4, \[Epsilon]_:0.01 , Tc_:10.^9]:=Module[{vw0,  rs1, \[Eta]1,vwIa0, vwc, vwIa},
vw0=(vweffStarImp[t]^2+vweffMSPImp[t,\[Epsilon]msp, Lsd]^2)^(1/2);
vwIa0=vweffIaImp[t, \[Epsilon]Ia];

{vwc, vwIa, Sqrt[vw0^2+vwIa^2+vwc^2], (vwIa-vwIa0 E^(-radiusIaImp[Mbh, t]/rs[Mbh,Sqrt[vw0^2+vwc^2+vwIa^2]]))/vwIa}/.FindRoot[{vwc==vwComptonImp[Mbh,t, Sqrt[vw0^2+vwc^2+vwIa^2], \[CapitalGamma], \[Epsilon], Tc], vwIa==vwIa0 E^(-radiusIaImp[Mbh, t]/rs[Mbh,Sqrt[vw0^2+vwc^2+vwIa^2]])}, {vwc,vw0}, {vwIa, vwIa0}]

 ]


(*Generozov's law and gas properties at rs*)
\[Chi][Mbh_,vw_]:=(1+0.12(Mbh/(10.^8 MS))^0.39 (vw/(5. 10^7))^-2)^(1/2)
(*Halo-MBH relation from Bandara et al. 2009*)

rs[Mbh_,vw_]:=3.5 G Mbh/(vw^2 \[Chi][Mbh,vw]^2)
A[\[CapitalGamma]_]:=(4. ad-(ad-1)(1+\[CapitalGamma]))/(2(ad -1))


tempRs[Mbh_,vw_]:=(ad-1)/ad*\[Mu]*mp*vw^2*\[Chi][Mbh,vw]^2./kb

rhoStarRs[Mbh_, vw_, \[CapitalGamma]_:1.]:=Mbh/((4.*\[Pi]) rinf[Mbh]^3)*(2.-\[CapitalGamma])*(rs[Mbh,vw]/rinf[Mbh])^(-1.-\[CapitalGamma])
mencRs[Mbh_,vw_, \[CapitalGamma]_:1.]:=Mbh*(rs[Mbh,vw]/rinf[Mbh])^(2.-\[CapitalGamma])

(*accretion rate onto BH for a given solution, assuming eta=1*)
mdotEdd[Mbh_]=(4 \[Pi] G Mbh me )/(\[Sigma]Thomson 0.1 c);
mdotsol[Mbh_, vw_, \[CapitalGamma]_:1., \[Eta]_:1.]:=\[Eta] mencRs[Mbh,vw,\[CapitalGamma]]/th
mdotIA[Mbh_, rIa_, \[CapitalGamma]_:1., \[Eta]_:1.]:=\[Eta] Mbh/th (rIa/rinf[Mbh])^(2.-\[CapitalGamma])

mdotsolHalo[M_,\[CapitalGamma]_:1]:=mdotsol[Mbh[M], vweffTot[M], \[CapitalGamma], \[Eta][M]]

qRs[Mbh_,vw_,  \[CapitalGamma]_:1, \[Eta]_:1]:=\[Eta] rhoStarRs[Mbh,vw, \[CapitalGamma]]/th

tff[r_, Mbh_]:=r^1.5/(G*Mbh)^0.5

rhoRs[Mbh_,vw_, \[CapitalGamma]_:1, \[Eta]_:1]:=mdotsol[Mbh,vw, \[CapitalGamma]]*tff[rs[Mbh,vw], Mbh]/(4.\[Pi]/3.*rs[Mbh,vw]^3.)


heatingRs[Mbh_,vw_, \[CapitalGamma]_:1, \[Eta]_:1]:=0.5 qRs[Mbh,vw, \[CapitalGamma],\[Eta]] vw^2  \[Chi][Mbh,vw]^2
coolingRs[Mbh_, vw_, \[CapitalGamma]_:1., \[Eta]_:1.]:=(rhoRs[Mbh, vw, \[CapitalGamma], \[Eta]]/(\[Mu]*mp))^2*lambdaC[tempRs[Mbh, vw]]
hc[Mbh_,vw_, \[CapitalGamma]_:1., \[Eta]_:1.]:=heatingRs[Mbh,vw, \[CapitalGamma], \[Eta]]/coolingRs[Mbh,vw,\[CapitalGamma], \[Eta]]

(*Maximum Eddington ratio which may be achieved before we obtain thermal instability assuming a critical heating to cooling ratio of of hcCrit*) 
eddrMaxCool[Mbh_, \[CapitalGamma]_:1, \[Eta]_:1, hcCrit_:1]:=6.58 10^-5 0.43^(3.7-6.29/(3.7 - \[CapitalGamma])) (\[Eta]/0.02)^(1.7/(3.7-\[CapitalGamma])) (Mbh/(10.^8 MS))^((1.36-0.68 \[CapitalGamma])/(3.7-\[CapitalGamma])) (675.938-337.969\[CapitalGamma])^(-((2 (2-\[CapitalGamma]))/(-7.4+2 \[CapitalGamma]))) hcCrit^((2 (2-\[CapitalGamma]))/(-7.4`+ 2.` \[CapitalGamma]))



EndPackage[ ]
