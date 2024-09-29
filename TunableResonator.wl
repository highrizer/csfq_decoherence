(* ::Package:: *)

(* Wolfram Language Package *)

BeginPackage["TunableResonator`"]
(* Exported symbols added here with SymbolName::usage *)  

qubitResonatorCouplingCoeffs::usage
TResVacuumFluctuations::usage
resonatorFrequency::usage
rfsquidLeff::usage
ICircfromIBias::usage

Begin["`Private`"] (* Begin Private Context *) 

phi00 = 3.2897327`*^-16;
hbar = 1.0545606529268985`*^-34;

(*The following four functions are copied from Design/Various.../QubitReadout.m*)
rfsquidClassicalPhaseSols[iB0_, f0_, beta0_] :=
    Module[ {iB = iB0, f = f0, beta = beta0, gammaJ, sols, listroots},
        sols = NSolve[
          iB == Sin[gammaJ] + 1/beta*(gammaJ + 2*Pi*f) && 
           gammaJ >= -2*Pi*f + beta*iB - beta && 
           gammaJ <= -2*Pi*f + beta*iB + beta, {gammaJ}, Reals];
        listroots = (gammaJ /. sols) // Flatten;
        listroots
    ]
  

rfsquidClassicalStablePhaseSols[iB0_, f0_, beta0_] :=
    Module[ {iB = iB0, f = f0, beta = beta0, classicalSolutions,filterStableSolutions,classicalStableSolutions,
        energyExpression,classicalStableSolutionsOrdered,classicalStableSolutionsOrderedEnergies},
        classicalSolutions = rfsquidClassicalPhaseSols[iB, f, beta];
        filterStableSolutions = ((Cos[#] + 1/beta) > 0)&;
        classicalStableSolutions = Select[classicalSolutions,filterStableSolutions];
        energyExpression = (-Cos[#] + 1/2/beta * (# + 2*Pi*f - beta*iB)^2 - (beta*iB)^2/2/beta)&;
        classicalStableSolutionsOrdered = SortBy[classicalStableSolutions, energyExpression];
        classicalStableSolutionsOrderedEnergies = energyExpression/@classicalStableSolutionsOrdered;
        <|"solutions"->classicalStableSolutionsOrdered,"energies"->classicalStableSolutionsOrderedEnergies|>
    ]


rfsquidLeff[IB0_, f0_, beta0_,Ic0_,phi00_] :=
    Module[ {IB = IB0, f = f0, beta = beta0, Ic = Ic0, phi0 = phi00,lowestEnergygammaJ,Leff},
        lowestEnergygammaJ = rfsquidClassicalStablePhaseSols[IB/Ic, f, beta]["solutions"][[1]];
        Leff = phi0/Ic/(Cos[lowestEnergygammaJ] + 1/beta);
        Leff
    ]

resonanceFrequency[c0_,l0_,Z00_,L0_] :=
    Module[ {c = c0,l = l0,Z0 = Z00,L = L0,sols,omegar,root},
        sols = FindRoot[
          Exp[2*I*omegar/c*l]==-(1-I*omegar*L/Z0)/(1+I*omegar*L/Z0), {omegar,0,Pi*c/2/l}];
        root = (omegar /. sols);
        (root//Abs)/2/Pi
    ]



ICircfromIBias[IB0_?NumericQ,f0_,beta0_,Ic0_] :=
    Module[ {IB = IB0, f = f0, beta = beta0, Ic = Ic0, lowestEnergygammaJ,Icirc},
        lowestEnergygammaJ = rfsquidClassicalStablePhaseSols[IB/Ic, f, beta]["solutions"][[1]];
        Icirc = IB-Ic*Sin[lowestEnergygammaJ];
        Icirc
    ]


ICircfromIbiasPolyCoeffs[fr_,beta_,Ic_] :=
    Module[ {interpFunc, zerothOrderCoeff,firstOrderCoeff,secondOrderCoeff},
        interpFunc = Interpolation[Table[{ib,ICircfromIBias[ib,fr,beta,Ic]},{ib,-Ic,Ic,Ic/1000}]];
        zerothOrderCoeff = interpFunc[0];
        firstOrderCoeff = interpFunc'[0];
        secondOrderCoeff = interpFunc''[0]/2;
        {zerothOrderCoeff, firstOrderCoeff,secondOrderCoeff}
    ]

(*current going through SQUID as written as 
c*I*(adagger-a) where c is a function of resonator properties, 
given in the coefficient given in this funcion
resonator is assumed to be lambda/4 terminated by rf-SQUID
rf-SQUID is modelled as an ideal inductor with inductance L in the resonator.
*)
(*see derivation in QEOWaterlooSoftware\TunableResonatorInteraction\TransmissionLineQuantization.nb*)
IResonator[\[Omega]_,c_,Z0_,l_,L_] :=
    (c Z0  Sqrt[hbar \[Omega]] )/(Sqrt[c Z0  ] Sqrt[c L Z0+l (Z0^2+L^2 \[Omega]^2)])
    

(*Assuming resonator interact with Qubit inductively, 
via the circulating current of the rf-SQUID,
 interaction is (c0+c1(adagger-a)+c2(adagger-a)^2 * dH/df, 
 where f is the reduced flux of the system the SQUID is coupled to.
 Units: dimensionless. The coefficients form an expression such as g1*(a+adag)*dHqb/dfz. 
 If fz is dimensionless (frustration), then the Hamiltonian unit is that of Hqb.
*)
qubitResonatorCouplingCoeffs[mrqz_, cr_,Z0r_,lr_,betar_,Icr_,fr_] :=
    Module[{omegar, effectiveL, biasCurrent,zerothOrderCoupling,
        firstOrderCoupling,secondOrderCoupling,circulatingCurrentCoeffs},
        effectiveL = rfsquidLeff[0,fr, betar,Icr,phi00];
        omegar = resonanceFrequency[cr,lr,Z0r,effectiveL]*2*Pi;
        biasCurrent = IResonator[omegar,cr,Z0r,lr,effectiveL];
        circulatingCurrentCoeffs = ICircfromIbiasPolyCoeffs[fr,betar,Icr];
        zerothOrderCoupling = mrqz * circulatingCurrentCoeffs[[1]] / phi00/2/Pi;
        firstOrderCoupling = mrqz * circulatingCurrentCoeffs[[2]] * biasCurrent/ phi00/2/Pi;
        secondOrderCoupling = mrqz * circulatingCurrentCoeffs[[3]] * biasCurrent^2/ phi00/2/Pi;
        {zerothOrderCoupling, firstOrderCoupling, secondOrderCoupling}
    ];

(*Tunable resonator vaccum fluctuation of circulating current in the SQUID.
includes 0, 1st and 2nd order fluctuations.
in units of Ampere. 
*)
TResVacuumFluctuations[cr_,Z0r_,lr_,betar_,Icr_,fr_] :=
    Module[ {omegar, effectiveL, biasCurrent,zerothOrderCoupling,
        firstOrderCoupling,secondOrderCoupling,circulatingCurrentCoeffs},
        effectiveL = rfsquidLeff[0,fr, betar,Icr,phi00];
        omegar = resonanceFrequency[cr,lr,Z0r,effectiveL]*2*Pi;
        biasCurrent = IResonator[omegar,cr,Z0r,lr,effectiveL];
        circulatingCurrentCoeffs = ICircfromIbiasPolyCoeffs[fr,betar,Icr];
        zerothOrderCoupling = circulatingCurrentCoeffs[[1]];
        firstOrderCoupling = circulatingCurrentCoeffs[[2]] * biasCurrent;
        secondOrderCoupling = circulatingCurrentCoeffs[[3]] * biasCurrent^2;
        {zerothOrderCoupling, firstOrderCoupling, secondOrderCoupling}
    ];

(*resonator frequency of a wavegide ended by rf-squid*)
resonatorFrequency[cr_,Z0r_,lr_,betar_,Icr_,fr_] :=
    Module[ {L, freq},
        L = rfsquidLeff[0, fr, betar, Icr, phi00];
        freq = resonanceFrequency[cr, lr, Z0r, L]
    ];

End[] (* End Private Context *)


EndPackage[]
