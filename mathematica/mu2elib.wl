(* ::Package:: *)

(* mu2elib.m  from mu2e_eft_6_2_20022.nb *)

(* Library of functions for Mu2E calculations *)

(* Some unit constants so we can use Mathematica utilities like IsotopeData *)
MeV = Quantity[1, "MeV"];

ElasticDir = Environment["MU2E_ELASTIC"];
If[ ElasticDir == $Failed ,
	Print["env MU2E_ELASTIC should be set to path to elastic density matrix directory"];
	Print["Obtain it from https://github.com/Berkeley-Electroweak-Physics/Elastic"]
	Throw["Missing Elastic Dir"]
];
If[ ! DirectoryQ[ElasticDir] ,
	Print["Elastic density matrix directory ", ElasticDir, ", specified by env MU2E_ELASTIC does not exist!"]
	Print["Obtain it from https://github.com/Berkeley-Electroweak-Physics/Elastic"]
	Throw["Missing Elastic Dir"]
];
Print["Found Elastic Dir at ", ElasticDir];
(* Construct file name for elastic results files *)
DMfn[s_]:= ElasticDir <> "/" <> s <> ".txt";

(*Directory of input DM files; [i,j,k] i is the nucleus (Ncode); j is the isotope; k is the interaction*)
DMstring[1,1,1]=DMfn["C12_ck"];
DMstring[1,2,1]=DMfn["C13_ck"];
DMstring[2,1,1]=DMfn["O16_bw"];
DMstring[2,1,2]=DMfn["O16_bw"];
DMstring[2,1,3]=DMfn["O16_bw"];
DMstring[2,1,4]=DMfn["O16_4hw"];
DMstring[2,2,1]=DMfn["O18_bw"];
DMstring[2,2,2]=DMfn["O18_usda"];
DMstring[2,2,3]=DMfn["O18_usdb"];
DMstring[2,2,4]=DMfn["O18_2hw"];
DMstring[3,1,1]=DMfn["F19_bw"];
DMstring[3,1,2]=DMfn["F19_usda"];
DMstring[3,1,3]=DMfn["F19_usdb"];
DMstring[4,1,1]=DMfn["Na23_bw"];
DMstring[4,1,2]=DMfn["Na23_usda"];
DMstring[4,1,3]=DMfn["Na23_usdb"];
DMstring[5,1,1]=DMfn["Al27_bw"];
DMstring[5,1,2]=DMfn["Al27_usda"];
DMstring[5,1,3]=DMfn["Al27_usdb"];
DMstring[6,1,1]=DMfn["Si28_bw"];
DMstring[6,1,2]=DMfn["Si28_usda"];
DMstring[6,1,3]=DMfn["Si28_usdb"];
DMstring[6,2,1]=DMfn["Si29_bw"];
DMstring[6,2,2]=DMfn["Si29_usda"];
DMstring[6,2,3]=DMfn["Si29_usdb"];
DMstring[6,3,1]=DMfn["Si30_bw"];
DMstring[6,3,2]=DMfn["Si30_usda"];
DMstring[6,3,3]=DMfn["Si30_usdb"];
DMstring[7,1,1]=DMfn["S32_bw"];
DMstring[7,1,2]=DMfn["S32_usda"];
DMstring[7,1,3]=DMfn["S32_usdb"];
DMstring[7,2,1]=DMfn["S33_bw"];
DMstring[7,2,2]=DMfn["S33_usda"];
DMstring[7,2,3]=DMfn["S33_usdb"];
DMstring[7,3,1]=DMfn["S34_bw"];
DMstring[7,3,2]=DMfn["S34_usda"];
DMstring[7,3,3]=DMfn["S34_usdb"];
DMstring[8,1,1]=DMfn["Ca40_kbp"];
DMstring[8,1,2]=DMfn["Ca40_kbp"];
DMstring[8,1,3]=DMfn["Ca40_kbp"];
DMstring[8,2,1]=DMfn["Ca42_kbp"];
DMstring[8,2,2]=DMfn["Ca42_GX1A"];
DMstring[8,2,3]=DMfn["Ca42_KB3G"];
DMstring[8,3,1]=DMfn["Ca44_kbp"];
DMstring[8,3,2]=DMfn["Ca44_GX1A"];
DMstring[8,3,3]=DMfn["Ca44_KB3G"];
DMstring[9,1,1]=DMfn["Ti46_kbp"];
DMstring[9,1,2]=DMfn["Ti46_GX1A"];
DMstring[9,1,3]=DMfn["Ti46_KB3G"];
DMstring[9,2,1]=DMfn["Ti47_kbp"];
DMstring[9,2,2]=DMfn["Ti47_GX1A"];
DMstring[9,2,3]=DMfn["Ti47_KB3G"];
DMstring[9,3,1]=DMfn["Ti48_kbp"];
DMstring[9,3,2]=DMfn["Ti48_GX1A"];
DMstring[9,3,3]=DMfn["Ti48_KB3G"];
DMstring[9,4,1]=DMfn["Ti49_kbp"];
DMstring[9,4,2]=DMfn["Ti49_GX1A"];
DMstring[9,4,3]=DMfn["Ti49_KB3G"];
DMstring[9,5,1]=DMfn["Ti50_kbp"];
DMstring[9,5,2]=DMfn["Ti50_GX1A"];
DMstring[9,5,3]=DMfn["Ti50_KB3G"];
DMstring[10,1,1]=DMfn["Fe54_kbp"];
DMstring[10,1,2]=DMfn["Fe54_kbp"];
DMstring[10,1,3]=DMfn["Fe54_kbp"];
DMstring[10,2,1]=DMfn["Fe56_kbp"];
DMstring[10,2,2]=DMfn["Fe56_GX1A"];
DMstring[10,2,3]=DMfn["Fe56_KB3G"];
DMstring[10,3,1]=DMfn["Fe57_kbp"];
DMstring[10,3,2]=DMfn["Fe57_GX1A"];
DMstring[10,3,3]=DMfn["Fe57_KB3G"];
DMstring[10,4,1]=DMfn["Fe58_kbp"];
DMstring[10,4,2]=DMfn["Fe58_GX1A"];
DMstring[10,4,3]=DMfn["Fe58_KB3G"];
DMstring[11,1,1]=DMfn["Cu63_GCN2850"];
DMstring[11,1,2]=DMfn["Cu63_jj44b"];
DMstring[11,1,3]=DMfn["Cu63_JUN45"];
DMstring[11,2,1]=DMfn["Cu65_GCN2850"];
DMstring[11,2,2]=DMfn["Cu65_jj44b"];
DMstring[11,2,3]=DMfn["Cu65_JUN45"];

(* Map element names/symbols/numbers to element numbers *)
elementmap = <|
   "1" -> 1, "C" -> 1, "Carbon" -> 1,
   "2" -> 2,  "O" -> 2, "Oxygen" -> 2,
   "3" -> 3,  "F" -> 3, "Fluorine" -> 3,
   "4" -> 4,  "Na" -> 4, "NA" -> 4, "Sodium" -> 4,
   "5" -> 5,  "Al" -> 5, "AL" -> 5, "Aluminum" -> 5,
   "6" -> 6,  "Si" -> 6, "SI" -> 6, "Silicon" -> 6,
   "7" -> 7,  "S" -> 7, "Sulfur" -> 7, "Sulpher" -> 7,
   "8" -> 8,  "Ca" -> 8, "CA" -> 8, "Calcium" -> 8,
   "9" -> 9,  "Ti" -> 9, "TI" -> 9, "Titanium" -> 9,
   "10" -> 10, "Fe" -> 10, "FE" -> 10, "Iron" -> 10,
   "11" -> 11,  "Cu" -> 11, "CU" -> 11, "Copper" -> 11
    |>;
elementhelp = 
"Use the number, element name, or symbol to select an element \n 1:Carbon(C), 2:Oxygen(O), 3:Fluorine(F), 4:Sodium(Na) \n 5:Aluminum(Al), 6:Silicon(Si), 7:Sulfur(S), 8:Calcium(Ca) \n 9:Titanium(Ti), 10:Iron(Fe), 11:Copper(Cu) \nYou may also optionally prefix the Name or symbol with 0 for an average isotope, or the specific isotope like 40Ca or 18Oxygen \n";

(* Map element number to standard name *)
Nucleus={"Carbon","Oxygen","Fluorine","Sodium","Aluminum","Silicon","Sulfur","Calcium","Titanium","Iron","Copper"};
(* Charge *)
ZvalA={6,8,9,11,13,14,16,20,22,26,29};
(* Adjustment to Z for muon wave function (eqn 26). *)
ZeffA={5.7030,7.4210,8.2298,9.8547,11.3086,12.0009,13.1839,15.6916,16.6562,18.6028,19.8563};
(* EbindA is the muon binding energy in 1S (depends on Z) ignore small A effect *)
EbindA={.099974,.177455,.224170,.3333684,.462953,.534646,0.692428,1.058480,1.261486,1.718157,2.088443};
(* momentum transfer to electron adjusted for coulomb interaction to get WF *)
qeffA={108.398,109.158,109.438,110.247,110.809,111.030,111.56,112.279,112.429,113.157,113.503};
(* Ratio of average muon Dirac components <f>/<g>  for each element *)
fgAvgA={-.01450,-.01895,-.02104,-.02542,-.02949,-.03141,-.03519,-.04221,-.04510,-.05148,-.05534};

(*Input ordinary muon capture rates in units of 1/sec; default values taken from the uncertainty-weighted average of the data compiled by Suzuki et al*)
MuCapRate={.0378,0.1024,0.2292,0.3773,0.6982,0.8652,1.351,2.525,2.592,4.411,5.673}*10^6;

(* Cataloged isotope options specifying A values for each element number *)
AOptions={{12,13},
{16,18},
{19},
{23},
{27},
{28,29,30},
{32,33,34},
{40,42,44},
{46,47,48,49,50},
{54,56,57,58},
{63,65}};

(* 
 * Try to get "Delta" from IsotopeData.  18O gets wrong answer that
 * doesn't match nuclear data cards or NuDat.  Other isotopes are correct.
 * Need to report to Wolfram.
 * Don't use this!
 *)
MassExcess = Table[ 
	Table[
		IsotopeData[{Nucleus[[i]], AOptions[[i]][[j]]}, "MassExcess"]/MeV, 
		 {j, 1, Length[AOptions[[i]]]}],
	{i, 1, Length[AOptions]}];

(*Input masses; mN is the nucleon mass;Delta is the mass defect in MeV from Nuclear Wallet Cards, 931.494028 is the 12C mass unit; Masses is the resulting array of nuclear masses in MeV*)
mN=939.;
Delta={{0.0,3.125},
{-4.737,-0.783},
{-1.487},
{-9.530},
{-17.197},
{-21.493,-21.895,-24.433},
{-26.015,-26.586,-29.931},
{-34.846,-38.547,-41.468},
{-44.128,-44.937,-48.493,-48.564,-51.432},
{-56.255,-60.607,-60.182,-62.155},
{-65.580,-67.264}};
Masses= Delta + 931.494028*AOptions;

(* Total Ground State J for each isotope *)
JS={{0,1/2},
       {0,0},
       {1/2},
       {3/2},
       {5/2},
       {0,1/2,0},
       {0,3/2,0},
       {0,0,0},
       {0,5/2,0,7/2,0},
       {0,0,1/2,0},
       {3/2,3/2}};
(* Isospin *)
TMTS={{{0,0},{1/2,-1/2}},
         {{0,0},{1,-1}},
         {{1/2,-1/2}},
         {{1/2,-1/2}},
         {{1/2,-1/2}},
         {{0,0},{1/2,-1/2},{1,-1}}, 
         {{0,0},{1/2,-1/2},{1,-1}},
         {{0,0},{1,-1},{2,-2}},
         {{1,-1},{3/2,-3/2},{2,-2},{5/2,-5/2},{3,-3}},
         {{1,-1},{2,-2},{5/2,-5/2},{3,-3}},
         {{5/2,-5/2},{7/2,-7/2}}};
(*  Abundances[[Ncode]][[isotopeidx]] gives fractional abundance for element Ncode *)
Abundances={{.9893,.0107},{.99757,.00205},{1.},{1.},{1.},{.92223,.04685,.03092},{.9499,.0075,.0425},{.9694,.00647,.0209},{.0825,.0744,.7372,.0541,.0518},{.05845,.91754,.02119,.00282},{0.6915,.3085}};

(* Number of Isotopes cataloged for each element number *)
NumAs = Table[Length[AOptions[[i]]], {i, 1, Length[AOptions]}];

(* Used for plots, which go to a file from the testing script *)
(* for script file ca40.yaml, sbf will be "ca40" *)
scriptBaseFile = None;
setBaseFile[sbf_] := Global`scriptBaseFile = sbf;

setElementData[Ncode_] := Module[{mbar, jj},
	Print["Setting up for ", Nucleus[[Ncode]]];
	Global`Zval = ZvalA[[Ncode]];
	Global`qeff = qeffA[[Ncode]];
	Global`fgAvg = fgAvgA[[Ncode]];
	Global`Ebind=EbindA[[Ncode]];
	Global`Zeff=ZeffA[[Ncode]];
	Global`RZ2= (Zeff  alpha  mu)^3/Pi;
	Global`mumass=105.6583755; (* Muon Mass *)
	Global`emass=0.511;   (* FIXME:   0.51099895000 MeV *)
	Global`alpha=1/(137.035999082);
	Global`hbar=6.5821*10^(-22);
	Global`Nmass=939.;     (* FIXME:  Average nucleon mass? *)
						   (*  p 938.27208816 n 939.56542052 *)

	(* this is the weak scale in MeV, sqrt[1/(Sqrt[2]GF)], 
	 * our normalization for the EFT's LECs
	 *)
	Global`mv=246200;

	Print["qeff  = ", qeff];
	Print["fgAvg = ", fgAvg];
	Print["Ebind = ", Ebind];
	Print["Zeff  = ", Zeff];

	(* Calculate Mbar, the average nuclear mass *)
	mbar=0;
	Do[ mbar=mbar+ AbNorm[jj] Masses[[Ncode,jj]],{jj,1,NumAs[[Ncode]]}];
	Print["Average Nuclear Mass = ",mbar," MeV"];
	Global`Mbar = mbar;

	Global`mu=mumass Mbar/(mumass+Mbar);
	Global`qval=Sqrt[Mbar/(mumass+Mbar) ((mumass-Ebind)^2-emass^2)];

	Print["  Reduced mass for nucleus+muon = ", mu];
	Print["  Z = ",Zval,"   Zeff = ",Zeff,"   RZ2 = ",RZ2];
	Print["  Momentum transfer to electron = ",qval,"MeV     qeff = ",qeff," MeV"];
	Print["  muon binding energy = ",Ebind," MeV"];
	Print["  Ratio of average muon Dirac components <f>/<g> = ", fgAvg];
];

Astring=" AOption=0 yields average over isotopes, \n    weighted by respective abundances, \n    including all isotopes with abundances>0.2%, \n    and renormalized to 1 to correct for \n    missing isotopes of abundance <0.2% \n \n AOption != 0 allows an individual isotope to be  \n    selected with unit abundance assumed \n    as in an isotopically pure target \n \n";

(* ai can be None, then ask *)
getAOption[ai_] := Module[{s, ii, jj, asum, abnorm, abar, pos},
	s = Astring;
	Do[ s=StringJoin[s," AOption=",ToString[ijk]," A=",ToString[AOptions[[Ncode,ijk]]],"  \n \n"],{ijk,1,NumAs[[Ncode]]}];
	If[SameQ[ai, None],
		AString=StringJoin[s,"Input integer AOption now"];
		Global`AOption=Input[s];
	,
		If[ai == 0,
			Global`AOption=0;
		,
			(* Search for isotope in AOptions[[Ncode]] to get isotope position *)
			pos = Position[AOptions[[Ncode]], ai];
			If[Length[pos] != 1, Throw["Unknown Isotope " <> ToString[ai] <> " for " Nucleus[[Ncode]] ]];
			AOption = pos[[1]][[1]];
		];
	];

	If[AOption==0,
		asum=0;
		Do[asum=asum+Abundances[[Ncode,ii]],{ii,1,NumAs[[Ncode]]}];
		Do[Global`AbNorm[ii]=Abundances[[Ncode,ii]]/asum,{ii,1,NumAs[[Ncode]]}]
	,
		Do[Global`AbNorm[ii]=0,{ii,1,NumAs[[Ncode]]}];
		Global`AbNorm[AOption]=1
	];
	(*Calculate Abar, the average A*) 
	abar=0; 
	Do[
		abar=abar+ AbNorm[jj] AOptions[[Ncode,jj]]
	,{jj,1,NumAs[[Ncode]]}];
	Global`Abar = abar;
];
reportAOption[] := Module[{},
	Do[
		Print["Target is ", AOptions[[Ncode,ii]], Nucleus[[Ncode]]," 
			at abundance ", AbNorm[ii]]
	,{ii,1,NumAs[[Ncode]]}];
	Print["Weighted Average A = ",Abar];
];

setElement[ae_] := Module[{e, s},
	s = StringTrim[ae];
	e = StringCases[s, RegularExpression["[a-zA-Z]+$"]];
	If[Length[e] != 1, Throw["No element specification"]];
	e = e[[1]];
	If[! KeyExistsQ[elementmap, e],
		Print["Unknown element ", e];
		Throw["Unknown element " <> e];
	];
	Global`Ncode = elementmap[e];
	Print["Target nucleus is ",Nucleus[[Ncode]]];
	a = StringCases[s, RegularExpression["^[0-9]+"]];
	If[Length[a] == 1,
		getAOption[FromDigits[a[[1]]]]
	,
		getAOption["0"] (* default to average of isotopes *)
	];
	setElementData[Ncode];
	reportAOption[];
];

(* Ask User for element selection *)
getElement[] := Module[{e, s, a},
	s = StringTrim[InputString[elementhelp]];
	e = StringCases[s, RegularExpression["[a-zA-Z]+$"]];
	If[Length[e] != 1, Throw["No element specification"]];
	e = e[[1]];
	If[! KeyExistsQ[elementmap, e],
		Print["Unknown element ", e];
		Throw["Unknown element " <> e];
	];
	Global`Ncode = elementmap[e];
	Print["Target nucleus is ",Nucleus[[Ncode]]];
	a = StringCases[s, RegularExpression["^[0-9]+"]];
	If[Length[a] == 1,
		getAOption[FromDigits[a[[1]]]]
	,
		(* Have to ask *)
		getAOption[None]
	];
	setElementData[Ncode];
	reportAOption[];
];

(*
 * Set Interaction with string name for interaction
 * Note that interaction index depends on Ncode, they
 * are not global.
 *)
setInteraction[aints_] := Module[{ism},
	ints = ToLowerCase[aints];
	Print["Looking for interaction ", ints, " for nucleus ", Nucleus[[Ncode]]];
	ism = -1; (* to catch error *)
	If[Ncode == 1,
		ism = 1;
		Print["Interaction used is Cohen and Kurath"];
	];
	If[Ncode == 2,
		tab = <| "bw" -> 1, "usda" -> 2, "usdb" -> 3, "2hw" -> 4, "4hw" -> 4 |>;
		ism = tab[ints];
		If[ism==1,Print["ISM=1:  Brown-Wildenthal"]];
		If[ism==2,Print["ISM=2:  USDA"]];
		If[ism==3,Print["ISM=3:  USDB"]];
		If[ism==4,Print["ISM=4: 4hw(16O) and 2hw(18O) HJ wave functions"]];
		If[ism<4,Print["16O treated as a closed shell if ISM=1,2,3"]];
	];
	If[2<Ncode<8,
		tab = <| "bw" -> 1, "usda" -> 2, "usdb" -> 3 |>;
		ism = tab[ints];
		If[ism==1,Print["ISM=1:  Brown-Wildenthal"]];
		If[ism==2,Print["ISM=2:  USDA"]];
		If[ism==3,Print["ISM=3:  USDB"]];
	];
	If[7<Ncode<11,
		tab = <| "kbp" -> 1, "gx1a" -> 2, "kb3g" -> 3 |>;
		ism = tab[ints];
		If[ism==1,Print["ISM=1:  KBP"]];
		If[ism==2,Print["ISM=2:  GX1A"]];
		If[ism==3,Print["ISM=3:  KB3G"]];
	];
	If[Ncode==11,
		Print["Ncode is 11"];
		tab = <| "gcn2850" -> 1, "jj44b" -> 2, "jun45" -> 3 |>;
		ism = tab[ints];
		If[ism==1,Print["ISM=1:  GCN2850"]];
		If[ism==2,Print["ISM=2:  jj44b"]];
		If[ism==3,Print["ISM=3:  JUN45"]];
	];
	If[ism < 0, Throw["Unknown interaction " <> aints <> " for nucleus " <> Nucleus[Ncode]]];
	Global`ISM = ism;
	Print["Set ISM=", ISM];
];

(* 
 * Interaction choices depend on Ncode
 * Orig code reused indices for diff interactions based on Ncode
 *)
getInteraction[]:= Module[{ism},
	ism=1;
	If[Ncode==1,Print["Interaction used is Cohen and Kurath"]];
	If[Ncode==2,
		ism=Input[" Interaction choices are \n   ISM=1: Brown-Wildenthal \n   ISM=2: USDA \n   ISM=3: USDB \n   ISM=4: 4hw(16O) and 2hw(18O) HJ wave functions \n \n Input choice now"];
		If[ism==1,Print["ISM=1:  Brown-Wildenthal"]];
		If[ism==2,Print["ISM=2:  USDA"]];
		If[ism==3,Print["ISM=3:  USDB"]];
		If[ism==4,Print["ISM=4: 4hw(16O) and 2hw(18O) HJ wave functions"]];
		If[ism<4,Print["16O treated as a closed shell if ISM=1,2,3"]];
	];
	If[2<Ncode<8,
		ism=Input[" Interaction choices are \n   ISM=1: Brown-Wildenthal \n   ISM=2: USDA \n   ISM=3: USDB  \n \n Input choice now"];
		If[ism==1,Print["ISM=1:  Brown-Wildenthal"]];
		If[ism==2,Print["ISM=2:  USDA"]];
		If[ism==3,Print["ISM=3:  USDB"]];
	];

	If[7<Ncode<11,
		ism=Input[" Interaction choices are \n   ISM=1: KBP \n   ISM=2: GX1A \n   ISM=3: KB3G \n 40Ca treated as a closed shell \n \n  Input choice now"];
		If[ism==1,Print["ISM=1:  KBP"]];
		If[ism==2,Print["ISM=2:  GX1A"]];
		If[ism==3,Print["ISM=3:  KB3G"]];
		If[Ncode==8,
			Print["40Ca always treated as a closed shell, independent of ISM"]
		];
	];

	If[Ncode==11,
		ism=Input[" Interaction choices are \n   ISM=1: GCN2850 \n   ISM=2: jj44b \n   ISM=3: JUN45 \n \n Input choice now"];
		If[ism==1,Print["ISM=1:  GCN2850"]];
		If[ism==2,Print["ISM=2:  jj44b"]];
		If[ism==3,Print["ISM=3:  JUN45"]];
	];
	Global`ISM = ism;
];

(* Set Default Oscillator parameter *)
setDefaultB[] := Module[{bb},
	(*Input oscillator parameter*)
	bb=Sqrt[41.467/(45Abar^(-1/3) -25Abar^(-2/3))];
	If[Ncode==1,bb=1.77,True];
	If[Ncode==3,bb=1.833,True];
	Global`bfm = bb;
];

(* b is scaled by hbar c to yield MeV^-1 *)
scaleB[] := Module[{},
	(* FIXME:  Should not be hardcoded \[HBar]c -  slightly off  *)
	Global`b = bfm / 197.3269718; 
];

setB[abfm_]:= Module[{},
	If[abfm == 0, 
		setDefaultB[]
	,
		bfm = abfm
	];
	scaleB[];
	Print["Oscillator parameter = ",b," fermis"];
];

(* Set Oscillator parameter *)
getB[] := Module[{as, b0fm},
	setDefaultB[];
	as=StringJoin[" Internally set, using average A \n Oscillator parameter b = ",ToString[bfm]," fermis \n \n"," If you want to manually reset b, enter value now \n Otherwise enter 0 \n \n"];
	b0fm=Input[as]; 
	If[b0fm==0,True,bfm=b0fm];
	Print["Oscillator parameter = ",bfm," fermis"];
	scaleB[];
];

readDM[] := Module[{DMin, ii, jj, kk, Lc, LcMax},
	Print["Reading Interaction Files for selected Isotopes"];
	Clear[Global`DM, Global`JF, Global`TF, Global`JI, Global`TI, 
	      Global`J0, Global`T0, Global`NB, Global`JB, Global`NK, Global`JK, Global`DmVal];
	Do[  (* For each isotope, can be multiple with Aoption=0 *)
		DMin=ReadList[DMstring[Ncode,ii,ISM],Number,RecordLists->True];
		Print["For ", AOptions[[Ncode]][[ii]], Nucleus[[Ncode]]];
		LcMax=Length[DMin];
		jj=0;
		Do[
			If[Length[DMin[[Lc]]]==6,
			    jj=jj+1;
				Global`NumDMs[ii]=jj;
				kk=0;
				Global`JF[ii]=DMin[[Lc,1]]/2;
				Global`TF[ii]=DMin[[Lc,2]]/2;
				Global`JI[ii]=DMin[[Lc,3]]/2;
				Global`TI[ii]=DMin[[Lc,4]]/2;
				Global`J0[ii,jj]=DMin[[Lc,5]]/2;
				Global`T0[ii,jj]=DMin[[Lc,6]]/2;
				If[JF[ii]==JS[[Ncode,ii]],True,Abort[]];
				If[TF[ii]==TMTS[[Ncode,ii,1]],True,Abort[]];
				If[JI[ii]==JS[[Ncode,ii]],True,Abort[]];
				If[TI[ii]==TMTS[[Ncode,ii,1]],True,Abort[]];
			];
			If[Length[DMin[[Lc]]]==5,
				kk=kk+1;
				Global`LDM[ii,jj]=kk;
				Global`NB[ii,jj,kk]=DMin[[Lc,1]];
				Global`JB[ii,jj,kk]=DMin[[Lc,2]]/2;
				Global`NK[ii,jj,kk]=DMin[[Lc,3]];
				Global`JK[ii,jj,kk]=DMin[[Lc,4]]/2;
				Global`DMval[ii,jj,kk]=DMin[[Lc,5]];
			]
		, {Lc,1,LcMax}];
		Do[
			Global`DM[ii,jj]=Table[{{NB[ii,jj,kk],JB[ii,jj,kk]},{NK[ii,jj,kk],JK[ii,jj,kk]},DMval[ii,jj,kk]} ,{kk,1,LDM[ii,jj]}];
		,{jj,1,NumDMs[ii]}];
		Print["  Number of DMs = ",NumDMs[ii]];
		Do[Print["   Lengths =",LDM[ii,jj]],{jj,1,NumDMs[ii]}];
		Do[
			Print["  JF = ",JF[ii],"  TF = ",TF[ii],
			   "  JI = ",JI[ii],"  TI = ",TI[ii],
			   "  J0 = ",J0[ii,jj],"  T0= ",T0[ii,jj]];
			Print["    ", MatrixForm[DM[ii,jj]]]
		,{jj,1,NumDMs[ii]}];
	,{ii,1,NumAs[[Ncode]]}];
	Print["Finished with Interactions"];
];


(* Test interactions by computing simple properties *)
calcSumRules[] := Module[{},
	Print["Computing A (T0=0) and Z-N (T0=1) Sum Rules"];
	(*Check J=0 T=0,1 DMs by evaluating the A and Z-N sum rules*)
	Do[
		Print[AOptions[[Ncode]][[ii]], Nucleus[[Ncode]]];
		Do[
			If[J0[ii,jj]==0,
				Valu=Sqrt[4 Pi/(2JI[ii]+1)]FM[ii,jj,0.00000001];
				Print["  J0 = ",J0[ii,jj],"  T0 = ",T0[ii,jj],"  Charge = ",Valu];
			]
		, {jj,1,NumDMs[ii]}];
	, {ii,1,NumAs[[Ncode]]}];
];

calcDecayRate[] := Module[{qm, y},
	qm=qeff/Nmass;
	y=(qeff b/2)^2;
	Global`Decayrate=Re[1/(mv^4 2  Pi) RZ2  qeff^2/(1+qval/Mbar)Heff[qm,y,cs,bs]/hbar];
	Print["  Decay rate = ",Decayrate,"/sec"];
];

setMCR[v_] := Module[{},
	If[v == 0,
		Global`MCR=MuCapRate[[Ncode]];
		Print["  Ordinary Muon Capture rate = ",MCR," /sec (default)"];
	,
		Global`MCR = v;
		Print["  Ordinary Muon Capture rate = ",MCR," /sec (override)"];
	];
];

getMCR[] := Module[{},
	Global`MCR=MuCapRate[[Ncode]];
	Astring=StringJoin[" Ordinary Muon capture rate from weighted average of Suzuki et al. data ",ToString[MCR]," /sec \n \n"," If you want to manually reset MCR, enter value now \n Otherwise enter 0 \n \n"];
	MCR0=Input[Astring];
	If[MCR0!=0,
		Global`MCR=MCR0
		Print["  Ordinary Muon Capture rate = ",MCR," /sec (override)"];
	,
		Print["  Ordinary Muon Capture rate = ",MCR," /sec (default)"];
	];
];

setML[amL_] := Module[{},
	Global`mL = amL;
	Print["mL = ", mL, " MeV"];
];

getML[] := Module[{as, v, found},
	as = " If you want to use relativistic rather than standard nonrelativistic operators, enter a nonzero positive value now for the leptonic scale mL in MeV \n Otherwise enter 0 \n \n";
	found = False;
	While[!found,
		v = InputString[as];
		v = StringTrim[v];
		If[v != "",
			v = ToExpression[v];
			If[NumericQ[v],
				found = True;
				Global`mL = v;
			,
				Print["Non numeric value ", v, " - try again"];
			];
		,
			Print["Empty input, try again"];
		];
	];
	Print["mL = ", mL, " MeV"];
];

setiLower[v_] := Module[{},
	If[ mL > 0,
		Global`iLower = v;
	,
		Global`iLower = 0;
	];
	If[iLower > 0, Print["Lower-component muon contributions included"]];
];

(* Decide if lower-component muon contributions are to be calculated *)
getILower[] := Module[{as},
	Global`iLower = 0;
	If[ mL > 0,
		as = "Enter 1 to include muon lower-component contributions \nEnter 0 to ignore muon lower-component contributions \n";
		Global`iLower = Input[as];
	];
];

clearCDBs[]:= Module[{i},
	Clear[Global`cs, Global`ds, Global`bs];
	For[i = 1, i <= 32, i++,
		Global`ds[i,0] = 0;
		Global`ds[i,1] = 0;
	];
	For[i = 1, i <= 16, i++,
		Global`cs[i,0] = 0;
		Global`cs[i,1] = 0;
		Global`bs[i,0] = 0;
		Global`bs[i,1] = 0;
	];
];

reportCS[] := Module[{i},
	Print["Nonzero cs values"];
	Do[
		If[cs[i,0] != 0 || cs[i,1] != 0,
			Print["   i ",i,"  {",cs[i,0],",",cs[i,1],"}"];
		];
	,{i,1,16}];
];

setCS[cslist_] := Module[{ci, idx, c0, c1},
	If[ mL > 0, Throw["cs not used with relativistic coefficients"]];
	For[ci = 1, ci <= Length[cslist], ci++,
		{idx, c0, c1} = cslist[[ci]];
		If[idx < 1 || idx > 16, Throw["cs index out of range 1..16"]];
		Global`cs[idx, 0] = c0;
		Global`cs[idx, 1] = c1;
	];
	reportCS[];
];

readCS[] := Module[{as, A, i},
	If[mL==0,
		as ="Input nonzero nonrelativistic ci coefficients as {i,ci[0],ci[1]} successively, \n for i=1,3,4,5,....16 as needed \n Inputting i=0 terminates the input \n Input can be numerical or symbolic \n The first coupling is isoscalar, the second is isovector \n Thus {i,ci[0],0} is a symbolic S.D. isoscalar coupling of strength ci[0] \n {i,0,ci[1]} is a symbolic S.D. isovector coupling of strength ci[1] \n {i,cip/2,cip/2} is a symbolic S.D. proton only coupling of stength cip  \n {i,cin/2,-cin/2} is a symbolic S.D. neutron only coupling of strength cin \n \n Input desired {i,ci[0],ci[1]} now";
		A = Input[as];
		i=A[[1]];
		Global`cs[i,0]=A[[2]];
		Global`cs[i,1]=A[[3]];
		Do[
			A=Input["Input {i,ci[0],ci[1]} now"];
			i=A[[1]];
			If[i==0,Break[]];
			Global`cs[i,0]=A[[2]];
			Global`cs[i,1]=A[[3]]
		,{jj,2,16}];
		reportCS[];
	];
];

reportDS[] := Module[{i},
	Print["nonzero ds values for scalar-mediated operators"];
	Do[If[ds[i,0]==0&&ds[i,1]==0,True,Print["   d(",i,"):  {",ds[i,0],",",ds[i,1],"}"]],{i,1,4}];
	Print["nonzero ds values for vector-mediated operators"];
	Do[If[ds[i,0]==0&&ds[i,1]==0,True,Print["   d(",i,"):  {",ds[i,0],",",ds[i,1],"}"]],{i,5,20}];
	Print["nonzero ds values for tensor-mediated operators"];
	Do[If[ds[i,0]==0&&ds[i,1]==0,True,Print["   d(",i,"):  {",ds[i,0],",",ds[i,1],"}"]],{i,21,32}]
];

(* Set relativistic coefficients *)
setDS[dslist_] := Module[{di, idx, d0, d1},
	If[mL == 0, Throw["ds coefficients don't make sense with mL == 0"]];
	For[di = 1, di <= Length[dslist], di++,
		{idx, d0, d1} = dslist[[di]];
		If[idx < 1 || idx > 32, Throw["ds index out of range 1..32"]];
		Global`ds[idx, 0] = d0;
		Global`ds[idx, 1] = d1;
	];
	reportDS[];
];

(* Read relativistic ds coefficients *)
readDS[] := Module[{as, i, A, jj},
	as = "Input nonzero relativistic di coefficients as {i,di[0],di[1]}
successively, \n for i=1,3,4,5,....32 \n Inputting i=0 terminates the input \n Input can be numerical or symbolic \n The first coupling is isoscalar, the second is isovector \n Thus {di[0],0} is a symbolic S.D. isoscalar coupling of strength ci[0] \n {0,di[1]} is a symbolic S.D. isovector coupling of strength di[1] \n {dip/2,dip/2} is a symbolic S.D. proton only coupling of stength dip  \n {din/2,-din/2} is a symbolic S.D. neutron only coupling of strength din \n \n
Input desired {i,di[0],di[1]} now";
	If[mL>0,
		A=Input[as];
		i=A[[1]];
		Global`ds[i,0]=A[[2]];
		Global`ds[i,1]=A[[3]];
		Do[
			A=Input["Input {i,di[0],di[1]} now"];
			i=A[[1]];
			If[i==0,Break[]];
			Global`ds[i,0]=A[[2]];
			Global`ds[i,1]=A[[3]]
		,{jj,2,32}];
		reportDS[];
	];
];

updateCS[] := Module[{},
	If[mL <= 0, Return[]];
	Global`cs[1,0]  += ds[1,0];
	Global`cs[1,1]  += ds[1,1];
	Global`cs[10,0] += ds[2,0]*qval/(2*mN);
	Global`cs[10,1] += ds[2,1]*qval/(2*mN);
	Global`cs[11,0] += -ds[3,0];
	Global`cs[11,1] += -ds[3,1];
	Global`cs[6,0]  += -ds[4,0]*qval/(2*mN);
	Global`cs[6,1]  += -ds[4,1]*qval/(2*mN);
	Global`cs[1,0]  += ds[5,0];
	Global`cs[1,1]  += ds[5,1];
	Global`cs[2,0]  += I*ds[5,0];
	Global`cs[2,1]  += I*ds[5,1];
	Global`cs[5,0]  += -ds[5,0];
	Global`cs[5,1]  += -ds[5,1];
	Global`cs[4,0]  += -ds[5,0] qval/(2 mN);
	Global`cs[4,1]  += -ds[5,1] qval/(2mN);
	Global`cs[6,0]  += -ds[5,0] qval/(2 mN);
	Global`cs[6,1]  += -ds[5,1] qval/(2mN);
	Global`cs[4,0]  += ds[6,0] qval/mN;
	Global`cs[4,1]  += ds[6,1] qval/mN;
	Global`cs[6,0]  += ds[6,0] qval/ mN;
	Global`cs[6,1]  += ds[6,1] qval/mN;

	Global`cs[7,0]  += ds[7,0];
	Global`cs[7,1]  += ds[7,1];
	Global`cs[10,0] += I*ds[7,0];
	Global`cs[10,1] += I*ds[7,1] ;
	Global`cs[9,0]  += -ds[7,0];
	Global`cs[9,1]  += -ds[7,1];
	Global`cs[10,0] += -ds[8,0]*qval/mN;
	Global`cs[10,1] += -ds[8,1]*qval/mN;
	Global`cs[1,0]  += -ds[9,0]*qval/mL;
	Global`cs[1,1]  += -ds[9,1]*qval/mL;
	Global`cs[5,0]  += -ds[9,0]*qval/mL;
	Global`cs[5,1]  += -ds[9,1]*qval/mL;
	Global`cs[4,0]  += -ds[9,0]*qval^2/(2*mN*mL);
	Global`cs[4,1]  += -ds[9,1]*qval^2/(2*mN*mL);
	Global`cs[6,0]  += -ds[9,0]*qval^2/(2*mN*mL);
	Global`cs[6,1]  += -ds[9,1]*qval^2/(2*mN*mL);
	Global`cs[4,0]  += ds[10,0]*qval^2/(mN*mL);
	Global`cs[4,1]  += ds[10,1]*qval^2/(mN*mL);
	Global`cs[6,0]  += ds[10,0]*qval^2/(mN*mL);
	Global`cs[6,1]  += ds[10,1]*qval^2/(mN*mL);
	Global`cs[7,0]  += -ds[11,0]*qval/mL;
	Global`cs[7,1]  += -ds[11,1]*qval/mL;
	Global`cs[9,0]  += -ds[11,0]*qval/mL;
	Global`cs[9,1]  += -ds[11,1]*qval/mL;
	Global`cs[10,0] += ds[12,0]*qval^2/(mN*mL);
	Global`cs[10,1] += ds[12,1]*qval^2/(mN*mL);

	Global`cs[11,0] += -I*ds[13,0];
	Global`cs[11,1] += -I*ds[13,1];
	Global`cs[8,0]  += -ds[13,0];
	Global`cs[8,1]  += -ds[13,1];
	Global`cs[9,0]  += -ds[13,0]*qval/(2*mN);
	Global`cs[9,1]  += -ds[13,1]*qval/(2*mN);
	Global`cs[9,0]  += ds[14,0]*qval/mN;
	Global`cs[9,1]  += ds[14,1]*qval/mN;
	Global`cs[14,0] += -I*ds[15,0];
	Global`cs[14,1] += -I*ds[15,1];
	Global`cs[4,0]  += -ds[15,0];
	Global`cs[4,1]  += -ds[15,1];
	Global`cs[6,0]  += I*ds[16,0]*qval/mN;
	Global`cs[6,1]  += I*ds[16,1]*qval/mN;
	Global`cs[11,0] += -ds[17,0]*qval/mL;
	Global`cs[11,1] += -ds[17,1]*qval/mL;
	Global`cs[8,0]  += -I*ds[17,0]*qval/mL;
	Global`cs[8,1]  += -I*ds[17,1]*qval/mL;
	Global`cs[9,0]  += -I*ds[17,0]*qval^2/(2 mN mL);
	Global`cs[9,1]  += -I*ds[17,1]*qval^2/(2 mN mL);
	Global`cs[16,0] += -I*ds[17,0]*qval/mL;
	Global`cs[16,1] += -I*ds[17,1]*qval/mL;
	Global`cs[9,0]  += I*ds[18,0]*qval^2/(mN mL);
	Global`cs[9,1]  += I*ds[18,1]*qval^2/(mN mL);
	Global`cs[14,0] += -ds[19,0]*qval/mL;
	Global`cs[14,1] += -ds[19,1]*qval/mL;
	Global`cs[4,0]  += -I*ds[19,0]*qval/mL;
	Global`cs[4,1]  += -I*ds[19,1]*qval/mL;
	Global`cs[6,0]  += -I*ds[19,0]*qval/mL;
	Global`cs[6,1]  += -I*ds[19,1]*qval/mL;
	Global`cs[6,0]  +=  ds[20,0]*qval^2/(mN mL);
	Global`cs[6,1]  +=  ds[20,1]*qval^2/(mN mL);

	(* New Entries *)
	Global`cs[1,0]  += -ds[21,0]*qval/mN- ds[22,0]*qval/mN+4*ds[23,0]*qval/mN+ds[29,0]*qval^2/(2*mL*mN);
	Global`cs[1,1]  += -ds[21,1]*qval/mN- ds[22,1]*qval/mN+4*ds[23,1]*qval/mN+ds[29,1]*qval^2/(2*mL*mN);
	Global`cs[3,0]  += -2*ds[21,0]+ ds[29,0]*qval/mL;
	Global`cs[3,1]  += -2*ds[21,1]+ ds[29,1]*qval/mL;
	Global`cs[4,0]  +=  2*ds[21,0]-4*I*ds[24,0]*qval/mN+ds[29,0]*qval/mL;
	Global`cs[4,1]  +=  2*ds[21,1]-4*I*ds[24,1]*qval/mN+ds[29,1]*qval/mL;
	Global`cs[6,0]  += -4*I*ds[24,0]*qval/mN+ ds[29,0]*qval/mL + 4*I*ds[31,0]*qval/mL;
	Global`cs[6,1]  += -4*I*ds[24,1]*qval/mN+ ds[29,1]*qval/mL + 4*I*ds[31,1]*qval/mL;
	Global`cs[9,0]  += -2*I*ds[25,0]-4*ds[28,0]*qval/mN + ds[30,0]*qval/mL;
	Global`cs[9,1]  += -2*I*ds[25,1]-4*ds[28,1]*qval/mN + ds[30,1]*qval/mL;
	Global`cs[10,0] += -2*ds[25,0]-4*ds[32,0]*qval/mL;
	Global`cs[10,1] += -2*ds[25,1]-4*ds[32,1]*qval/mL;
	Global`cs[11,0] +=  ds[25,0]*qval/mN+ ds[26,0]*qval/mN-4*ds[27,0]*qval/mN-I*ds[30,0]*qval^2/(2*mL*mN);
	Global`cs[11,1] +=  ds[25,1]*qval/mN+ ds[26,1]*qval/mN-4*ds[27,1]*qval/mN-I*ds[30,1]*qval^2/(2*mL*mN);
	Global`cs[12,0] += -2*ds[25,0] + 4*ds[32,0]*qval/mL;
	Global`cs[12,1] += -2*ds[25,1] + 4*ds[32,1]*qval/mL;
	Global`cs[13,0] += -2*I*ds[21,0] + 4*ds[31,0]*qval/mL;
	Global`cs[13,1] += -2*I*ds[21,1] + 4*ds[31,1]*qval/mL;
	Global`cs[15,0] += -I*ds[30,0]*qval/mL + 4*ds[32,0]*qval/mL;
	Global`cs[15,1] += -I*ds[30,1]*qval/mL + 4*ds[32,1]*qval/mL;

	Print["nonzero computed cs values"];
	Do[If[cs[i,0]==0&&cs[i,1]==0,True,
		Print["   i ",i,"  {",cs[i,0],",",cs[i,1],"}"]],{i,1,16}];
];

(* Compute bs values from relativistic ds values *)
updateBS[] := Module[{},
	If[mL <= 0 || iLower != 1, 
		Print["All bs values are 0"];
		Return[]
	];
	(* put in mult * to make python match easier *)
	Global`bs[2,0]=I*ds[1,0] - I*ds[5,0] - I*ds[9,0]*qval/mL;
	Global`bs[2,1]=I*ds[1,1] - I*ds[5,1] - I*ds[9,1]*qval/mL;
	Global`bs[3,0]= -ds[1,0] + ds[5,0]+ ds[9,0]*qval/mL;
	Global`bs[3,1]= -ds[1,1] + ds[5,1]+ ds[9,1]*qval/mL;
	Global`bs[5,0]= ds[15,0] - I*ds[19,0]*qval/mL + 2*ds[21,0] - ds[29,0]*qval/mL;
	Global`bs[5,1]= ds[15,1] - I*ds[19,1]*qval/mL + 2*ds[21,1] - ds[29,1]*qval/mL;
	Global`bs[7,0]= I*ds[3,0] + ds[13,0] + I*ds[17,0]*qval/mL;
	Global`bs[7,1]= I*ds[3,1] + ds[13,1] + I*ds[17,1]*qval/mL;
	Global`bs[8,0]= -ds[7,0] + ds[11,0]*qval/mL + 2*I*ds[25,0] + ds[30,0]*qval/mL;
	Global`bs[8,1]= -ds[7,1] + ds[11,1]*qval/mL + 2*I*ds[25,1] + ds[30,1]*qval/mL;
	Global`bs[12,0]= -I*ds[7,0] + I*ds[11,0]*qval/mL - 2*ds[25,0] + I*ds[30,0]*qval/mL;
	Global`bs[12,1]= -I*ds[7,1] + I*ds[11,1]*qval/mL - 2*ds[25,1] + I*ds[30,1]*qval/mL;
	Global`bs[13,0]= I*ds[15,0] + ds[19,0]*qval/mL + 2*I*ds[21,0] - I*ds[29,0]*qval/mL;
	Global`bs[13,1]= I*ds[15,1] + ds[19,1]*qval/mL + 2*I*ds[21,1] - I*ds[29,1]*qval/mL;
	Global`bs[14,0]= I*ds[15,0] + 2*I*ds[21,0] - 4*ds[31,0]*qval/mL;
	Global`bs[14,1]= I*ds[15,1] + 2*I*ds[21,1] - 4*ds[31,1]*qval/mL;
	Global`bs[15,0]= I*ds[11,0]*qval/mL + I*ds[30,0]*qval/mL - 4*ds[32,0]*qval/mL;
	Global`bs[15,1]= I*ds[11,1]*qval/mL + I*ds[30,1]*qval/mL - 4*ds[32,1]*qval/mL;
	Global`bs[16,0]= ds[11,0]*qval/mL + ds[30,0]*qval/mL + 4*I*ds[32,0]*qval/mL;
	Global`bs[16,1]= ds[11,1]*qval/mL + ds[30,1]*qval/mL + 4*I*ds[32,1]*qval/mL;

	Print["Nonzero computed bs values:"];
	Do[
		If[bs[i,0]==0&&bs[i,1]==0,
			True
		,
			Print["   i ",i,"  {",bs[i,0],",",bs[i,1],"}"]
		]
	,{i,1,16}];
];


(* The Seven Operators package *)
minus=-1;
JNorm[j_]:=2j+1; (* number of z projection states for angular mom j *)
QNorm[j_]:=Sqrt[2 j+1];
(* Triangular constraint on angular momentum addition *)
Triangular[x_,y_,z_]:=Abs[x-y]<=z<=x+y;
(* condition on the z-projection of the momentum, m *)
RangeM[j_,m_]:=-Abs[j]<=m<=Abs[j];
(* Conditions for a nonzero 6J symbol *)
SixJConditionTriad[j1_,j2_,j3_]:=IntegerQ[j1+j2+j3]&&Triangular[j1,j2,j3]

SixJCondition[j1_,j2_,j3_,J1_,J2_,J3_]:=
    SixJConditionTriad[j1,j2,j3] && SixJConditionTriad[j1,J2,J3]&&
    SixJConditionTriad[J1,j2,J3] && SixJConditionTriad[J1,J2,j3];
(* impose the condition explicitly in the definition of the six-j symbol 
   (the condition is already included in Mathematica,but imposing it 
   explicitly has the advantage that the six-j symbol is evaluated only 
   if it is nonzero, and so warning messages are avoided) *)
SixJ[{j1_,j2_,j3_},{J1_,J2_,J3_}]:=
    Which[SixJCondition[j1,j2,j3,J1,J2,J3],SixJSymbol[{j1,j2,j3},{J1,J2,J3}],True,0];
    
(* conditions to have a nonzero 3-j symbol: *)
ThreeJCondition[{j1_,m1_},{j2_,m2_},{j3_,m3_}]:=
    Triangular[j1,j2,j3] && RangeM[j1,m1] &&
    RangeM[j2,m2] && RangeM[j3,m3] && m1+m2+m3==0;

ThreeJ[{j1_,m1_},{j2_,m2_},{j3_,m3_}]:=
    Which[ThreeJCondition[{j1,m1},{j2,m2},{j3,m3}],
        ThreeJSymbol[{j1,m1},{j2,m2},{j3,m3}],True,0];
        
(* Nine-j symbol: *)
(* choose the cutoff of the summation in the expression of NineJSymbol: *)
cutsum=50;
SummandNineJ[j1_,j2_,j12_,j3_,j4_,j34_,j13_,j24_,j_,g_]:=
    minus^(2*g) * JNorm[g]* SixJ[{j1,j2,j12},{j34,j,g}] *
    SixJ[{j3,j4,j34},{j2,g,j24}] * SixJ[{j13,j24,j},{g,j1,j3}];


NineJSymbol[{j1_,j2_,j12_},{j3_,j4_,j34_},{j13_,j24_,j_}]:=
    Sum[SummandNineJ[j1,j2,j12,j3,j4,j34,j13,j24,j,g],{g,0,cutsum+1/2,1/2}];
(* several tests were performed.As benchmark,we used the results from the 
       website:http://www-stone.ch.cam.ac.uk/wigner.shtml *)

(* Todo:   Better base implementation *)
(* Smarter about sum range *)
(* Andrei Derevianko July,1997 *)
NineJSymbolBetter[{a_, b_, c_}, {d_, e_, f_}, {g_, h_, j_}] :=
 Module[{xLo, xUp}, 
  xLo = Max[Abs[a - j], Abs[b - f], Abs[d - h]];
  xUp = Min[a + j, b + f, d + h];
  r = Sum[(-1)^(2 x) * (2 x + 1) * SixJ[{a, b, c}, {f, j, x}] *
    SixJ[{d, e, f}, {b, x, h}] *
    SixJ[{g, h, j}, {x, a, d}], 
    {x, xLo, xUp}];
  Return[r];
];

(* Auxiliary Functions *)
Print["Loading Auxiliary Functions"];
BesselFactor1[y_,{np_,lp_},{n_,l_},lcap_]:=
    (2^lcap)/(JNorm[lcap]!!) y^(lcap/2) Exp[-y] Sqrt[(np-1)! (n-1)!];

BesselFactor2[y_,{np_,lp_},{n_,l_},lcap_]:=Sqrt[Gamma[np+lp+1/2]*Gamma[n+l+1/2]];

Summand1[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=minus^(m+mp)/(m! mp! (n-1-m)! (np-1-mp)!);

Summand2[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=
    Gamma[(l+lp+lcap+2*m+2*mp+3)/2]/(Gamma[l+m+3/2]*Gamma[lp+mp+3/2]);

Summand3[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=
    Hypergeometric1F1[(lcap-lp-l-2*mp-2*m)/2,lcap+3/2,y];

BesselFactor3[y_,{np_,lp_},{n_,l_},lcap_]:=
    Sum[(Summand1[y,{np,lp,mp},{n,l,m},lcap]*
         Summand2[y,{np,lp,mp},{n,l,m},lcap]*
         Summand3[y,{np,lp,mp},{n,l,m},lcap]),
        {m,0,n-1},{mp,0,np-1}];


BesselElement[y_,{np_,lp_},{n_,l_},lcap_]:=
    BesselFactor1[y,{np,lp},{n,l},lcap] *
    BesselFactor2[y,{np,lp},{n,l},lcap] *
    BesselFactor3[y,{np,lp},{n,l},lcap];


BesselFactor1A[y_,{np_,lp_},{n_,l_},lcap_]:=
    (2^(lcap-1))/(JNorm[lcap]!!) y^((lcap-1)/2) Exp[-y] Sqrt[(np-1)! (n-1)!];

BesselFactor2[y_,{np_,lp_},{n_,l_},lcap_]:=Sqrt[Gamma[np+lp+1/2]*Gamma[n+l+1/2]];

(**)


Summand1[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=(-1)^(m+mp)/(m! mp! (n-1-m)! (np-1-mp)!);

Summand2A[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=
    Gamma[(l+lp+lcap+2*m+2*mp+2)/2]/(Gamma[l+m+3/2]*Gamma[lp+mp+3/2]);

Summand3A[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=
    -(l+lp+lcap+2*m+2*mp+2)/2*Hypergeometric1F1[(lcap-lp-l-2*mp-2*m-1)/2,lcap+3/2,y] +
    2*m*Hypergeometric1F1[(lcap-lp-l-2*mp-2*m+1)/2,lcap+3/2,y];


BesselFactor3A[y_,{np_,lp_},{n_,l_},lcap_]:=
    Sum[(Summand1[y,{np,lp,mp},{n,l,m},lcap]*
        Summand2A[y,{np,lp,mp},{n,l,m},lcap]*
        Summand3A[y,{np,lp,mp},{n,l,m},lcap]),
        {m,0,n-1},{mp,0,np-1}];

(**)

BesselElementMinus[y_,{np_,lp_},{n_,l_},lcap_]:=
    BesselFactor1A[y,{np,lp},{n,l},lcap] *
    BesselFactor2[y,{np,lp},{n,l},lcap] *
    BesselFactor3A[y,{np,lp},{n,l},lcap];


(**)

Summand4A[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=
    -(l+lp+lcap+2*m+2*mp+2)/2*Hypergeometric1F1[(lcap-lp-l-2*mp-2*m-1)/2,lcap+3/2,y] +
    (2*l+2*m+1)*Hypergeometric1F1[(lcap-lp-l-2*mp-2*m+1)/2,lcap+3/2,y];

BesselFactor4A[y_,{np_,lp_},{n_,l_},lcap_]:=
    Sum[(Summand1[y,{np,lp,mp},{n,l,m},lcap] *
         Summand2A[y,{np,lp,mp},{n,l,m},lcap] *
         Summand4A[y,{np,lp,mp},{n,l,m},lcap]),
         {m,0,n-1},{mp,0,np-1}];

BesselElementPlus[y_,{np_,lp_},{n_,l_},lcap_]:=
    BesselFactor1A[y,{np,lp},{n,l},lcap] *
    BesselFactor2[y,{np,lp},{n,l},lcap] *
    BesselFactor4A[y,{np,lp},{n,l},lcap];

(* RELATIVISTIC Section *)
Print["Define lower-component Bessel function matrix elements (See Eq. D.16, D.26, and D.27 in Evan's thesis)"];
(* Define lower-component Bessel function matrix elements (See Eq. D.16, D.26, and D.27 in Evan's thesis)  *)
IMat[L_,m_,y_]:=
	Sqrt[Pi]/4*y^(L/2)*Exp[-y]*Gamma[1/2*(L+m+1)]/Gamma[L+3/2]*Hypergeometric1F1[1+(L-m)/2,L+3/2,y];

BesselElementLower1[y_,{np_,lp_},{n_,l_},lcap_]:=
	1/(2*Sqrt[y])*Sqrt[2*Gamma[np]*2*Gamma[n]/(Gamma[np+lp+1/2]*Gamma[n+l+1/2])] *
	 Sum[Binomial[np+lp-1/2,np-i-1]*Binomial[n+l-1/2,n-j-1]*(-1)^(i+j)/(i!*j!)*IMat[lcap,1+2*i+2*j+lp+l,y],{i,0,np-1},{j,0,n-1}];

BesselElementLower2[y_,{np_,lp_},{n_,l_},lcap_]:=
	Sqrt[2*Gamma[np]*2*Gamma[n]/(Gamma[np+lp+1/2]*Gamma[n+l+1/2])]*
	Sum[
		Binomial[np+lp-1/2,np-i-1]*Binomial[n+l-1/2,n-j-1]*(-1)^(i+j)/(i!*j!)/(2*lcap+1)
		*(lcap*IMat[lcap-1,2+2*i+2*j+lp+l,y]-(lcap+1)*IMat[lcap+1,2+2*i+2*j+lp+l,y])
	,{i,0,np-1},{j,0,n-1}];

(* END RELATIVISTIC Section *)
    
Print["Define Lnumber and Nodal"];
(* define n and l as functions of N,J,j *)
Lnumber[NPrincipal_,j_]:=Which[EvenQ[NPrincipal-(j+1/2)],(j+1/2),True,j-1/2];

Nodal[NPrincipal_,j_]:=(NPrincipal-Lnumber[NPrincipal,j])/2+1;


Print["Define parity and parity conservation support"];
(* Parity and parity conservation *)
(* define parities of an individual state and of the operators: *)
ParityState[NPrincipal_,j_]:=minus^(Lnumber[NPrincipal,j]);
(* Normal Parity: *)
ParityNormal[jcap_]:=minus^(jcap);
ParityConsNormal[NPrincipalp_,jp_,NPrincipal_,j_,jcap_]:=
    ParityState[NPrincipal,j]*ParityState[NPrincipalp,jp]*ParityNormal[jcap];

(* Function that that checks that 
    1. parity is conserved 
    2. the nodal numbers are >0 , 
    3. the momenta J,j and j' satisfy the triangular inequality 
 *)
PhysicalConditionsNormal[NPrincipalp_,jp_,NPrincipal_,j_,jcap_]:=
    ParityConsNormal[NPrincipal,j,NPrincipalp,jp,jcap]==1 && 
    Nodal[NPrincipal,j]>0 && Nodal[NPrincipalp,jp]>0 && Abs[j-jp]<=jcap<=j+jp;


(* Abnormal parity: *)
ParityAbnormal[jcap_]:=minus^(jcap+1);
ParityConsAbnormal[NPrincipalp_,jp_,NPrincipal_,j_,jcap_]:=
    ParityState[NPrincipal,j]*ParityState[NPrincipalp,jp]*ParityAbnormal[jcap];
    
(* Function that checks that
    1. parity is conserved 
    2. the nodal numbers are >0 ; 
    3. the momenta J,j and j' satisfy the triangular inequality
 *)
 PhysicalConditionsAbnormal[NPrincipalp_,jp_,NPrincipal_,j_,jcap_]:=
     ParityConsAbnormal[NPrincipal,j,NPrincipalp,jp,jcap]==1 &&
     Nodal[NPrincipal,j]>0&&Nodal[NPrincipalp,jp]>0&&Abs[j-jp]<=jcap<=j+jp;

(* Operators: normal parity *)
Print["Operators: normal parity"];
(* MJ: *)
mjelement[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    minus^(1/2+j+jcap) *
    Sqrt[(JNorm[j]*JNorm[jp]*JNorm[l]*JNorm[lp]*JNorm[jcap])/(4 Pi)] * 
    ThreeJ[{lp,0},{jcap,0},{l,0}]*SixJ[{lp,jp,1/2},{j,l,jcap}] *
    BesselElement[y,{np,lp},{n,l},jcap];
(* write it in terms of N,N',j,j': *)
MJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    mjelement[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},
        {Nodal[ncap,j],Lnumber[ncap,j],j},jcap];
(* We also impose explicitly that the function is =0 if the physical conditions are not satisfied  *)
MJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],
        MJE[y,{ncapp,jp},{ncap,j},jcap],0];
(* SigmaJ *)
MJLSigma[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=
    minus^lp *
    Sqrt[(JNorm[j]*JNorm[jp]*JNorm[l]*JNorm[lp]*JNorm[jcap]*JNorm[lcap])/(4 Pi)] *
    Sqrt[6]*
    ThreeJ[{lp,0},{lcap,0},{l,0}]*
    NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}]*
    BesselElement[y,{np,lp},{n,l},lcap];
(* write it in terms of N, N', j, j' *)
SigmaJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    MJLSigma[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},
               {Nodal[ncap,j],Lnumber[ncap,j],j},jcap,jcap];
(* We also impose explicitly that the function is =0 if the physical conditions are not satisfied *)
SigmaJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],
        SigmaJE[y,{ncapp,jp},{ncap,j},jcap],0];


(* DeltaPJ: *)
MJLDivQoverall[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=
    minus^(lcap+j+1/2)*
    QNorm[lp]*QNorm[jp]*QNorm[j]*QNorm[jcap]*QNorm[lcap] * 
    SixJ[{lp,jp,1/2},{j,l,jcap}]/Sqrt[4 Pi];

MJLDivQsummand1[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=
    -Sqrt[l+1] QNorm[l+1] *
    SixJ[{lcap,1,jcap},{l,lp,l+1}] *
    ThreeJ[{lp,0},{lcap,0},{l+1,0}] *
    BesselElementMinus[y,{np,lp},{n,l},lcap]

MJLDivQsummand2[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=
    Sqrt[l] * QNorm[l-1] * 
    SixJ[{lcap,1,jcap},{l,lp,l-1}] *
    ThreeJ[{lp,0},{lcap,0},{l-1,0}] *
    BesselElementPlus[y,{np,lp},{n,l},lcap];

MJLDivQ[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=
    MJLDivQoverall[y,{np,lp,jp},{n,l,j},jcap,lcap] *
    (MJLDivQsummand1[y,{np,lp,jp},{n,l,j},jcap,lcap]+
     MJLDivQsummand2[y,{np,lp,jp},{n,l,j},jcap,lcap]);

DeltaPrime[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    -Sqrt[jcap]/QNorm[jcap]*
    MJLDivQ[y,{np,lp,jp},{n,l,j},jcap,jcap+1] +
    Sqrt[jcap+1]/QNorm[jcap]*MJLDivQ[y,{np,lp,jp},{n,l,j},jcap,jcap-1];

(* write it in terms of N,N',j,j': *)
DeltaPrimeJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    DeltaPrime[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},
                 {Nodal[ncap,j],Lnumber[ncap,j],j},jcap];
                 
(* We also impose explicitly that the function is =0 if the physical conditions are not satisfied  *)
DeltaPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],
       DeltaPrimeJE[y,{ncapp,jp},{ncap,j},jcap], 0];


(* DeltaPPJ and DeltaTPPJ: *)
DeltaPrimePrime[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    Sqrt[jcap+1]/QNorm[jcap] * MJLDivQ[y,{np,lp,jp},{n,l,j},jcap,jcap+1] +
    Sqrt[jcap]/QNorm[jcap] * MJLDivQ[y,{np,lp,jp},{n,l,j},jcap,jcap-1];
(* write it in terms of N,N',j,j': *)
DeltaPrimePrimeJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    DeltaPrimePrime[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},
                      {Nodal[ncap,j],Lnumber[ncap,j],j},jcap];
(* We also impose explicitly that the function is =0 if the physical conditions are not satisfied  *)
DeltaPPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],
        DeltaPrimePrimeJE[y,{ncapp,jp},{ncap,j},jcap], 0];
(* define DeltaTPPJ *)
DeltaTPPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    DeltaPPJ[y,{ncapp,jp},{ncap,j},jcap] - MJ[y,{ncapp,jp},{ncap,j},jcap]/2;


(* PhiPPJ *)
PhiPPoverall[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    minus^(lp+1)*6 *QNorm[lp]*QNorm[jp]*QNorm[j]/Sqrt[4 Pi];

PhiPPsummand1[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    QNorm[l+1]*Sqrt[l+1] *
    ThreeJ[{lp,0},{jcap+1,0},{l+1,0}] *
    BesselElementMinus[y,{np,lp},{n,l},jcap+1]*
    Sum[minus^(jcap-lcap+1)*JNorm[lcap] * 
        SixJ[{jcap+1,1,lcap},{1,jcap,1}] *
        SixJ[{jcap+1,1,lcap},{l,lp,l+1}] *
        NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],
       {lcap,jcap,jcap+1}];

PhiPPsummand2[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    If[l==0,0,
        QNorm[l-1]* Sqrt[l]*
        ThreeJ[{lp,0},{jcap+1,0},{l-1,0}]*
        BesselElementPlus[y,{np,lp},{n,l},jcap+1]*
        Sum[
            minus^(jcap-lcap)*JNorm[lcap]*
            SixJ[{jcap+1,1,lcap},{1,jcap,1}]*
            SixJ[{jcap+1,1,lcap},{l,lp,l-1}]*
            NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],
           {lcap,jcap,jcap+1}]];

PhiPPsummand3[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    If[jcap==0,0,
        QNorm[l+1]*Sqrt[l+1]*
        ThreeJ[{lp,0},{jcap-1,0},{l+1,0}]*
        BesselElementMinus[y,{np,lp},{n,l},jcap-1]*
        Sum[
            minus^(jcap-lcap+1)*JNorm[lcap]*
            SixJ[{jcap-1,1,lcap},{1,jcap,1}]*
            SixJ[{jcap-1,1,lcap},{l,lp,l+1}]*
            NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],
           {lcap,jcap-1,jcap}]];

PhiPPsummand4[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    If[jcap==0,0,
        If[l==0,0,
            QNorm[l-1]* Sqrt[l]*
            ThreeJ[{lp,0},{jcap-1,0},{l-1,0}]*
            BesselElementPlus[y,{np,lp},{n,l},jcap-1]*
            Sum[
                minus^(jcap-lcap)*JNorm[lcap]*
                SixJ[{jcap-1,1,lcap},{1,jcap,1}]*
                SixJ[{jcap-1,1,lcap},{l,lp,l-1}]*
                NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],
               {lcap,jcap-1,jcap}]]];

PhiPPJX[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    PhiPPoverall[y,{np,lp,jp},{n,l,j},jcap]*
    (QNorm[jcap+1]*Sqrt[jcap+1]*
        (PhiPPsummand1[y,{np,lp,jp},{n,l,j},jcap]+
            PhiPPsummand2[y,{np,lp,jp},{n,l,j},jcap]
         ) +
        If[jcap==0,0,
            QNorm[jcap-1]*Sqrt[jcap]*
            (PhiPPsummand3[y,{np,lp,jp},{n,l,j},jcap]+PhiPPsummand4[y,{np,lp,jp},{n,l,j},jcap])
        ]
     );
(* Write PhiPP  in terms of N,N',j,j': *)
 PhiPPJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
     PhiPPJX[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},
             {Nodal[ncap,j],Lnumber[ncap,j],j},jcap];
(* We also impose explicitly that the function is = 0 if the physical conditions are not satisfied *)
PhiPPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],
        PhiPPJE[y,{ncapp,jp},{ncap,j},jcap],0];


(* PhiPJ and PhiTPJ: *)
PhiPJX[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    PhiPPoverall[y,{np,lp,jp},{n,l,j},jcap]*
    (-QNorm[jcap+1]*Sqrt[jcap]*
        (PhiPPsummand1[y,{np,lp,jp},{n,l,j},jcap]+PhiPPsummand2[y,{np,lp,jp},{n,l,j},jcap])
        +If[jcap==0,0,
            QNorm[jcap-1]*Sqrt[jcap+1]*
            (PhiPPsummand3[y,{np,lp,jp},{n,l,j},jcap]+PhiPPsummand4[y,{np,lp,jp},{n,l,j},jcap])]);
(* Write PhiP  in terms of N,N',j,j': *)
PhiPJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    PhiPJX[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},
             {Nodal[ncap,j],Lnumber[ncap,j],j},jcap];
(* We also impose explicitly that the function is = 0 if the physical conditions are not satisfied *)
PhiPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],
        PhiPJE[y,{ncapp,jp},{ncap,j},jcap],0];
(* Now evaluate PhiPT to have simple turn-around properties *)
PhiTPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    PhiPJ[y,{ncapp,jp},{ncap,j},jcap]+SigmaJ[y,{ncapp,jp},{ncap,j},jcap]/2;


(* Operators: abnormal parity *)
Print["Operators: abnormal parity"];
(* DeltaJ: *)
Deltaop[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=MJLDivQ[y,{np,lp,jp},{n,l,j},jcap,jcap];
(* write it in terms of N,N',j,j': *)
DeltaJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
     Deltaop[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},
             {Nodal[ncap,j],Lnumber[ncap,j],j},jcap];
(* We also impose explicitly that the function is=0 if the physical conditions are not satisfied *)
DeltaJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
     If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],
     DeltaJE[y,{ncapp,jp},{ncap,j},jcap],0];


(* SigmaPPJ: *)
SigmaSecond[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    (Sqrt[jcap+1] MJLSigma[y,{np,lp,jp},{n,l,j},jcap,jcap+1] +
    Sqrt[jcap] MJLSigma[y,{np,lp,jp},{n,l,j},jcap,jcap-1])/QNorm[jcap];
(* write it in terms of N,N',j,j': *)
SigmaSecondJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    SigmaSecond[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},
                  {Nodal[ncap,j],Lnumber[ncap,j],j},jcap];
(* We also impose explicitly that  the function is =0 if the physical conditions are not satisfied  *)
SigmaPPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],
        SigmaSecondJE[y,{ncapp,jp},{ncap,j},jcap],0];


(* SigmaPJ: *)
SigmaPrime[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
	-Sqrt[jcap]/QNorm[jcap]*MJLSigma[y,{np,lp,jp},{n,l,j},jcap,jcap+1]+
	Sqrt[jcap+1]/QNorm[jcap]*MJLSigma[y,{np,lp,jp},{n,l,j},jcap,jcap-1];
SigmaPrimeJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
	SigmaPrime[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},
				 {Nodal[ncap,j],Lnumber[ncap,j],j},jcap];
SigmaPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
	If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],
		SigmaPrimeJE[y,{ncapp,jp},{ncap,j},jcap],0]


(* OmegaJ and OmegaTJ: *)
MJLDivSigoverall[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    minus^(lp) *
    QNorm[lp]*QNorm[jp]*QNorm[j]*QNorm[2 j-l]*QNorm[jcap] *
    SixJ[{lp,jp,1/2},{j,2j-l,jcap}] *
    ThreeJ[{lp,0},{jcap,0},{2j-l,0}] / Sqrt[4 Pi];

MJLDivSigsummand1[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    -KroneckerDelta[j,l+1/2]*BesselElementMinus[y,{np,lp},{n,l},jcap];

MJLDivSigsummand2[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    KroneckerDelta[j,l-1/2]*BesselElementPlus[y,{np,lp},{n,l},jcap];

MJLDivSig[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    MJLDivSigoverall[y,{np,lp,jp},{n,l,j},jcap]*
    (MJLDivSigsummand1[y,{np,lp,jp},{n,l,j},jcap]+MJLDivSigsummand2[y,{np,lp,jp},{n,l,j},jcap]);

(* write Omega  in terms of N,N',j,j': *)
OmegaJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    MJLDivSig[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},
                {Nodal[ncap,j],Lnumber[ncap,j],j},jcap];

(* We also impose explicitly that  the function is =0 if the physical conditions are not satisfied *)
OmegaJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],
        OmegaJE[y,{ncapp,jp},{ncap,j},jcap],0];

(* define OmegaTJ *)
OmegaTJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
    OmegaJ[y,{ncapp,jp},{ncap,j},jcap] + SigmaPPJ[y,{ncapp,jp},{ncap,j},jcap]/2;



(* PhiJ and PhiTJ: *)
Phioverall[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    minus^(lp+1)*6 *QNorm[lp]*QNorm[jp]*QNorm[j]*JNorm[jcap]/Sqrt[4 Pi];

Phisummand1[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
    QNorm[l+1]*Sqrt[l+1] *
    ThreeJ[{lp,0},{jcap,0},{l+1,0}] *
    BesselElementMinus[y,{np,lp},{n,l},jcap] * 
    Sum[
        minus^(jcap-lcap+1) * JNorm[lcap] *
        SixJ[{jcap,1,lcap},{1,jcap,1}] * 
        SixJ[{jcap,1,lcap},{l,lp,l+1}] *
        NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],
       {lcap,Max[0,jcap-1],jcap+1}];

Phisummand2[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
	If[l==0,0,
		QNorm[l-1]* Sqrt[l] *
		ThreeJ[{lp,0},{jcap,0},{l-1,0}] *
		BesselElementPlus[y,{np,lp},{n,l},jcap] *
		Sum[
			minus^(jcap-lcap)*JNorm[lcap]*
			SixJ[{jcap,1,lcap},{1,jcap,1}]*
			SixJ[{jcap,1,lcap},{l,lp,l-1}]*
			NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],
		   {lcap,Max[0,jcap-1],jcap+1}]];

PhiJX[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
	Phioverall[y,{np,lp,jp},{n,l,j},jcap]*
	(Phisummand1[y,{np,lp,jp},{n,l,j},jcap]+Phisummand2[y,{np,lp,jp},{n,l,j},jcap]);


(* write Phi  in terms of N,N',j,j': *)
PhiJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
	PhiJX[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap];
(* We also impose explicitly that the function is = 0 if the physical conditions are not satisfied *)
PhiJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
	If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],PhiJE[y,{ncapp,jp},{ncap,j},jcap],0];
(* Now evaluate PhiT to have simple turn-around properties *)
PhiTJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=
	PhiJ[y,{ncapp,jp},{ncap,j},jcap] - SigmaPJ[y,{ncapp,jp},{ncap,j},jcap]/2;

(* RELATIVISTIC Section *)
Print["Lower Component Operators: normal parity"];
(* Lower Component Operators: normal parity *)
(* M_J^{1)} : Eq. D20 in Evan's thesis *)
m1jelement[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
	minus^(1/2+j+jcap)*Sqrt[jcap*(jcap+1)]*Sqrt[(JNorm[j]*JNorm[jp]*JNorm[l]*JNorm[lp]*JNorm[jcap])/(4*Pi)]*
	ThreeJ[{lp,0},{jcap,0},{l,0}] *
	SixJ[{lp,jp,1/2},{j,l,jcap}] *
	BesselElementLower1[y,{np,lp},{n,l},jcap];

Print["M1JE next"];

M1JE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=m1jelement[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap];

M1J[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],M1JE[y,{ncapp,jp},{ncap,j},jcap],0];

(* M_J^{2)} : Eq. D21 in Evan's thesis *)
Print["m2jelement"];
m2jelement[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
	minus^(1/2+j+jcap)*Sqrt[(JNorm[j]*JNorm[jp]*JNorm[l]*JNorm[lp]*JNorm[jcap])/(4*Pi)] *
	ThreeJ[{lp,0},{jcap,0},{l,0}] *
	SixJ[{lp,jp,1/2},{j,l,jcap}] *
	BesselElementLower2[y,{np,lp},{n,l},jcap];

M2JE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=m2jelement[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap];

M2J[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],M2JE[y,{ncapp,jp},{ncap,j},jcap],0];

Print["Lower Component Operators: abnormal parity"];
(* Lower Component Operators: abnormal parity *)
(* \Sigma_J^{\prime (0)} : Eq. D.23 in Evan's thesis *)
MJLSigma0[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=
	minus^lp Sqrt[(JNorm[j]*JNorm[jp]*JNorm[l]*JNorm[lp]*JNorm[jcap]*JNorm[lcap])/(4 Pi)] *
	Sqrt[6] *
	ThreeJ[{lp,0},{lcap,0},{l,0}] *
	NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}] *
	BesselElement[y,{np,lp},{n,l},jcap];

SigmaPrime0[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
	Sqrt[jcap]/QNorm[jcap] *MJLSigma0[y,{np,lp,jp},{n,l,j},jcap,jcap+1] +
	 Sqrt[jcap+1]/QNorm[jcap]*MJLSigma0[y,{np,lp,jp},{n,l,j},jcap,jcap-1];

SigmaPrime0JE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=SigmaPrime0[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap];

SigmaP0J[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],SigmaPrime0JE[y,{ncapp,jp},{ncap,j},jcap],0];


(* \Sigma_J^{\prime\prime (0)} : *)
Print["\\Sigma_J^{\\prime\\prime (0)}"];
SigmaPrimePrime0[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
	-Sqrt[jcap+1]/QNorm[jcap]*MJLSigma0[y,{np,lp,jp},{n,l,j},jcap,jcap+1] +
	  Sqrt[jcap]/QNorm[jcap]*MJLSigma0[y,{np,lp,jp},{n,l,j},jcap,jcap-1];

SigmaPrimePrime0JE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=SigmaPrimePrime0[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap];

SigmaPP0J[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],SigmaPrimePrime0JE[y,{ncapp,jp},{ncap,j},jcap],0];

(* \Sigma_J^{\prime (2)} : *)
Print["\\Sigma_J^{\\prime (2)}"];

MJLSigma2[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=
	minus^lp Sqrt[(JNorm[j]*JNorm[jp]*JNorm[l]*JNorm[lp]*JNorm[jcap]*JNorm[lcap])/(4 Pi)] *
	 Sqrt[6] *
	 ThreeJ[{lp,0},{lcap,0},{l,0}] *
	 NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}] *
	 BesselElementLower2[y,{np,lp},{n,l},lcap];

SigmaPrime2[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
	-Sqrt[jcap]/QNorm[jcap]*MJLSigma2[y,{np,lp,jp},{n,l,j},jcap,jcap+1] +
	 Sqrt[jcap+1]/QNorm[jcap]*MJLSigma2[y,{np,lp,jp},{n,l,j},jcap,jcap-1];

SigmaPrime2JE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=SigmaPrime2[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap];

SigmaP2J[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],SigmaPrime2JE[y,{ncapp,jp},{ncap,j},jcap],0];

(* \Sigma_J^{\prime\prime (2)} : *)
Print["\\Sigma_J^{\\prime\\prime (2)}"];

SigmaPrimePrime2[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=
	Sqrt[jcap+1]/QNorm[jcap]*MJLSigma2[y,{np,lp,jp},{n,l,j},jcap,jcap+1] +
	Sqrt[jcap]/QNorm[jcap]*MJLSigma2[y,{np,lp,jp},{n,l,j},jcap,jcap-1];

SigmaPrimePrime2JE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=SigmaPrimePrime2[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap];

SigmaPP2J[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],SigmaPrimePrime2JE[y,{ncapp,jp},{ncap,j},jcap],0];

(* Functions to sum over density matrices *)
Print["Functions to sum over density matrices"];
ISOT[TT_,MT_,T0_]:=KroneckerDelta[T0,0]Sqrt[2]/QNorm[TT]+KroneckerDelta[T0,1]  If[MT==0,0,MT Sqrt[6 /((2TT+1)(TT+1)TT)]];
FM[ii_,jj_,y_]:=FM[ii,jj,y]=
	Simplify[
		ISOT[TMTS[[Ncode,ii,1]],TMTS[[Ncode,ii,2]],T0[ii,jj]]*
		Sum[
			DM[ii,jj][[kk,3]] MJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]]
		,{kk,1,LDM[ii,jj]}]];

FPhiPP[ii_,jj_,y_]:=FPhiPP[ii,jj,y]=
	Simplify[ISOT[TMTS[[Ncode,ii,1]],TMTS[[Ncode,ii,2]],T0[ii,jj]]*
	   Sum[
	   	DM[ii,jj][[kk,3]] PhiPPJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]]
	   ,{kk,1,LDM[ii,jj]}]];

(* Now for spin dependent cases *)
FSigmaPP[ii_,jj_,y_]:=FSigmaPP[ii,jj,y]=
	Simplify[
		ISOT[TMTS[[Ncode,ii,1]],TMTS[[Ncode,ii,2]],T0[ii,jj]]*
		Sum[
			DM[ii,jj][[kk,3]] SigmaPPJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]]
		,{kk,1,LDM[ii,jj]}]];
FDelta[ii_,jj_,y_]:=FDelta[ii,jj,y]=
	Simplify[
		ISOT[TMTS[[Ncode,ii,1]],TMTS[[Ncode,ii,2]],T0[ii,jj]]*
		Sum[
			DM[ii,jj][[kk,3]] DeltaJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]]
		,{kk,1,LDM[ii,jj]}]];
FSigmaP[ii_,jj_,y_]:=FSigmaP[ii,jj,y]=
	Simplify[
		ISOT[TMTS[[Ncode,ii,1]],TMTS[[Ncode,ii,2]],T0[ii,jj]]*
		Sum[
			DM[ii,jj][[kk,3]] SigmaPJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]]
		,{kk,1,LDM[ii,jj]}]];
FPhiTP[ii_,jj_,y_]:=FPhiTP[ii,jj,y]=
	Simplify[
		ISOT[
			TMTS[[Ncode,ii,1]],TMTS[[Ncode,ii,2]],T0[ii,jj]] *
			Sum[
				DM[ii,jj][[kk,3]] PhiTPJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]]
			,{kk,1,LDM[ii,jj]}]];

(* RELATIVISTIC Section *)
(* Lower component operators *)
FM1[ii_,jj_,y_]:=FM1[ii,jj,y]=Simplify[ISOT[TMTS[[Ncode,ii,1]],TMTS[[Ncode,ii,2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] M1J[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]];
FM2[ii_,jj_,y_]:=FM2[ii,jj,y]=Simplify[ISOT[TMTS[[Ncode,ii,1]],TMTS[[Ncode,ii,2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] M2J[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]];
FSigmaP0[ii_,jj_,y_]:=FSigmaP0[ii,jj,y]=Simplify[ISOT[TMTS[[Ncode,ii,1]],TMTS[[Ncode,ii,2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] SigmaP0J[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]];
FSigmaPP0[ii_,jj_,y_]:=FSigmaPP0[ii,jj,y]=Simplify[ISOT[TMTS[[Ncode,ii,1]],TMTS[[Ncode,ii,2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] SigmaPP0J[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]];
FSigmaP2[ii_,jj_,y_]:=FSigmaP2[ii,jj,y]=Simplify[ISOT[TMTS[[Ncode,ii,1]],TMTS[[Ncode,ii,2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] SigmaP2J[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]];
FSigmaPP2[ii_,jj_,y_]:=FSigmaPP2[ii,jj,y]=Simplify[ISOT[TMTS[[Ncode,ii,1]],TMTS[[Ncode,ii,2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] SigmaPP2J[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]];


Print["Defining Lower Component operators Eq. B4 in long paper"];
(* Lower Component operators Eq. B4 in long paper *)
RM1[i_,j_,bs_]:=bs[3,i]Conjugate[bs[3,j]] +bs[7,i] Conjugate[bs[7,j]];
RM2[i_,j_,bs_]:=bs[2,i]Conjugate[bs[2,j]] +bs[7,i] Conjugate[bs[7,j]];
RSigmaP0[i_,j_,bs_]:=(bs[12,i]-bs[15,i])Conjugate[bs[12,j]-bs[15,j]] +bs[14,i] Conjugate[bs[14,j]];
RSigmaP2[i_,j_,bs_]:=(bs[13,i]-bs[14,i])Conjugate[bs[13,j]-bs[14,j]] +bs[15,i] Conjugate[bs[15,j]];
RSigmaPP0[i_,j_,bs_]:=bs[8,i]Conjugate[bs[8,j]]+bs[13,i] Conjugate[bs[13,j]];
RSigmaPP2[i_,j_,bs_]:=(bs[13,i]-bs[14,i])Conjugate[bs[13,j]-bs[14,j]] +bs[16,i] Conjugate[bs[16,j]];
RSigmaP0SigmaP2[i_,j_,bs_]:=Re[bs[14,i]Conjugate[bs[13,j]-bs[14,j]] +(bs[12,i]-bs[15,i])Conjugate[bs[15,j]]];
RSigmaPP0SigmaPP2[i_,j_,bs_]:=Re[bs[8,i]Conjugate[bs[16,j]]+bs[13,i]Conjugate[bs[13,j]-bs[14,j]]];

(* Upper/lower interference: *)
Print["Defining RMM2 .. RSigmaPP2SigmaPP"];
RMM2[i_,j_,cs_,bs_]:=Im[cs[1,i]Conjugate[bs[2,j]] -cs[11,i] Conjugate[bs[7,j]]];
RPhiPPM2[i_,j_,cs_,bs_]:=Im[(cs[12,i]-cs[15,i])Conjugate[bs[7,j]] +cs[3,i] Conjugate[bs[2,j]]];
RSigmaP0SigmaP[i_,j_,cs_,bs_]:=Im[bs[14,i]Conjugate[cs[4,j]] -(bs[12,i]-bs[15,i]) Conjugate[cs[9,j]]];
RSigmaP2SigmaP[i_,j_,cs_,bs_]:=Im[-(bs[13,i]-bs[14,i]) Conjugate[cs[4,j]]+bs[15,i]Conjugate[cs[9,j]]];
RDeltaSigmaP0[i_,j_,cs_,bs_]:=Im[cs[5,i]Conjugate[bs[14,j]]-cs[8,i]Conjugate[bs[12,j]-bs[15,j]]];
RDeltaSigmaP2[i_,j_,cs_,bs_]:=Im[cs[8,i]Conjugate[bs[15,j]]-cs[5,i]Conjugate[bs[13,j]-bs[14,j]]];
RPhiTPM1[i_,j_,cs_,bs_]:=Im[cs[12,i]Conjugate[bs[7,j]]+cs[13,i]Conjugate[bs[3,j]]];
RSigmaPP0SigmaPP[i_,j_,cs_,bs_]:=Im[bs[13,i] Conjugate[cs[4,j]-cs[6,j]]-bs[8,i]Conjugate[cs[10,j]]];
RSigmaPP2SigmaPP[i_,j_,cs_,bs_]:=Im[(bs[13,i]-bs[14,i])* Conjugate[cs[4,j]-cs[6,j]]-bs[16,i]Conjugate[cs[10,j]]];

(* End new RELATIVISTIC Section *)
Print["Defining RM .. RDeltaSigmaP"];

(*Now define the effective theory input, with the isospin indices {isoscalar,isovector} specified by {i,j}*) 
RM[i_,j_,cs_]:= cs[1,i]Conjugate[cs[1,j]] + cs[11,i] Conjugate[cs[11,j]];
RPhiPP[i_,j_,cs_]:=  cs[3,i] Conjugate[cs[3,j]] + (cs[12,i]- cs[15,i]) Conjugate[cs[12,j]- cs[15,j]];
RPhiPPM[i_,j_,cs_]:=Re[ cs[3,i]Conjugate[cs[1,j]] - (cs[12,i]-cs[15,i]) Conjugate[cs[11,j]]];
RPhiTP[i_,j_,cs_]:= cs[12,i] Conjugate[cs[12,j]] + cs[13,i] Conjugate[cs[13,j]];
RSigmaPP[i_,j_,cs_]:=cs[10,i] Conjugate[cs[10,j]] +(cs[4,i]-cs[6,i])Conjugate[ cs[4,j]-cs[6,j]] ;
RSigmaP[i_,j_,cs_]:=cs[4,i] Conjugate[cs[4,j]]+  cs[9,i] Conjugate[cs[9,j]] ;
RDelta[i_,j_,cs_]:= cs[5,i] Conjugate[cs[5,j]] +cs[8,i] Conjugate[cs[8,j]];
RDeltaSigmaP[i_,j_,cs_]:= Re[cs[5,i] Conjugate[cs[4,j]]+ cs[8,i]Conjugate[cs[9,j]]];

Print["Got to WM .. WDelta definitions"];

(*We now calculate the total response functions, summed over isotopes,
  weighted by the abundance,again with very small kinematic differences 
  between isotopes due to masses ignored, to be later approximated by an 
  effective average nuclear mass
 *)
WM[qm_,y_,cs_]:=Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FM[ii,jj,y] FM[ii,jjp,y]RM[T0[ii,jj],T0[ii,jjp],cs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WPhiPP[qm_,y_,cs_]:=qm^2 Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FPhiPP[ii,jj,y] FPhiPP[ii,jjp,y]RPhiPP[T0[ii,jj],T0[ii,jjp],cs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WPhiPPM[qm_,y_,cs_]:=-2 qm Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FPhiPP[ii,jj,y] FM[ii,jjp,y]RPhiPPM[T0[ii,jj],T0[ii,jjp],cs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WPhiTP[qm_,y_,cs_]:=qm^2 Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FPhiTP[ii,jj,y] FPhiTP[ii,jjp,y]RPhiTP[T0[ii,jj],T0[ii,jjp],cs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WSigmaPP[qm_,y_,cs_]:=Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FSigmaPP[ii,jj,y] FSigmaPP[ii,jjp,y]RSigmaPP[T0[ii,jj],T0[ii,jjp],cs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WDelta[qm_,y_,cs_]:=qm^2 Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FDelta[ii,jj,y] FDelta[ii,jjp,y]RDelta[T0[ii,jj],T0[ii,jjp],cs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WSigmaP[qm_,y_,cs_]:=Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FSigmaP[ii,jj,y] FSigmaP[ii,jjp,y]RSigmaP[T0[ii,jj],T0[ii,jjp],cs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WDeltaSigmaP[qm_,y_,cs_]:=-2 qm Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FDelta[ii,jj,y] FSigmaP[ii,jjp,y]RDeltaSigmaP[T0[ii,jj],T0[ii,jjp],cs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];

(* RELATIVISTIC Section *)
Print["Lower component operators: Eq. B3, B5 in long paper"];
(* Lower component operators: Eq. B3, B5 in long paper *)
WM1[qm_,y_,bs_]:=Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FM1[ii,jj,y] FM1[ii,jjp,y]RM1[T0[ii,jj],T0[ii,jjp],bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WM2[qm_,y_,bs_]:=Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FM2[ii,jj,y] FM2[ii,jjp,y]RM2[T0[ii,jj],T0[ii,jjp],bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WSigmaP0[qm_,y_,bs_]:=Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FSigmaP0[ii,jj,y] FSigmaP0[ii,jjp,y]RSigmaP0[T0[ii,jj],T0[ii,jjp],bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WSigmaP2[qm_,y_,bs_]:=Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FSigmaP2[ii,jj,y] FSigmaP2[ii,jjp,y]RSigmaP2[T0[ii,jj],T0[ii,jjp],bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WSigmaPP0[qm_,y_,bs_]:=Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FSigmaPP0[ii,jj,y] FSigmaPP0[ii,jjp,y]RSigmaPP0[T0[ii,jj],T0[ii,jjp],bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WSigmaPP2[qm_,y_,bs_]:=Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FSigmaPP2[ii,jj,y] FSigmaPP2[ii,jjp,y]RSigmaPP2[T0[ii,jj],T0[ii,jjp],bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WMM2[qm_,y_,cs_,bs_]:=2 Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FM[ii,jj,y] FM2[ii,jjp,y]RMM2[T0[ii,jj],T0[ii,jjp],cs,bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WPhiPPM2[qm_,y_,cs_,bs_]:=-2qm Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FPhiPP[ii,jj,y] FM2[ii,jjp,y]RPhiPPM2[T0[ii,jj],T0[ii,jjp],cs,bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WSigmaP0SigmaP[qm_,y_,cs_,bs_]:=-2Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FSigmaP0[ii,jj,y] FSigmaP[ii,jjp,y]RSigmaP0SigmaP[T0[ii,jj],T0[ii,jjp],cs,bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WSigmaP2SigmaP[qm_,y_,cs_,bs_]:=-2Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FSigmaP2[ii,jj,y] FSigmaP[ii,jjp,y]RSigmaP2SigmaP[T0[ii,jj],T0[ii,jjp],cs,bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WSigmaP0SigmaP2[qm_,y_,bs_]:=-2Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FSigmaP0[ii,jj,y] FSigmaP2[ii,jjp,y]RSigmaP0SigmaP2[T0[ii,jj],T0[ii,jjp],bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WDeltaSigmaP0[qm_,y_,cs_,bs_]:=-2qm*Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FDelta[ii,jj,y] FSigmaP0[ii,jjp,y]RDeltaSigmaP0[T0[ii,jj],T0[ii,jjp],cs,bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WDeltaSigmaP2[qm_,y_,cs_,bs_]:=-2qm*Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FDelta[ii,jj,y] FSigmaP2[ii,jjp,y]RDeltaSigmaP2[T0[ii,jj],T0[ii,jjp],cs,bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WPhiTPM1[qm_,y_,cs_,bs_]:=-2qm*Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FPhiTP[ii,jj,y] FM1[ii,jjp,y]RPhiTPM1[T0[ii,jj],T0[ii,jjp],cs,bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WSigmaPP0SigmaPP[qm_,y_,cs_,bs_]:=2Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FSigmaPP0[ii,jj,y] FSigmaPP[ii,jjp,y]RSigmaPP0SigmaPP[T0[ii,jj],T0[ii,jjp],cs,bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WSigmaPP2SigmaPP[qm_,y_,cs_,bs_]:=2Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FSigmaPP2[ii,jj,y] FSigmaPP[ii,jjp,y]RSigmaPP2SigmaPP[T0[ii,jj],T0[ii,jjp],cs,bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];
WSigmaPP0SigmaPP2[qm_,y_,bs_]:=2Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FSigmaPP0[ii,jj,y] FSigmaPP2[ii,jjp,y]RSigmaPP0SigmaPP2[T0[ii,jj],T0[ii,jjp],bs],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}];

(* End RELATIVISTIC Section *)

Heff[qm_,y_,cs_,bs_]:=
	If[mL>0&&iLower==1,
		Simplify[Expand[WM[qm,y,cs]+WPhiPP[qm,y,cs]+WPhiPPM[qm,y,cs]+WPhiTP[qm,y,cs]
		   +WSigmaPP[qm,y,cs]+WDelta[qm,y,cs]+WSigmaP[qm,y,cs]+WDeltaSigmaP[qm,y,cs]
       	   +fgAvg^2*(WM1[qm,y,bs]+WM2[qm,y,bs]+WSigmaP0[qm,y,bs]+WSigmaP2[qm,y,bs]
		      +WSigmaPP0[qm,y,bs]+WSigmaPP2[qm,y,bs]+WSigmaP0SigmaP2[qm,y,bs]
			  +WSigmaPP0SigmaPP2[qm,y,bs])
		   +fgAvg*(WMM2[qm,y,cs,bs]+WPhiPPM2[qm,y,cs,bs]+WPhiTPM1[qm,y,cs,bs]
		      +WSigmaP0SigmaP[qm,y,cs,bs]+WSigmaP2SigmaP[qm,y,cs,bs]
			  +WDeltaSigmaP0[qm,y,cs,bs]+WDeltaSigmaP2[qm,y,cs,bs]
			  +WSigmaPP0SigmaPP[qm,y,cs,bs]+WSigmaPP2SigmaPP[qm,y,cs,bs])]]
,
		Simplify[Expand[WM[qm,y,cs]+WPhiPP[qm,y,cs]+WPhiPPM[qm,y,cs]
			+WPhiTP[qm,y,cs]+WSigmaPP[qm,y,cs]+WDelta[qm,y,cs]+WSigmaP[qm,y,cs]
			+WDeltaSigmaP[qm,y,cs]]]
	];

(* 
 * Response function reporting support 
 * The original code only defines these functions when 
 * either reporting analytic results for response functions or
 * when plotting them.   Other than a little runtime, there is
 * no reason not to always define the functions.
 *
 * One issue is that the functions are specialized for the value of Ncode
 * enabling them to be simplified.
 *)
defineResponseFunctions[Ncode_, isod_] := Module[{ResMV, ResSigmaPPV, ResSigmaPV, ResSDV, ResDeltaV, ResPhiPPV, ResPhiTPV },
	ResMV=Simplify[Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FM[ii,jj,y] FM[ii,jjp,y]isod[[T0[ii,jj]+1]] isod[[T0[ii,jjp]+1]],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}]];
	Clear[ResM];
	Global`ResM[yy_]:=ResMV/.y->yy;

	ResSigmaPPV=Simplify[Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FSigmaPP[ii,jj,y] FSigmaPP[ii,jjp,y]isod[[T0[ii,jj]+1]] isod[[T0[ii,jjp]+1]],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}]];
	Clear[ResSigmaPP];
	Global`ResSigmaPP[yy_]:=ResSigmaPPV /. y->yy;

	ResSigmaPV=Simplify[Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FSigmaP[ii,jj,y] FSigmaP[ii,jjp,y]isod[[T0[ii,jj]+1]] isod[[T0[ii,jjp]+1]],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}]];
	Clear[ResSigmaP];
	Global`ResSigmaP[yy_]:=ResSigmaPV /. y->yy;

	ResSDV=Simplify[ResSigmaPP[y]+ResSigmaP[y]];
	Clear[ResSD];
	Global`ResSD[yy_]:=ResSDV /.y->yy;

	ResDeltaV=Simplify[Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FDelta[ii,jj,y] FDelta[ii,jjp,y]isod[[T0[ii,jj]+1]] isod[[T0[ii,jjp]+1]],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}]];
	Clear[ResDelta];
	Global`ResDelta[yy_]:=ResDeltaV /. y->yy;

	ResPhiPPV=Simplify[Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FPhiPP[ii,jj,y] FPhiPP[ii,jjp,y]isod[[T0[ii,jj]+1]] isod[[T0[ii,jjp]+1]],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}]];
	Clear[ResPhiPP];
	Global`ResPhiPP[yy_]:= ResPhiPPV /. y->yy;

	ResPhiTPV=Simplify[Sum[AbNorm[ii](4 Pi)/(2JI[ii]+1)Sum[If[J0[ii,jj]==J0[ii,jjp],FPhiTP[ii,jj,y] FPhiTP[ii,jjp,y]isod[[T0[ii,jj]+1]] isod[[T0[ii,jjp]+1]],0],{jjp,1,NumDMs[ii]},{jj,1,NumDMs[ii]}],{ii,1,NumAs[[Ncode]]}]];
	Clear[ResPhiTP];
	Global`ResPhiTP[yy_]:=ResPhiTPV /. y->yy;
];

(* Report single response, ask gives user the option to report *)
reportResponse[ask_, name_, f_, y_] := Module[{ichoice},
	If[ask, 
		ichoice = InputString[" " <> name <> "\n \n Press Enter or 1 then Enter to calculate, 0 otherwise", ""];
		ichoice = StringTrim[ichoice];
	,
		ichoice = "1";
	];

	If[ichoice == "1" || ichoice == "",
		Print[name <> ":"];
		Print[f[y]];
	];
];

(* Report response functions.   ask gives user the option for each one, False -> all *)
reportResponses[ask_, y_]:= Module[{},
	reportResponse[ask, "Vector charge response (M)", ResM, y];
	reportResponse[ask, "Axial longitudinal spin response (Sigma'')", ResSigmaPP, y];
	reportResponse[ask, "Axial transverse spin response (Sigma')", ResSigmaP, y];
	reportResponse[ask, "Standard spin-dependent response \n (Sum of longitudinal and tranverse spin responses)", ResSD, y];
	reportResponse[ask, "Vector transverse magnetic response (Delta)", ResDelta, y];
	reportResponse[ask, "Vector longitudinal response (Phi'')", ResPhiPP, y];
	reportResponse[ask, "Vector transverse electric response (PhiT')", ResPhiTP, y];
];

(* Plot a single response function.   ask gives user the option *)
plotResponse[ask_, res_, name_, f_, y_] := Module[{plt, ichoice},
	ichoice = "1";
	If[ask,
		ichoice = InputString[" " <> name <> "\n \n Press Enter or 1 then Enter to plot, 0 otherwise", ""];
		ichoice = StringTrim[ichoice];
	];
	If[ichoice == "1" || ichoice == "",
		Print["Plot of " <> name <> " as function of y:"];
		plt = Plot[f[y], {y, 0, 2}, PlotRange->All, AxesOrigin->{0,0}];
		If[SameQ[scriptBaseFile, None],
			Print[plt];
		,
			fn = scriptBaseFile <> "_" <> res <> ".pdf";
			Print["Exporting plot ", fn];
			Export[fn, plt];
		];
	];
];

(* These are also called from Yaml->Mathematica generated script *)
(* Note that scriptBaseFile is set in the scripts, but not in the interactive flow *)
plotResponseVcrm[ask_] := Module[{},
	Print["In plotResponseVcrm with ask=", ask];
	plotResponse[ask, "vcrm", "Vector charge response (M)", ResM, y];
];
plotResponseAlsr[ask_] := plotResponse[ask, "alsr", "Axial longitudinal spin response (Sigma'')", ResSigmaPP, y];
plotResponseAtsr[ask_] := plotResponse[ask, "Atsr", "Axial transverse spin response (Sigma')", ResSigmaP, y];
plotResponseSsd[ask_] := plotResponse[ask, "ssd",
	"Standard spin-dependent response \n (Sum of longitudinal and tranverse spin responses)", ResSD, y];
plotResponseVtmr[ask_] := plotResponse[ask, "vtmr", "Vector transverse magnetic response (Delta)", ResDelta, y];
plotResponseVlr[ask_] := plotResponse[ask, "vlr", "Vector longitudinal response (Phi'')", ResPhiPP, y];
plotResponseVter[ask_] := plotResponse[ask, "vter", "Vector transverse electric response (PhiT')", ResPhiTP, y];

(* Plot response functions.   ask gives user the interactive option for each one, False -> all *)
plotResponses[ask_]:= Module[{},
	plotResponseVcrm[ask];
	plotResponseAlsr[ask];
	plotResponseAtsr[ask];
	plotResponseSsd[ask];
	plotResponseVtmr[ask];
	plotResponseVlr[ask];
	PlotResponseVter[ask];
];

plotResponseD[data_] := Module[{plist, si, s},
	Print["In plotResponseD"];
	If[! KeyExistsQ[data, "plots"], Return[]]; (* nothing to do *)
	plist = data["plots"];
	Print["plotResponseD with ", plist];
	If[plist == "none", Return[]];
	If[plist == "all",
		plotResponses[False]; (* False says don't ask, just do *)
		Return[];
	];
	For[si = 1, si <= Length[plist], si++,
		s = plist[[si]];
		Print["Specific Plot ", s];
		If[s == "vcrm", plotResponseVcrm[False]];
		If[s == "alsr", plotResponseAlsr[False]];
		If[s == "atsr", plotResponseAtsr[False]];
		If[s == "ssd",  plotResponseSsd[False]];
		If[s == "vtmr", plotResponseVtmr[False]];
		If[s == "vlr",  plotResponseVlr[False]];
		If[s == "vter", plotResponseVter[False]];
	];
];

Print["Loaded mu2elib.m"];

