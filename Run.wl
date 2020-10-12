(* ::Package:: *)

(* ::Subsubsection::Closed:: *)
(*Importing Modules*)


(* ::Text:: *)
(*Imports ModuleNaming module and FixedPoint library*)


gridSize = 120;
zNormalization = 0;


arguments = $CommandLine;
If[MemberQ[$CommandLine,"-script"], debug = False, debug = True];
SetAttributes[debug,Protected];

If[debug,
	SetDirectory["/home/andrzej/FunctionalRenormalization"];
	Print["Debug mode is active!"],
	
	notebookPath = arguments[[3]];
	slashes = StringPosition[notebookPath,"/"];
	If[Length[slashes]!=0,
		projectPath = StringTake[notebookPath, slashes[[-1,1]]];
		SetDirectory[projectPath];
	];];

Get["ModuleNaming.wl"];
Get["FixedPoint.wl"];
Get["Plots.wl"];


(* ::Subsubsection::Closed:: *)
(*Initialization*)


(* ::Text:: *)
(*Serializes the command line arguments*)


If[debug,
	(* Debug run configuration *)
	task = "ON"; (* ON/ Z4 / tricritical/ regulator*)
	configuration = 1;
	ON = True;
	regulatorName = "smooth",
	
	(* Release run configuration *)
	task = Null;
	anisotropySymmetry = 4;
	regulatorName = Null;
	configuration = Null;
	Block[{i},
	For[i = 1, i<=Length[arguments], i++,
		If[arguments[[i]] == "-task", task = arguments[[i+1]]];
		If[arguments[[i]] == "-regulator", regulatorName = arguments[[i+1]]];
		If[arguments[[i]] == "-conf", configuration = ToExpression[arguments[[i+1]]]];
		];
	];
];
ON = MemberQ[ONTASKS, task];

(* Check if all required options are specified *)
If[task == Null || regulatorName == Null || (anisotropySymmetry == "ON" && nn == Null),
	Print["Specify all required options: task, regulator"];
	Quit[];
	];

(* Select regulator function *)
SelectRegulator[regulatorName];
cacheDirectory = GetCacheDirectory[task, 0, regulatorName];

(* Select run configuration file *)
runconf = Null;
If[task=="ON", runconf = ONRUNCONF];
If[task=="Z4", runconf = Z4RUNCONF];
If[task=="tricritical", runconf = TRICRITICALRUNCONF];
If[task=="lowd", runconf = LOWDRUNCONF];
If[task=="regulator", runconf = ""];

If[runconf == Null, 
	Print["Invalid task specified"]; 
	Quit[]];

If[task == "tricritical",
	searchTricritical = True;
	guessfile = TRICRITICALGUESSFILE,
	If[task == "lowd",
		searchTricritical = False;
		guessfile = LOWDGUESSFILE,
	
		searchTricritical = False;
		guessfile = INITIALGUESSFILE
		];
	];


(* ::Subsubsection::Closed:: *)
(*Importing equations*)


(* ::Text:: *)
(*Imports cached flow equations*)


isotropicEquationsFile = GetEquationsFileName[task, 0, 0];
If[FileExistsQ[isotropicEquationsFile], 
	integrandsListIsotropic = Import[isotropicEquationsFile],
	Print["No isotropic equations found"]; Quit[]];
If[ON,
integrandsListAnisotropic = {},
anisotropicEquationsFile = GetEquationsFileName[task, 0, 1];
If[FileExistsQ[anisotropicEquationsFile], 
	integrandsListAnisotropic = Import[anisotropicEquationsFile],
	Print["No anisotropic equations found"]; Quit[]];
];


(* ::Subsubsection::Closed:: *)
(*O(N) dimension dependence*)


If[ task == "ON" || task == "tricritical",
{guess, dim, \[Rho]Max} = Import[guessfile];
guess = AdjustNumberOfPoints[guess, gridSize,zNormalization];

configurations = Import[runconf];

If[!Element[configuration,Integers] || configuration<1 || configuration > Length[configurations],
	Print["Choose proper run cofiguration id"];
	Quit[]];
constants = configurations[[configuration]];

If[MemberQ[Keys[constants], "N"], 
	If[constants["N"] > 4,
		eqlistIntermediate = IntegrateEquations[integrandsListIsotropic /. n -> (constants["N"]+2)/2.];
		{guess, \[Rho]Max} = \[Rho]MaxScan[eqlistIntermediate, guess, dim, \[Rho]Max, 2., 0.1*\[Rho]Max, 20];
	];
	integrandsListIsotropic = (integrandsListIsotropic /. n -> constants["N"]);
	integrandsListAnisotropic = (integrandsListAnisotropic /. n -> constants["N"])];

eqlist = IntegrateEquations[integrandsListIsotropic];

scan = \[Rho]MaxScan[eqlist, guess, dim, \[Rho]Max, 2., 0.1*\[Rho]Max, 20];
If[Length[scan]!=2, Print["Initial fixed point search failed"]; Quit[];];
{newGuess, \[Rho]Max} = scan;

fixedPointDict = FixedPointDimScan[eqlist, newGuess, constants, dim, \[Rho]Max];

exponents = FindExponents[integrandsListIsotropic, integrandsListAnisotropic, fixedPointDict, {"dim","\[Rho]Max"}, constants];

If[ON,
	nestedFPD = NestAssociation[fixedPointDict];
	fixedPointDict = UnnestAssociation[Association[constants[["N"]]->nestedFPD]];];

If[ON,
	exponents = UnnestAssociation[TransposeAssociation[Association[constants[["N"]]->exponents]]],
	exponents = UnnestAssociation[exponents]];
Export[cacheDirectory<>FIXEDPOINTSFILES[configuration],{fixedPointDict,exponents}];
	
If[FileExistsQ[cacheDirectory<>FIXEDPOINTSFILE],
	{allFPs, allExponents} = Import[cacheDirectory<>FIXEDPOINTSFILE];
	allFPs = Join[allFPs, fixedPointDict];
	allExponents = Join[allExponents, exponents],
	
	allFPs = fixedPointDict;
	allExponents = exponents];
	
Export[cacheDirectory<>FIXEDPOINTSFILE,{allFPs,allExponents}];
];


(* ::Subsubsection::Closed:: *)
(*Low-d task computation*)


If[task=="lowd",
{guess, dim, \[Rho]Max} = Import[guessfile];
guess = AdjustNumberOfPoints[guess, gridSize,zNormalization];

{keyLabels, configurations} = Import[runconf];

If[!Element[configuration,Integers] || configuration<1 || configuration > Length[configurations],
	Print["Choose proper run cofiguration id"];
	Quit[]
	];
{keys, constants} = configurations[[configuration]];

If[MemberQ[Keys[constants], "N"], 
	integrandsListIsotropic = (integrandsListIsotropic /. n -> constants["N"]);
	integrandsListAnisotropic = (integrandsListAnisotropic /. n -> constants["N"])];

eqlist = IntegrateEquations[integrandsListIsotropic];

{newGuess, new\[Rho]Max} = \[Rho]MaxScan[eqlist /. n->2, guess, dim, \[Rho]Max, 2., 0.025*\[Rho]Max, 20];
constants[["\[Rho]Max"]] = new\[Rho]Max;

fixedPointDict = FindFixedPointDict[eqlist, newGuess, keys, keyLabels, constants];
	
AppendTo[keyLabels,"\[Rho]Max"];
constants = Drop[constants, Position[Keys[constants],"\[Rho]Max"][[1]]];
exponents = FindExponents[integrandsListIsotropic, integrandsListAnisotropic, fixedPointDict, keyLabels, constants];

If[ON && MemberQ[Keys[constants],"N"],
	nestedFPD = NestAssociation[fixedPointDict];
	fixedPointDict = UnnestAssociation[Association[constants[["N"]]->nestedFPD]];
	
	exponents = UnnestAssociation[TransposeAssociation[Association[constants[["N"]]->exponents]]],
	exponents = UnnestAssociation[exponents]];

Export[cacheDirectory<>FIXEDPOINTSFILES[configuration],{fixedPointDict,exponents}];
	
If[FileExistsQ[cacheDirectory<>FIXEDPOINTSFILE],
	{allFPs, allExponents} = Import[cacheDirectory<>FIXEDPOINTSFILE];
	allFPs = Join[allFPs, fixedPointDict];
	allExponents = Join[allExponents, exponents],
	
	allFPs = fixedPointDict;
	allExponents = exponents];
	
Export[cacheDirectory<>FIXEDPOINTSFILE,{allFPs,allExponents}];

If[ON,
	plotExp = NestAssociation[FindPMSValues[allExponents]];
	For[i=1, i<=Length[plotExp], i++,
		PlotDataSeries[plotExp[[i]], {"d", Keys[plotExp][[i]]}, "N", cacheDirectory];
	];]
];


(* ::Subsubsection::Closed:: *)
(*Cardy line plotting*)


If[ task == "ON"||task == "tricritical",
onDirectory = GetCacheDirectory["ON", 0, regulatorName];
triDirectory = GetCacheDirectory["tricritical", 0, regulatorName];

points[cache_] := Block[{ns,dims,nested, fps,exp},
	{fps,exp} = Import[cache];
	nested = NestAssociation[fps];
	ns = Keys[nested];
	dims = Keys[Values[nested]][[All,-1]];
	Return[Map[{#[[2]],#[[1]]}&,Sort[Transpose[{ns,dims}]]]];
];

onPoints = points[onDirectory<>FIXEDPOINTSFILE];
triPoints = points[triDirectory<>FIXEDPOINTSFILE];

cacheDirectory = GetCacheDirectory[task, 0, regulatorName];

min = Min[triPoints[[All,1]],onPoints[[All,1]]];
max = Max[triPoints[[All,1]],onPoints[[All,1]]];
plot = Show[ListLinePlot[{onPoints,triPoints},PlotLegends->PointLegend[{"Critical","Tricritical"},LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 20]],PlotMarkers->{Automatic, 7}],
	Plot[2+ (d-2)*Pi^2/4,{d,min,max},PlotLegends->LineLegend[{"Cardy"},LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 20]],PlotStyle->Red],
AxesLabel->{"d","N"}, LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 16],
ImageSize->Full];
Export["cardy.png",plot];
];


(* ::Subsubsection::Closed:: *)
(*Critical exponents dimension dependence plotting*)


If[ task == "ON" || task == "tricritical",
{fps,allExponents} = Import[cacheDirectory<>FIXEDPOINTSFILE];
allExponents = NestAssociation[allExponents];
PlotDimDependence[exp_] := Block[{data,min,max,p},
data = KeySort[allExponents[exp]];
If[exp==e1, data["\[Infinity]"] = Association[2->0,3->1]];
min = Min[Values[data]];
If[exp==\[Nu], min=Max[0,min]];
max = Max[Values[data]];
If[exp==\[Nu], data["\[Infinity]"] = Association[Table[(3.-i/100)->1./(1-i/100),{i,0,99}]]];
p = ListLinePlot[Values[data],PlotLegends->PointLegend[Keys[data], LegendLabel->"N", LegendFunction->"Panel",
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 20]], LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 16],
	AxesLabel->{"d",exp}, ImageSize->Full, PlotRange->{{1,3},{Min[min 1.2, 0.8 min, 0], Max[max 1.2, 0.8 max, 0]}}];
Export[cacheDirectory<>ToString[exp]<>".png",p];];
Map[PlotDimDependence, Keys[allExponents]];
];


(* ::Subsubsection::Closed:: *)
(*Regulator dependence analysis*)


If[task=="regulator",
{fps,exps} = Import[cacheDirectory<>FIXEDPOINTSFILE];
nestedfps = NestAssociation[fps];

nestedfps = KeySort[Association[Table[N[Keys[nestedfps][[i]]] ->nestedfps[[i]],{i,1,Length[nestedfps]}]]];
	
ns = {1.,2.,3.,4.,6.,8.};
alphasLow = Table[2-i 0.25,{i,0,6}];
alphasHigh = Table[2+i 0.25,{i,1,12}];

nn = ns[[configuration]];

dims = Table[Keys[nestedfps[nn]][[1+4i]],{i,0,4}];

If[!Element[configuration,Integers] || configuration<1 || configuration > Length[ns],
	Print["Choose proper run cofiguration id"];
Quit[]];

integrandsListIsotropic = (integrandsListIsotropic /. n -> nn);

eqlist = IntegrateEquations[integrandsListIsotropic];

GetFPs[dim_]:=Block[{fixedPointDictHigh,fixedPointDictLow,constants,guess},
	guess = nestedfps[nn,dim][[1]];
	constants = Association["dim"->dim,"\[Rho]Max"->Keys[nestedfps[nn,dim]][[1]]];
	fixedPointDictHigh = FindFixedPointDict[eqlist, guess, alphasHigh, {"alpha"}, constants];
	fixedPointDictLow = FindFixedPointDict[eqlist, guess, alphasLow, {"alpha"}, constants];
	Return[Association[dim->KeySort[Join[fixedPointDictLow,fixedPointDictHigh]]]];
];
fixedPointDict = UnnestAssociation[Association[nn->Association[Map[GetFPs,dims]]]];

exponents = FindExponents[integrandsListIsotropic, {}, fixedPointDict, {"N","dim","alpha","\[Rho]Max"}, Association[]];

exponents = NestAssociation[UnnestAssociation[exponents]];

If[FileExistsQ[cacheDirectory<>REGULATORDEPENDENCEFILE],
	{allFPs, allExponents} = Import[cacheDirectory<>REGULATORDEPENDENCEFILE];
	allFPs = Join[allFPs, fixedPointDict];
	allExponents = Join[allExponents, exponents],
	
	allFPs = fixedPointDict;
	allExponents = exponents];

Export[cacheDirectory<>REGULATORDEPENDENCEFILE,{allFPs,allExponents}];

PlotRegulatorDependence[specs_] := Block[{exp,nn,data,name,p,maxy,miny,maxx,minx,xmargin,ymargin},
	{exp, nn} = specs;
	data = exponents[exp,nn];
	name = ToString[exp]<>"(N="<>ToString[nn]<>")";
	
	maxy = Max[0,Max[Values[Values[data]]]];
	miny = Min[0,Min[Values[Values[data]]]];
	ymargin = (maxy-miny)*0.05;
	
	maxx = Max[Keys[Values[data]]];
	minx = Min[Keys[Values[data]]];
	xmargin = (maxy-miny)*0.05;
	
	p = ListPlot[Values[data],
		PlotLegends->PointLegend[Keys[data], LegendLabel->"d", LegendFunction->"Panel",
			LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 25]],
		AxesLabel->{"\[Alpha]",exp}, ImageSize->Full, PlotLabel->name,
		PlotRange->{{minx-xmargin, maxx+xmargin},{miny-ymargin,maxy+ymargin}},
		LabelStyle->Directive[FontFamily->"Arial",FontSize->30]];
	Export[cacheDirectory<>name<>".png",p];
];

Map[PlotRegulatorDependence,Tuples[{{e1,\[Nu],\[Eta]},Keys[exponents[\[Nu]]]}]];
];


(* ::Subsubsection:: *)
(*Other*)


If[!debug, Quit[]];
(* Anything below this line won't be executed in command line run *)


{guess, dim, \[Rho]Max} = Import[guessfile];
guess = AdjustNumberOfPoints[guess, gridSize,zNormalization];

configurations = Import[runconf];

If[!Element[configuration,Integers] || configuration<1 || configuration > Length[configurations],
	Print["Choose proper run cofiguration id"];
	Quit[]];
constants = configurations[[configuration]];

eqlist = IntegrateEquations[integrandsListIsotropic];


{fps,exps} = Import[cacheDirectory<>FIXEDPOINTSFILE];
nestedfps = NestAssociation[fps];


eqlist = eqlist/.{n->2.,vd->1.};


Length[eqlist]


0.12/0.013


guess = nestedfps[2,3][[1]];


guess[[242;;-1,2]] *=0.013;


dev = eqlist/.guess/.ConstRep[3,9.4];


ConstRep[3,0.00156]


ListPlot[dev]


fp = FindFixedPoint[eqlist,guess,ConstRep[3,9.4]];


PlotFixedPoint[fp,9.4]


fpList[0] = guess;
\[Rho]MaxList[0] = 0.1;
dims={3.,2.8,2.6,2.4,2.2,2.,1.9,1.8,1.75,1.7,1.65,1.625,1.6,1.575,1.55,1.5375,1.525,1.5,1.475,1.45,1.425,1.4,1.3,1.2,1.1,1.05,1.};


For[i = 1,  i <= Length[dims], i++,
  {newfp, new\[Rho]Max} = \[Rho]MaxScan[eqlist, fpList[i - 1], dims[[i]], \[Rho]MaxList[i - 1], 2., \[Rho]MaxList[i - 1]*0.025, 20];
  fpList[i] = newfp;
  \[Rho]MaxList[i] = new\[Rho]Max;
  Print[i];];


For[i = 15,  i <= Length[dims], i++,
  {newfp, new\[Rho]Max2} = \[Rho]MaxScan[eqlist2, fpList2[i - 1], dims[[i]], \[Rho]MaxList2[i - 1], 2., \[Rho]MaxList2[i - 1]*0.025, 20];
  fpList2[i] = newfp;
  \[Rho]MaxList2[i] = new\[Rho]Max2;
  Print[i];];


eqlist2 = eqlist[[1;;242]];


Show[PlotFixedPoint[fpList[15],\[Rho]MaxList[15]],ListPlot[Table[{i/120*0.65, KroneckerDelta[76-i]},{i,0,120}]]]


h[y_] = Table[NumericDerivatives[integrandsListIsotropic[[2,3]]/.n->1 /. regulatorReplacement/.\[Rho]-> i eps]/.fpList[10]/.ConstRep[1.2,\[Rho]MaxList[10]],{i,0,120}];


h[0.0001]
