(* ::Package:: *)

(* ::Section:: *)
(*Initialization*)


(* ::Subsection::Closed:: *)
(*Setting the project path and importing the naming module*)


(* ::Text:: *)
(*Imports ModuleNaming module and sets the working directory*)


$HistoryLength = 10;
gridSize = 60;


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


(* ::Subsection::Closed:: *)
(*Serializing the command line/debug arguments*)


If[debug,
	(* Debug run configuration *)
	task = "ON"; (* ON/ Z4 / tricritical/ regulator /lowd*)
	configuration = 55;
	ON = True;
	regulatorName = "smooth";
	strict=False;
	LPA = False,
	
	(* Release run configuration *)
	LPA = False;
	task = Null;
	anisotropySymmetry = 4;
	regulatorName = Null;
	configuration = Null;
	strict = False;
	Block[{i},
	For[i = 1, i<=Length[arguments], i++,
		If[arguments[[i]] == "-task", task = arguments[[i+1]]];
		If[arguments[[i]] == "-regulator", regulatorName = arguments[[i+1]]];
		If[arguments[[i]] == "-conf", configuration = ToExpression[arguments[[i+1]]]];
		If[arguments[[i]] == "LPA", LPA = True];
		If[arguments[[i]] == "strict", strict = True];
		];
	];
];
ON = MemberQ[ONTASKS, task];

(* Check if all required options are specified *)
If[task == Null || regulatorName == Null || (anisotropySymmetry == "ON" && nn == Null),
	Print["Specify all required options: task, regulator"];
	Quit[];
	];

If[LPA,
	zNormalization = Floor[5*gridSize/6],
	zNormalization = 0];


(* ::Subsection::Closed:: *)
(*Choosing the configurations files and importing modules*)


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

Get["FixedPoint.wl"];
Get["Plots.wl"];

(* Select regulator function *)
SelectRegulator[regulatorName];
If[strict,
	cacheDirectory = GetCacheDirectory[task, 0, regulatorName, LPA, "strict"],
	cacheDirectory = GetCacheDirectory[task, 0, regulatorName, LPA]];



(* ::Subsection::Closed:: *)
(*Importing equations*)


(* ::Text:: *)
(*Imports cached flow equations*)


If[strict,
	isotropicEquationsFile = GetEquationsFileName[task, 0, 0,"strict"],
	isotropicEquationsFile = GetEquationsFileName[task, 0, 0];];
If[FileExistsQ[isotropicEquationsFile], 
	integrandsListIsotropic = Import[isotropicEquationsFile];
	If[LPA, integrandsListIsotropic= MakeLPA[integrandsListIsotropic];],
	Print["No isotropic equations found"]; Quit[]];
If[ON,
	integrandsListAnisotropic = {},
	
	If[strict,
		anisotropicEquationsFile = GetEquationsFileName[task, 0, 1,"strict"],
		anisotropicEquationsFile = GetEquationsFileName[task, 0, 1];];
	If[FileExistsQ[anisotropicEquationsFile], 
		integrandsListAnisotropic = Import[anisotropicEquationsFile],
		Print["No anisotropic equations found"]; Quit[]];
	If[LPA, integrandsListAnisotropic= MakeLPA[integrandsListAnisotropic];]
];


(* ::Section:: *)
(*O(N) models*)


(* ::Subsection::Closed:: *)
(*O(N) dimension dependence*)


If[ task == "ON" || task == "tricritical",
{guess, dim, \[Rho]Max} = Import[guessfile];
guess = AdjustNumberOfPoints[guess, gridSize,zNormalization,\[Rho]Max];

configurations = Import[runconf];

If[!Element[configuration,Integers] || configuration<1 || configuration > Length[configurations],
	Print["Choose proper run cofiguration id"];
	Quit[]];
constants = configurations[[configuration]];

If[MemberQ[Keys[constants], "N"], 
	If[constants["N"] > 4,
		eqlistIntermediate = IntegrateEquations[integrandsListIsotropic /. n -> 4.];
		{guess, \[Rho]Max} = \[Rho]MaxScan[eqlistIntermediate, guess, dim, \[Rho]Max, 2., 0.025*\[Rho]Max, 40];
	];
	integrandsListIsotropic = (integrandsListIsotropic /. n -> constants["N"]);
	integrandsListAnisotropic = (integrandsListAnisotropic /. n -> constants["N"]);
	If[constants["N"]==1, 
		integrandsListIsotropic = integrandsListIsotropic[[1;;2]];
		guess = Join[guess[[1;;2*gridSize+1]],guess[[-1;;-1]]]]];

eqlist = IntegrateEquations[integrandsListIsotropic];
ScanCheck[scan_] := If[Length[scan]!=2, Print["Initial fixed point search failed"]; Quit[];];

scan = \[Rho]MaxScan[eqlist, guess, dim, \[Rho]Max, 2., 0.025*\[Rho]Max, 80];
ScanCheck[scan];
{newGuess, \[Rho]Max} = scan;
If[constants["alpha"]>3,
	scan2 = \[Rho]MaxScan[eqlist, newGuess, dim, 1.2*\[Rho]Max, 3., 0.0125*\[Rho]Max, 40];
	ScanCheck[scan2];
	{newGuess, \[Rho]Max} = scan2;];
	
scan3 = \[Rho]MaxScan[eqlist, newGuess, dim, 1.2*\[Rho]Max, constants["alpha"], 0.0125*\[Rho]Max, 40];
ScanCheck[scan3];
{newGuess, \[Rho]Max} = scan3;

{fixedPointDict,exponents} = FullDScan[eqlist, newGuess, dim, 0.025, \[Rho]Max, constants["alpha"]];
(*fixedPointDict = FixedPointDimScan[eqlist, newGuess, constants, dim, \[Rho]Max];
exponents = FindExponents[integrandsListIsotropic, integrandsListAnisotropic, fixedPointDict, {"dim","\[Rho]Max"}, constants];*)

If[ON,
	nestedFPD = Association[constants[["alpha"]] -> Association[constants[["N"]]->NestAssociation[fixedPointDict]]];
	fixedPointDict = UnnestAssociation[nestedFPD];];

If[ON,
	nestedExps = TransposeAssociation[Association[constants["alpha"]->TransposeAssociation[Association[constants[["N"]]->exponents]]]];
	exponents = UnnestAssociation[nestedExps],
	exponents = UnnestAssociation[exponents]];
Export[cacheDirectory<>FIXEDPOINTSFILES[configuration],{fixedPointDict,exponents}];

cache = Import[cacheDirectory<>FIXEDPOINTSFILE];
If[cache != $Failed,
	{allFPs, allExponents} = cache;
	allFPs = Join[allFPs, fixedPointDict];
	allExponents = Join[allExponents, exponents],
	
	allFPs = fixedPointDict;
	allExponents = exponents];
	
Export[cacheDirectory<>FIXEDPOINTSFILE,{allFPs,allExponents}];
];


(* ::Subsection::Closed:: *)
(*Low-d task computation*)


If[task=="lowd",
guesses = Import[guessfile];
{keyLabels, configurations} = Import[runconf];

If[!Element[configuration,Integers] || configuration < 1 || 
		configuration > Length[configurations] || configuration > Length[guesses],
	Print["Choose proper run cofiguration id"];
	Quit[]
	];

{guess, dim, \[Rho]Max} = guesses[[configuration]];
guess = AdjustNumberOfPoints[guess, gridSize,zNormalization,\[Rho]Max];

{keys, constants} = configurations[[configuration]];

If[MemberQ[Keys[constants], "N"], 
	integrandsListIsotropic = (integrandsListIsotropic /. n -> constants["N"]);
	integrandsListAnisotropic = (integrandsListAnisotropic /. n -> constants["N"])];

eqlist = IntegrateEquations[integrandsListIsotropic];

{newGuess, new\[Rho]Max} = \[Rho]MaxScan[eqlist /. n->keys[[1]], guess, dim, \[Rho]Max, 2., 0.05*\[Rho]Max, 20];
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
];



(* ::Subsection::Closed:: *)
(*Cardy line plotting*)


If[False && task == "ON"||task == "tricritical",
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
lauPoints = 1+{{0,0},{63,3},{105,13},{130,25},{159,50},{180,75},{197,100},{211,125},{222,150},{231,175},{239,200},{250,250}}/250.;
lauLine = Interpolation[lauPoints];

cacheDirectory = GetCacheDirectory[task, 0, regulatorName];

min = Min[onPoints[[All,1]]];
max = Max[onPoints[[All,1]]];
plot = Show[ListLinePlot[{onPoints[[3;;-1]]},PlotLegends->PointLegend[{"This work"},LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 20]],
		PlotMarkers->{Automatic, 7},PlotRange->All],
	Plot[2+ (d-2)*Pi^2/4,{d,min,max},PlotLegends->LineLegend[{"Cardy"},LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 20]],PlotStyle->Red],
	Plot[lauLine[d],{d,1,2},PlotLegends->LineLegend[{"Lau"},LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 20]],PlotStyle->Darker[Green]],
AxesLabel->{"d","N"}, LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 16], GridLines->Automatic,
ImageSize->Full];
Export["cardy.png",plot];
];


(* ::Subsection::Closed:: *)
(*Lau line plotting*)


If[ False,
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
lauPoints = 1+{{0,0},{63,3},{105,13},{130,25},{159,50},{180,75},{197,100},{211,125},{222,150},{231,175},{239,200},{250,250}}/250.;
lauLine = Interpolation[lauPoints];

cacheDirectory = GetCacheDirectory[task, 0, regulatorName];

plot = Show[ListLinePlot[{onPoints[[3;;-3]]},PlotLegends->PointLegend[{"This work"},LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 20]],
		PlotMarkers->{Automatic, 7}, PlotRange -> {{0.95,2.05},{0.95,2.05}}],
	Plot[2+ (d-2)*Pi^2/4,{d,1,2},PlotLegends->LineLegend[{"Cardy"},LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 20]],PlotStyle->Red],
	Plot[lauLine[d],{d,1,2},PlotLegends->LineLegend[{"Lau"},LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 20]],PlotStyle->Darker[Green]],
AxesLabel->{"d","N"}, LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 16], GridLines->Automatic,
ImageSize->Full];
Export["lau.png",plot];
];


(* ::Subsection::Closed:: *)
(*Plotting critical exponents dimension dependence *)


(*{fps,allExponents} = Import["ON_functional_smooth backup/"<>FIXEDPOINTSFILE];
expss = NestAssociation[allExponents];*)


(*{fps,allExponents} = Import["ON_functional_smooth backup/"<>FIXEDPOINTSFILE];
allExponents = NestAssociation[allExponents];
ns = {4,3,2.5,2}(*{6,4,3,2.5,2}*);
exps = {e1,\[Eta]};
allExponents = Association[Table[exps[[i]]-> KeySelect[allExponents[exps[[i]]],MemberQ[ns,#]&],{i,1,Length[exps]}]];
l1 = .025;
l2 = .001;
s = l1/2;
dash[1] = {};
dash[2] = {l1,s};
dash[3] = {l1,s,l2,s};
dash[4] = {l1, s, l2, s, l2, s};
dash[5] = {l2,s};
curves = Table[Directive[AbsoluteThickness[4],ColorData[97, "ColorList"][[i]]],{i,1,5}];
markers = Tuples[{{"\[FilledCircle]", "\[FilledSquare]", "\[FilledDiamond]", "\[FilledUpTriangle]",""},{30}}];

ticks = Association[\[Eta]-> {Table[2+0.25i,{i,0,4}],Table[0.1i,{i,0,4}]},e1-> {Table[2+0.25i,{i,0,4}],Table[0.5i,{i,0,4}]}];
PlotDimDependence[exp_] := Block[{data,min,max,p,valueRange},

data = Map[KeySelect[#,(#\[GreaterEqual]2)&]&, KeySort[allExponents[exp]]];

valueRange =  Association[e1 -> {-0.01, 1.55}, \[Nu] -> {0, 12}, \[Eta] -> {-0.001, 0.33}];
If[MemberQ[Keys[valueRange],exp],
	Map[(data[[#]] = DropOutliers[data[[#]],valueRange[exp][[1]],valueRange[exp][[2]]])&,Range[Length[data]]];
];

If[exp==e1, data["\[Infinity]"] = Association[2->0,3->1]];
If[exp==\[Eta], data["\[Infinity]"] = Association[2->0.001,3->0.001]];
min = Min[Values[data]];
max = Max[Values[data]];
If[MemberQ[Keys[valueRange],exp],
	min = valueRange[exp][[1]]; max = valueRange[exp][[2]],
	min = Min[min 1.05, 0.95 min, -0.01];
	max = Max[max 1.05, 0.95 max, 0];];
If[exp==\[Nu], min=Max[0,min]; max = Min[3,max];];
If[exp==\[Nu], data["\[Infinity]"] = Association[Table[(3.-i/100)->1./(1-i/100),{i,0,99}]]];
label = exp;
If[exp==e1, label=Superscript[\[Nu],-1]];
p = ListLinePlot[Values[data],PlotTheme->"Monochrome",PlotLegends->Placed[PointLegend[Map[("N="<>ToString[#])&,Keys[data]], LegendFunction->"Frame", LegendLayout->"Row",
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 40,Black],LegendMarkerSize->{40,20}],Bottom], LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 40],
	AxesLabel->{Style["d",FontSize->40,Black,Bold],Style[label,FontSize->40,Black,Bold]}, 
	ImageSize->Full, PlotRange->{{2,3},{min,max}},PlotMarkers->markers,
	PlotStyle->curves,Ticks->ticks[exp],GridLines->ticks[exp]];
Export[cacheDirectory<>ToString[exp]<>".png",p];];
Map[PlotDimDependence, Keys[allExponents]];*)


If[ False && task == "ON" || task == "tricritical",
{fps,allExponents} = Import["ON_functional_smooth backup/"<>FIXEDPOINTSFILE];
allExponents = NestAssociation[allExponents];
ns = {4,3,2.5,2}(*{6,4,3,2.5,2}*);
exps = {e1,e2,e3,\[Eta]};
allExponents = Association[Table[exps[[i]]-> KeySelect[allExponents[exps[[i]]],MemberQ[ns,#]&],{i,1,Length[exps]}]];
l1 = .025;
l2 = .001;
s = l1/2;
dash[1] = {};
dash[2] = {l1,s};
dash[3] = {l1,s,l2,s};
dash[4] = {l1, s, l2, s, l2, s};
dash[5] = {l2,s};
curves = Table[Directive[AbsoluteThickness[4],ColorData[97, "ColorList"][[i]]],{i,1,5}];
markers = Tuples[{{"\[FilledCircle]", "\[FilledSquare]", "\[FilledDiamond]", "\[FilledUpTriangle]",""},{30}}];

ticks = Association[\[Eta]-> {Table[2+0.25i,{i,0,4}],Table[0.1i,{i,0,4}]},
e1-> {Table[2+0.25i,{i,0,4}],Table[0.5i,{i,0,4}]},
e2-> {Table[2+0.25i,{i,0,4}],Table[-0.5i,{i,0,4}]},
e3-> {Table[2+0.25i,{i,0,4}],Table[-0.5i,{i,0,6}]}];
PlotDimDependence[exp_] := Block[{data,min,max,p,valueRange},

data = Map[KeySelect[#,(#>=2)&]&, KeySort[allExponents[exp]]];

valueRange =  Association[e1 -> {-0.01, 1.55}, \[Nu] -> {0, 12}, \[Eta] -> {-0.001, 0.33}, e2-> {-2, 0.01}, e3-> {-4, 0.01}];
If[MemberQ[Keys[valueRange],exp],
	Map[(data[[#]] = DropOutliers[data[[#]],valueRange[exp][[1]],valueRange[exp][[2]]])&,Range[Length[data]]];
];

If[exp==e1, data["\[Infinity]"] = Association[2->0,3->1]];
If[exp==\[Eta], data["\[Infinity]"] = Association[2->0.001,3->0.001]];
min = Min[Values[data]];
max = Max[Values[data]];
If[MemberQ[Keys[valueRange],exp],
	min = valueRange[exp][[1]]; max = valueRange[exp][[2]],
	min = Min[min 1.05, 0.95 min, -0.01];
	max = Max[max 1.05, 0.95 max, 0];];
If[exp==\[Nu], min=Max[0,min]; max = Min[3,max];];
If[exp==\[Nu], data["\[Infinity]"] = Association[Table[(3.-i/100)->1./(1-i/100),{i,0,99}]]];
label = exp;
If[exp==e1, label=Superscript[\[Nu],-1]];
p = ListLinePlot[Values[data],PlotTheme->"Monochrome",PlotLegends->Placed[PointLegend[Map[("N="<>ToString[#])&,Keys[data]], LegendFunction->"Frame", LegendLayout->"Row",
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 40,Black],LegendMarkerSize->{40,20}],Bottom], LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 40],
	AxesLabel->{Style["d",FontSize->40,Black,Bold],Style[label,FontSize->40,Black,Bold]}, 
	ImageSize->Full, PlotRange->{{2,3},{min,max}},PlotMarkers->markers,
	PlotStyle->curves,Ticks->ticks[exp],GridLines->ticks[exp]];
Export["~/Desktop/"<>ToString[exp]<>".png",p];];
Map[PlotDimDependence, Keys[allExponents]];
];


(* ::Section::Closed:: *)
(*Z4-symmetric anisotropy*)


(* ::Subsection::Closed:: *)
(*Z4 - anisotropy dimension dependence*)


If[ task == "Z4",
{guess, dim, \[Rho]Max} = Import[guessfile];
guess = AdjustNumberOfPoints[guess, gridSize,zNormalization,\[Rho]Max];
If[LPA,guess = LPApGuess[guess];];

configurations = Import[runconf];

If[!Element[configuration,Integers] || configuration<1 || configuration > Length[configurations],
	Print["Choose proper run cofiguration id"];
	Quit[]];
constants = configurations[[configuration]];
alpha = constants["alpha"];

eqlist = IntegrateEquations[integrandsListIsotropic];

CheckScan[scan_] := Block[{},
	If[Length[scan]!=2, Print["Initial fixed point search failed"]; Quit[];];
	{guess, \[Rho]Max} = scan;
	];


FindInitialGuess[a0_,amax_,da_]:= Block[{a},
	For[a =a0, Sign[da]*(amax-a)>0, a+=da,
	 scan = \[Rho]MaxScan[eqlist, guess, dim, \[Rho]Max, a, 0.05*\[Rho]Max, 40];
	 CheckScan[scan];
	 ];
	 scan = \[Rho]MaxScan[eqlist, guess, dim, \[Rho]Max, amax, 0.05*\[Rho]Max, 40];
	 CheckScan[scan];
];	 

scan = \[Rho]MaxScan[eqlist, guess, dim, \[Rho]Max, 2., 0.05*\[Rho]Max, 40];
CheckScan[scan];

If[alpha<1., FindInitialGuess[1., alpha, -0.2],
	If[alpha>2, FindInitialGuess[2., alpha, +0.75],
		scan = \[Rho]MaxScan[eqlist, guess, dim, \[Rho]Max, alpha, 0.05*\[Rho]Max, 40];
		CheckScan[scan];
		];
	];

fixedPointDict = FixedPointDimScan[eqlist, guess, constants, dim, \[Rho]Max];

selectedDims = {3.,2.9,2.8,2.7,2.6,2.5,2.4,2.3,2.2,2.15,2.1};
shortFixedPointDict = KeySelect[fixedPointDict,(Min[Abs[#[[1]]-selectedDims]]<=10^-4 || #[[1]]<=2.1)&&(#[[1]]>=2.-10^-4)&];

exponents = FindExponents[integrandsListIsotropic, integrandsListAnisotropic, shortFixedPointDict, {"dim","\[Rho]Max"}, constants];

fixedPointDict = UnnestAssociation[Association[alpha->fixedPointDict]];
exponents = UnnestAssociation[TransposeAssociation[Association[alpha->exponents]]];

Export[cacheDirectory<>FIXEDPOINTSFILES[configuration],{fixedPointDict,exponents}];
	
If[FileExistsQ[cacheDirectory<>FIXEDPOINTSFILE],
	{allFPs, allExponents} = Import[cacheDirectory<>FIXEDPOINTSFILE];
	allFPs = Join[allFPs, fixedPointDict];
	allExponents = Join[allExponents, exponents],
	
	allFPs = fixedPointDict;
	allExponents = exponents];
	
Export[cacheDirectory<>FIXEDPOINTSFILE,{allFPs,allExponents}];
];


(* ::Subsection::Closed:: *)
(*Z4 - anisotropy plotting regulator dependence*)


If[ task == "Z4",
{fps,allExponents} = Import[cacheDirectory<>FIXEDPOINTSFILE];
FindPMS[exps_,e_,returnVals_]:=Block[{pms=Association[], scanned, SinglePMS, j, vals, keys, noPMS},
	scanned = TransposeAssociation[exps[e]];
	SinglePMS[assoc_]:= Block[{interpolation, pos, root, found=False},
		der[f_,x_] := (f[[2;;-1]]-f[[1;;-2]])/(x[[2;;-1]]-x[[1;;-2]]);
		deri = der[Keys[assoc],Values[assoc]];
		For[i=1,i<Length[assoc]-1,i++,
			If[Sign[deri[[i]]]!=Sign[deri[[i+1]]],
				found = True;
				pos = i+1;
				Break[];
			];
		];
		If[!found, Return[{}]];
		Return[{Keys[assoc][[pos]],assoc[[pos]]}];
	];
	If[returnVals, j=2,j=1];
	keys = Table[Keys[scanned][[i]],{i,1,Length[scanned]}];
	vals = Map[SinglePMS[KeySort[#]]&,Values[scanned]];
	noPMS = Position[vals,{}];
	keys = Delete[keys, noPMS];
	vals = Delete[vals, noPMS];
	Return[AssociationThread[keys, vals[[All,j]]]];
];
allExponents = NestAssociation[allExponents];
exps = {\[Nu],\[Eta],e1,y};
(*ticks = Association[\[Eta]-> {Table[2+0.25i,{i,0,4}],Table[0.1i,{i,0,4}]},e1-> {Table[2+0.25i,{i,0,4}],Table[0.5i,{i,0,4}]}];*)
ydata = TransposeAssociation[allExponents[y]];
lowd = {2., 2.1,2.2,2.3};(*,2.2,2.3,2.4};*)
midd = {2.2,2.3,2.4,2.5};
hid = {2.5,2.6,2.7,2.8,2.9,3.};
datasets = {hid,lowd,midd};
names = {"high","low","middle"};

ranges = {{0,0.23},{-2,0.3},{-0.4,0.35}};
If[LPA, ranges = {{0,0.6},{-2.,0.35},{-0.7,0.6}}];

BigPlot[data_,legends_,axesLabels_,plotRange_:Automatic,ticks_:Automatic, PlotFunc_:ListPlot]:=Block[{p},
p = PlotFunc[data,PlotLegends->Placed[PointLegend[legends, LegendFunction->"Frame", LegendLayout->"Row",
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 40, Black],LegendMarkerSize->{40,20}],Bottom], 
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 40, Black],
	AxesLabel->{Style[axesLabels[[1]],FontSize->40,Black,Bold],Style[axesLabels[[2]],FontSize->40,Black,Bold]}, 
	ImageSize->Full, PlotRange->plotRange,PlotMarkers->{Automatic, 12},PlotStyle->AbsoluteThickness[4],
	Ticks->ticks,GridLines->ticks];
	Return[p];
];

PlotReg[dataset_,name_,range_]:= Block[{data,legends},
data = Map[SelectClosest[ydata,#]&,dataset];
legends = Map["d="<>ToString[#]&,dataset];
Export[cacheDirectory<>name<>".png",BigPlot[data,legends,{\[Alpha],Subscript[y, 4]},range]];];
MapThread[PlotReg,{datasets,names,ranges}];

];



(* ::Subsection::Closed:: *)
(*Regulator dependence analysis*)


If[task=="regulator",
{fps,exps} = Import[cacheDirectory<>FIXEDPOINTSFILE];
nestedfps = NestAssociation[fps];

nestedfps = KeySort[Association[Table[N[Keys[nestedfps][[i]]] ->nestedfps[[i]],{i,1,Length[nestedfps]}]]];
	
(*ns = {1.,2.,3.,4.,6.,8.};
nn = ns[[configuration]];*)
nn = 2.;

(*dims = Table[Keys[nestedfps[nn]][[1+4i]],{i,0,4}];*)
dims = Keys[nestedfps[nn]][[1;;-6]];
dd = dims[[configuration]];

alphasLow = Table[2-i 0.0625,{i,0,27}];
alphasHigh = Table[2+i 0.125,{i,1,8}];


If[!Element[configuration,Integers] || configuration<1 || configuration > Length[dims],
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
(*fixedPointDict = UnnestAssociation[Association[nn->Association[Map[GetFPs,dims]]]];*)
fixedPointDict = UnnestAssociation[Association[nn->GetFPs[dd]]];

exponents = FindExponents[integrandsListIsotropic, {}, fixedPointDict, {"N","dim","alpha","\[Rho]Max"}, Association[]];

exponents = UnnestAssociation[exponents];

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
(*Map[PlotRegulatorDependence,Tuples[{{e1,\[Nu],\[Eta]},Keys[exponents[\[Nu]]]}]];*)
];


If[task=="regulator",
FindPMS[exps_,e_,returnVals_]:=Block[{pms=Association[], scanned, SinglePMS, j},
	scanned = exps[e][[1]];
	SinglePMS[assoc_]:= Block[{interpolation,root},
		interpolation = Interpolation[Transpose[{Keys[assoc],Values[assoc]}]];
		Check[root = FindRoot[interpolation'[\[Alpha]],{\[Alpha],1.8},DampingFactor->0.5],root=Undefined,{InterpolatingFunction::dmval}];
		Return[{\[Alpha]/.root, interpolation[\[Alpha]/.root]}];
	];
	If[returnVals, j=2,j=1];
	Return[Association[Table[Keys[scanned][[i]]->SinglePMS[Values[scanned][[i]]][[j]],{i,1,Length[scanned]}]]];
];
{fps,allExponents} = Import[cacheDirectory<>REGULATORDEPENDENCEFILE];
allExponents = NestAssociation[allExponents];
exps = {\[Nu],\[Eta],e1};
(*ticks = Association[\[Eta]-> {Table[2+0.25i,{i,0,4}],Table[0.1i,{i,0,4}]},e1-> {Table[2+0.25i,{i,0,4}],Table[0.5i,{i,0,4}]}];*)
{fps2,exponents2} = Import[cacheDirectory<>FIXEDPOINTSFILE];
allExponents2 = NestAssociation[exponents2];

data = Association[Map[#->Quiet[FindPMS[allExponents,#,False]]&,exps[[1;;2]]]];
data2 = Association[Map[#->Quiet[FindPMS[allExponents,#,True]]&,exps]];
data3 = allExponents2[\[Eta],2.];

BigPlot[data_,legends_,axesLabels_,plotRange_:Automatic,ticks_:Automatic, PlotFunc_:ListPlot]:=Block[{p},
p = PlotFunc[data,PlotLegends->Placed[PointLegend[legends, LegendFunction->"Frame", LegendLayout->"Row",
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 40, Black],LegendMarkerSize->{40,20}],Bottom], 
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 40, Black],
	AxesLabel->{Style[axesLabels[[1]],FontSize->40,Black,Bold],Style[axesLabels[[2]],FontSize->40,Black,Bold]}, 
	ImageSize->Full, PlotRange->plotRange,PlotMarkers->{Automatic, 12},PlotStyle->AbsoluteThickness[4],
	Ticks->ticks,GridLines->ticks];
	Return[p];
];

Export[cacheDirectory<>"PMS.png",
	BigPlot[Values[data],Map[("PMS("<>ToString[#]<>")")&,Keys[data]],{"d",Subscript[\[Alpha],PMS]}, {{1.975,3.025},{1.4,2.25}}]];
Export[cacheDirectory<>"etaPMS.png",
	BigPlot[data2[\[Eta]],{},{"d",Subscript[\[Eta],PMS]},
		{{2.,3.01},{-0.005,0.305}},{{2.,2.2,2.4,2.6,2.8,3.},{0.,0.1,0.2,0.3}},ListLinePlot]];
Export[cacheDirectory<>"etaPMS2.png",
	Show[BigPlot[{data2[\[Eta]]},{Subscript[\[Eta],PMS]},{"d",\[Eta]},
		{{2.,3.01},{-0.005,0.305}},{{2.,2.2,2.4,2.6,2.8,3.},{0.,0.1,0.2,0.3}},ListLinePlot],
		ListPlot[data3,PlotLegends->Placed[PointLegend[{"\[Eta](\[Alpha]=2)"}, LegendFunction->"Frame", LegendLayout->"Row",
			LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 40, Black], LegendMarkerSize->{40,20}],Bottom],
			PlotMarkers->{"\[FilledSquare]",36},PlotStyle-> ColorData[97, "ColorList"][[2]]]
		]];
Export[cacheDirectory<>"nuPMS.png",
	BigPlot[data2[\[Nu]],{},{"d",Subscript[\[Nu],PMS]}]];
Export[cacheDirectory<>"e1PMS.png",
	BigPlot[data2[e1],{},{"d",1/\[Nu]}]];
Export[cacheDirectory<>"nu1.png",
	BigPlot[Values[allExponents[e1][[1,-5;;-3]]],Keys[allExponents[e1][[1]]][[-5;;-3]],{"\[Alpha]",1/\[Nu]}]];
Export[cacheDirectory<>"eta1.png",
	BigPlot[Values[allExponents[\[Eta]][[1,-3;;-1]]],Map["d="<>ToString[#]&,Keys[allExponents[\[Eta]][[1]]][[-3;;-1]]],{"\[Alpha]",\[Eta]},
		{{0.95,3.05},{0.25,0.322}},{{1.0,1.5,2.0,2.5,3.0},{0.26,0.28,0.30,0.32}}]];

];


(* ::Section:: *)
(*Other*)


If[!debug, Quit[]];
(* Anything below this line won't be executed in command line run *)


(* ::Subsection::Closed:: *)
(*Z4-symmetric anisotropy*)


{fps,allExponents} = Import[cacheDirectory<>FIXEDPOINTSFILE];


allExponents = Association[];
assocmap[i_] := Block[{fpss,expss,f},
	f = cacheDirectory<>FIXEDPOINTSFILES[i]; 
	If[FileExistsQ[f],
		{fpss,expss} = Import[f];
		allExponents = Join[expss,allExponents];]
		];
Map[assocmap,Range[57]]


allExponents = NestAssociation[allExponents];
exps = {\[Nu],\[Eta],e1,y};
(*ticks = Association[\[Eta]-> {Table[2+0.25i,{i,0,4}],Table[0.1i,{i,0,4}]},e1-> {Table[2+0.25i,{i,0,4}],Table[0.5i,{i,0,4}]}];*)


PlotSelected[exps_,dims_,title_,range_:{}]:= Block[{expons={},SelectStyle,p, legend, yrange, min, max},
SelectStyle[length_]:=Block[{},
	curves = Join[
	Table[Directive[AbsoluteThickness[0.01],ColorData[97, "ColorList"][[i]]],{i,1,length}],
	Table[Directive[AbsoluteThickness[2],ColorData[97, "ColorList"][[i]]],{i,1,length}]];
	markers = Join[Tuples[{Table["\[FilledCircle]",{i,1,length}],{10}}],Tuples[{Table["\[Cross]",{i,1,length}],{1}}]];];

expons = Re[Map[Values[KeySelect[TransposeAssociation[allExponents[#]],(Min[Abs[#-dims]]<0.001)&]]&,exps]];

If[Length[exps]==2,
	SelectStyle[Length[dims]];
	legend = Map["d="<>ToString[#]&,Sort[dims,#1>#2&]],
	If[Length[dims]<=2,
		expons = Transpose[expons];
		SelectStyle[Length[exps]];
		legend = Map[ToString[#]&,exps];,
		
		Print["Either exps or dims should have length 2"];
		Return[];
	];
];
expons = Flatten[expons];
If[Length[expons]==0,Return[]];
If[Length[range]==0,
min = Min[expons];
max = Max[expons];

yrange = {Min[min-(max-min)*0.05,-0.025],Max[max+(max-min)*0.05,0.025]},
yrange=range;
];
p = ListLinePlot[expons, PlotStyle->curves, PlotMarkers->markers, LabelStyle->Directive[Bold, 14], 
	PlotLegends->LineLegend[legend,LegendMarkerSize->{30,20}], PlotLabel->Style[ToString[title],FontSize->20],
	ImageSize->Large, AxesLabel->{"\[Alpha]"},PlotRange->{{0,6.1},yrange}];
	
Export["~/Desktop/"<>ToString[title]<>".png",p];
];


PlotSelected[{y1,e1},{2.,2.025,2.05,2.1,2.15,2.2,2.3},"First EV"]


PlotSelected[{y1,y2,y3,y4,y5,e1,e2},{2.},"Zoom In (d=2)",{1.025,-1.025}]


Map[PlotSelected[{y1,y2,y3,y4,y5,e1,e2},{#},"EVs (d="<>ToString[PaddedForm[#,{4,3}]]<>")",{2.025,-2}]&,{2.,2.025,2.05,2.1,2.15,2.2,2.3,2.4,2.5, 2.6,3.}]



SelectClosest[allExponents[y2][[-1]],2.2]


SelectClosest[allExponents[y2,1.8],2.]


PlotSelected[{\[Eta],\[Eta]},{2.,2.025,2.05,2.1,2.15,2.2,2.3},"eta"]


EncodeEigenvector[ev_]:=Block[{eval0 = {v,w,t}, labels={i,j},f,vec={},k},
If[Length[ev] == 7*gridSize+3,labels = {v,w,zs,zp,t,zs1,zp1};];
If[Length[ev] == 3*gridSize+1,labels = {v,zs,zp};];

k=1;
For[i=1,i<=Length[labels],i++,
	start = 1;
	f= labels[[i]];
	If[MemberQ[eval0,f],start=0];
	For[j=start,j<=gridSize,j++,
		AppendTo[vec,f[j]->Re[ev[[k]]]];
	];
];
Return[vec];
];


(* ::Subsection:: *)
(*Ising model in 2D*)


(* ::Subsubsection::Closed:: *)
(*Stability matrix*)


{fps,exps} = Import["/home/andrzej/FunctionalRenormalization/fixed_points1.wdx"];


packed = SelectClosest[NestAssociation[fps][1.],2.];
\[Rho]Max = Keys[packed][[1]];
fp=packed[[1]];


functions = {v,z,yy};
fp2 = Flatten[Map[Table[#[i]->(#[2*i]/.fp),{i,If[BooleanQ[#==z],1,0],gridSize}]&,functions]];
fp2 = Join[fp2,{\[Eta]->(\[Eta]/.fp)}];


AppendTo[integrandsListIsotropic[[1]],integrandsListIsotropic[[1,3]]/.\[Rho]->0];
AppendTo[integrandsListIsotropic[[2]],integrandsListIsotropic[[2,3]]/.\[Rho]->0];
AppendTo[integrandsListIsotropic[[3]],integrandsListIsotropic[[4,3]]];
integrandsListIsotropic = Drop[integrandsListIsotropic,{4}];


integrands = integrandsListIsotropic/.n->1.;


eqs = IntegrateEquations[integrandsListIsotropic/.n->1];


\[Rho]Max2 = 7.5;


fp2 = FindFixedPoint[eqs, fp, ConstRep[2,\[Rho]Max2,2]];


sm = StabilityMatrix[integrands,fp2,ConstRep[2,\[Rho]Max2,2]];


sm2 = sm[[1;;241,1;;241]];


es2 = EigenSys[sm2];


es2[[1;;6,1]]


(* ::Subsubsection::Closed:: *)
(*Flow fitting*)


data = Transpose[Import["~/Desktop/flow.csv"][[2;;-1]]];


start = 100;
cutoff = 2650;
sparsity = 25;
s = data[[1]][[start;;cutoff;;sparsity]];
k = data[[2]][[start;;cutoff;;sparsity]];
u = data[[3]][[start;;cutoff;;sparsity]];
Z = data[[8]][[start;;cutoff;;sparsity]];
sfull = data[[1]][[1;;-1;;2 sparsity]];
zfull = data[[8]][[1;;-1;;2 sparsity]];
kfull = data[[2]][[1;;-1;;2 sparsity]];
ufull = data[[3]][[1;;-1;;2 sparsity]];
(*ListPlot[Transpose[{s,k}],PlotRange->All]
ListPlot[Transpose[{s,Z}],PlotRange->All]
ListPlot[Transpose[{s,u}],PlotRange->All]*)


mu = NonlinearModelFit[Transpose[{s,u}], u0 + A Cos[\[Kappa] x + \[Phi]] Exp[-\[Omega] x],{u0,A,\[Kappa],\[Phi],\[Omega]},x];
mu["ParameterTable"]
limu = Limit[mu[x],x->Infinity];
Show[ListLogPlot[Transpose[{sfull,Abs[ufull-limu]}]],LogPlot[Abs[mu[x]-limu],{x,0,7.5},PlotPoints->200],PlotRange->All]


s[[1]]


s[[-1]]


mz = NonlinearModelFit[Transpose[{s,Z}],z0 + A  Exp[-\[Omega] x],{z0,A,\[Kappa],\[Phi],\[Omega]},x];
mz["ParameterTable"]
limz = Limit[mz[x],x->Infinity];
Show[ListLogPlot[Transpose[{sfull,Abs[zfull-limz]}]],LogPlot[Abs[mz[x]-limz],{x,0,7.5},PlotPoints->200],PlotRange->All]


Exp[5.]


mk = NonlinearModelFit[Transpose[{s,k}],k0 + A Cos[\[Kappa] x + \[Phi]] Exp[-\[Omega] x],{k0,A,\[Kappa],\[Phi],\[Omega]},x];
mk["ParameterTable"]
limk = Limit[mk[x],x->Infinity];
Show[ListLogPlot[Transpose[{sfull,Abs[kfull-limk]}]],LogPlot[Abs[mk[x]-limk],{x,0,7.5},PlotPoints->200],PlotRange->All]


(* ::Subsubsection::Closed:: *)
(*Scaling term fitting*)


scalingdata = Import["~/Desktop/scaling.csv"];


\[Phi]v60 = {0.8827095034203433`,0.8708502261534294`,0.8587396917620681`,0.8463636219411453`,0.8337061796964774`,0.8207496962929695`,0.8074743348326565`,0.7938576716197199`,0.7798741704920028`,0.7654945164666294`,0.7506847624467241`,0.7354052245641369`,0.7196090351426094`,0.7032402226963413`,0.6862311284103333`,0.6684988758265061`,0.649940463875216`,0.6304258159102967`,0.6097877219997297`,0.587806932920095`,0.5641894580731178`,0.5385308904366447`,0.5102582733577723`,0.478531254677223`,0.44206531447446673`,0.3987958519557004`,0.34519087106109875`,0.2747106857992594`,0.17395007304154927`,0.011696805717928563`,-0.29372973651091827`,-0.9057669898424977`,-1.639685898986093`,-2.055040680754708`,-2.2640110198971612`,-2.3869809418488117`,-2.4703537308862624`,-2.532946056239073`,-2.5834723074278587`,-2.6264449333688487`,-2.6644101862053424`,-2.698902497751375`,-2.7308944748754804`,-2.7610269717009963`,-2.7897348135122932`,-2.8173194356506763`,-2.8439929033230595`,-2.8699057279867692`,-2.895165114541177`,-2.919847341576633`,-2.9440064174927074`,-2.967680291564982`,-2.9908954043231013`,-3.0136700706929003`,-3.03601701499157`,-3.057945270981683`,-3.0794615952689184`,-3.100571503005896`,-3.1212800324583667`,3.141592653589793`};
\[Phi]z60 = {3.116447773714084`,3.1168761770453166`,3.117203558155534`,3.1174223031773183`,3.1175243398906556`,3.1175011357966915`,3.1173436752983372`,3.117042433670609`,3.1165873490903313`,3.1159677928568694`,3.115172537771974`,3.114189724620903`,3.1130068267073128`,3.111610612421656`,3.109987105865237`,3.1081215456164566`,3.105998341811612`,3.103601031831198`,3.100912235041334`,3.0979136072493376`,3.094585795807605`,3.0909083966577535`,3.0868599150672367`,3.082417732399838`,3.077558082006906`,3.072256038261206`,3.0664855239141517`,3.0602193423751465`,3.0534292432170833`,3.0460860312213076`,3.0381597315774447`,3.029619826392108`,3.020435580305019`,3.010576475520333`,3.0000127785327755`,2.9887162616486456`,2.9766611011857944`,2.963824969800647`,2.9501903312649405`,2.9357459305681037`,2.9204884489433818`,2.9044242614368936`,2.88757119451553`,2.8699601359205267`,2.8516363049005404`,2.8326599583883323`,2.8131063011517092`,2.793064399713679`,2.7726349812317066`,2.751927129923131`,2.731054060180138`,2.7101283159412364`,2.6892568780000365`,2.668536714212989`,2.6480512587022784`,2.627868160490826`,2.6080384343399907`,2.58859692949386`,2.569563855164189`,2.55094697361448`};
\[Phi]v120 = {0.8885490129479291`,0.8827095034203433`,0.876810433988333`,0.8708502261534294`,0.8648272247986086`,0.8587396917620681`,0.852585799130116`,0.8463636219411453`,0.8400711301443388`,0.8337061796964774`,0.8272665026777215`,0.8207496962929695`,0.8141532106063208`,0.8074743348326565`,0.8007101819821008`,0.7938576716197199`,0.786913510463003`,0.7798741704920028`,0.7727358641902947`,0.7654945164666294`,0.7581457327243755`,0.7506847624467241`,0.7431064575443378`,0.7354052245641369`,0.7275749696775275`,0.7196090351426094`,0.7115001256596971`,0.7032402226963413`,0.6948204844290949`,0.6862311284103333`,0.6774612933864892`,0.6684988758265061`,0.6593303356079516`,0.649940463875216`,0.6403121042216435`,0.6304258159102967`,0.6202594646300268`,0.6097877219997297`,0.5989814492774372`,0.587806932920095`,0.5762249289298974`,0.5641894580731178`,0.5516462732199534`,0.5385308904366447`,0.52476603278392`,0.5102582733577723`,0.4948935713520505`,0.478531254677223`,0.46099578662666285`,0.44206531447446673`,0.4214554522226191`,0.3987958519557004`,0.37359560308264544`,0.34519087106109875`,0.31266349872891347`,0.2747106857992594`,0.22942964518973355`,0.17395007304154927`,0.10378841269187176`,0.011696805717928563`,-0.11431914000569683`,-0.29372973651091827`,-0.5531518172778075`,-0.9057669898424977`,-1.300296788222306`,-1.639685898986093`,-1.8851823379616521`,-2.055040680754708`,-2.1752696401050975`,-2.2640110198971612`,-2.3323451163751776`,-2.3869809418488117`,-2.4320890968026725`,-2.4703537308862624`,-2.503562923773089`,-2.532946056239073`,-2.5593724713889574`,-2.5834723074278587`,-2.605712331047141`,-2.6264449333688487`,-2.6459406367709746`,-2.6644101862053424`,-2.6820198920433045`,-2.698902497751375`,-2.715165016499457`,-2.7308944748754804`,-2.746162185881948`,-2.7610269717009963`,-2.775537625384264`,-2.7897348135122932`,-2.803652563090022`,-2.8173194356506763`,-2.830759463507723`,-2.8439929033230595`,-2.857036848035907`,-2.8699057279867692`,-2.8826117246079725`,-2.895165114541177`,-2.9075745579337737`,-2.919847341576633`,-2.9319895852043936`,-2.9440064174927074`,-2.9559021269152805`,-2.967680291564982`,-2.9793438912230914`,-2.9908954043231013`,-3.002336891958165`,-3.0136700706929003`,-3.024896375636033`,-3.03601701499157`,-3.0470330171186464`,-3.057945270981683`,-3.0687545607541193`,-3.0794615952689184`,-3.0900670327397712`,-3.100571503005896`,-3.1109756167746405`,-3.1212800324583667`,-3.1314852448130672`,3.141592653589793`};
\[Phi]z120 = {3.116447773714084`,3.1166741473575756`,3.1168761770453166`,3.1170529601083397`,3.117203558155534`,3.1173270041326395`,3.1174223031773183`,3.117488432320947`,3.1175243398906556`,3.1175289448020345`,3.1175011357966915`,3.1174397706417523`,3.1173436752983372`,3.1172116430612484`,3.117042433670609`,3.116834772395718`,3.1165873490903313`,3.1162988172194552`,3.1159677928568694`,3.1155928536528688`,3.115172537771974`,3.114705342800085`,3.114189724620903`,3.1136240962615385`,3.1130068267073128`,3.112336239685932`,3.111610612421656`,3.1108281743601562`,3.109987105865237`,3.1090855368888604`,3.1081215456164566`,3.107093157089931`,3.105998341811612`,3.1048350143327337`,3.103601031831198`,3.1022941926840653`,3.100912235041334`,3.099452835408829`,3.0979136072493376`,3.096292099612738`,3.094585795807605`,3.0927921121287816`,3.0909083966577535`,3.0889319281550813`,3.0868599150672367`,3.0846894946731833`,3.082417732399838`,3.080041621339507`,3.077558082006906`,3.0749639623783516`,3.072256038261206`,3.06943101404772`,3.0664855239141517`,3.0634161335331225`,3.0602193423751465`,3.056891586683669`,3.0534292432170833`,3.0498286338606473`,3.0460860312213076`,3.042197665328869`,3.0381597315774447`,3.03396840005164`,3.029619826392108`,3.0251101643646114`,3.020435580305019`,3.0155922696191295`,3.010576475520333`,3.0053845101888585`,3.0000127785327755`,2.9944578047217156`,2.9887162616486456`,2.9827850034506427`,2.9766611011857944`,2.9703418817177365`,2.963824969800647`,2.957108333283784`,2.9501903312649405`,2.9430697649146773`,2.9357459305681037`,2.928218674537964`,2.9204884489433818`,2.9125563676753288`,2.9044242614368936`,2.8960947306096907`,2.88757119451553`,2.8788579354745925`,2.8699601359205267`,2.8608839067325915`,2.8516363049005404`,2.8422253386649947`,2.8326599583883323`,2.8229500316219553`,2.8131063011517092`,2.8031403252279166`,2.793064399713679`,2.7828914625016616`,2.7726349812317066`,2.7623088260578745`,2.751927129923131`,2.741504139458175`,2.731054060180138`,2.72059090008248`,2.7101283159412364`,2.6996794666874364`,2.6892568780000365`,2.6788723218615274`,2.668536714212989`,2.65826003308358`,2.6480512587022784`,2.637918336185403`,2.627868160490826`,2.6179065824954884`,2.6080384343399907`,2.5982675715551093`,2.58859692949386`,2.5790285889974154`,2.569563855164189`,2.560203328949507`,2.55094697361448`,2.5417942892211376`};
\[Phi]v = \[Phi]v120;
\[Phi]z = \[Phi]z120;


kv =scalingdata[[3;;-1]][[All,3]];kz = scalingdata[[3;;-1]][[All,4]];
sv = scalingdata[[3;;-1]][[All,1]];sz = scalingdata[[3;;-1]][[All,2]];
\[Rho]points = scalingdata[[1,1]]*Range[Length[kv]];
lowcut = 61;
highcut= 76;
\[Rho]pointsv = Join[\[Rho]points[[1;;lowcut]],\[Rho]points[[highcut;;-1]]];
\[Rho]pointsz = \[Rho]points[[2;;-1]];
\[Phi]v2 = Join[\[Phi]v[[1;;lowcut]],\[Phi]v[[highcut;;-1]]];
\[Phi]z2 = \[Phi]z[[2;;-1]];
kv2 = Join[kv[[1;;lowcut]],kv[[highcut;;-1]]];
sv2 = Join[sv[[1;;lowcut]],sv[[highcut;;-1]]];
sz2 = sz[[2;;-1]];


fitv = NonlinearModelFit[Transpose[{\[Phi]v2,sv2}], ((-\[Kappa]^2+\[Omega]^2) Cos[\[Phi]+\[Theta]]+2 \[Kappa] \[Omega] Sin[\[Phi]+\[Theta]])/(\[Omega] Cos[\[Phi]+\[Theta]]+\[Kappa] Sin[\[Phi]+\[Theta]]),{{\[Omega],1.5},{\[Kappa],0.3},{\[Phi],3}},{\[Theta]}]
Show[ListPlot[Transpose[{\[Rho]points,sv}]],ListLinePlot[Table[{\[Rho]points[[i]],fitv[\[Phi]v[[i]]]},{i,1,Length[\[Rho]points]}]]]
fitv["ParameterTable"]


fitz = NonlinearModelFit[Transpose[{\[Phi]z2,sz2}], ((-\[Kappa]^2+\[Omega]^2) Cos[\[Phi]+\[Theta]]+2 \[Kappa] \[Omega] Sin[\[Phi]+\[Theta]])/(\[Omega] Cos[\[Phi]+\[Theta]]+\[Kappa] Sin[\[Phi]+\[Theta]]),{{\[Omega],1.6},{\[Kappa],0.2},{\[Phi],3}},{\[Theta]}]
Show[ListPlot[Transpose[{\[Rho]points,sz}]],ListLinePlot[Table[{\[Rho]points[[i]],fitz[\[Phi]z[[i]]]},{i,1,Length[\[Rho]points]}]]]
fitz["ParameterTable"]


k0[\[Theta]_,\[Phi]_] := (\[Kappa]^2 (\[Kappa]^2+\[Omega]^2))/(\[Omega] Cos[\[Theta]+\[Phi]]+\[Kappa] Sin[\[Theta]+\[Phi]])^2/.{\[Kappa]->0.25,\[Omega]->1.6};
\[Phi]0=2.53;
fitz = Transpose[{\[Rho]points,k0[\[Phi]z,\[Phi]0]}];
fitv = Transpose[{\[Rho]points,k0[\[Phi]v,\[Phi]0]}];
Show[ListPlot[{Transpose[{\[Rho]points,kv}],Transpose[{\[Rho]points,kz}]}],ListLinePlot[{fitv,fitz}],PlotRange->{0,0.5}]


modz = NonlinearModelFit[Transpose[{\[Phi]z,kz}][[2;;-1]],(\[Kappa]^2 (\[Kappa]^2+\[Omega]^2))/(\[Omega] Cos[\[Theta]+\[Phi]]+\[Kappa] Sin[\[Theta]+\[Phi]])^2,{\[Kappa],\[Omega],\[Phi]},\[Theta]];
modz["ParameterTable"]
modv = NonlinearModelFit[Transpose[{\[Phi]v,kv}][[2;;-1]],(\[Kappa]^2 (\[Kappa]^2+\[Omega]^2))/(\[Omega] Cos[\[Theta]+\[Phi]]+\[Kappa] Sin[\[Theta]+\[Phi]])^2,{\[Kappa],\[Omega],\[Phi]},\[Theta]];
modv["ParameterTable"]


(* ::Subsubsection::Closed:: *)
(*Ising Eigenvalues*)


{fps,exps} = Import["/home/andrzej/FunctionalRenormalization/fixed_points1.wdx"];


nested = NestAssociation[exps];
nested2 = TransposeAssociation[nested][1.];
nested3 = TransposeAssociation[KeySelect[TransposeAssociation[nested2],#>=2&]];


selected = KeySelect[nested3,MemberQ[{e2,e3,e4},#]&];


ee2 = KeySort[Join[KeySelect[selected[e2],#<=2.4&],KeySelect[selected[e4],(2.6>=#&& #>2.4)&],KeySelect[selected[e3],2.6<#&]]];
ee3 = KeySort[Join[KeySelect[selected[e3],#<=2.4&],KeySelect[selected[e2],(2.6>=#&& #>2.4)&],KeySelect[selected[e2],2.6<#&]]];
ee4 = KeySort[Join[KeySelect[selected[e4],#<=2.4&],KeySelect[selected[e3],(2.6>=#&& #>2.4)&],KeySelect[selected[e4],2.6<#&]]];


split= Flatten[Map[{Re[#],Im[#]}&,{ee2, ee3,ee4}]];


split= Flatten[Map[{Re[#],Im[#]}&,Values[selected]]];


split= Flatten[Re[Values[selected]]];


curves = Table[Directive[AbsoluteThickness[4],ColorData[97, "ColorList"][[Ceiling[i/2]]]],{i,1,6}];
markers = Tuples[{{"\[FilledCircle]", "\[FilledDiamond]", "\[FilledCircle]", "\[FilledDiamond]", "\[FilledCircle]", "\[FilledDiamond]"},{13}}];




p = ListPlot[split,PlotMarkers->markers,PlotStyle->curves, AxesLabel->{"d","Eigenvalue"},LabelStyle->Large, Ticks->Automatic,GridLines->Automatic,ImageSize->Large];
Export["ising.png", p]



(* ::Subsubsection:: *)
(*O(N)-models*)


allExps = Association[];
allFPS = Association[];


ReadCache[i_] := Block[{fps,exps},
	{fps,exps} = Import["/home/andrzej/FunctionalRenormalization/ON_functional_smooth/fixed_points"<>ToString[i]<>".wdx"];
	allExps = Join[allExps,exps];
	allFPS = Join[allFPS,fps];
];


Map[ReadCache,Range[147]];


{fps,exps} = Import["/home/andrzej/FunctionalRenormalization/ON_functional_smooth/fixed_points.wdx"];


fp = SelectClosest[NestAssociation[fps][1.,1.],2.][[1]];


p=ListPlot[{es[[1,2]],Re[es[[2,2]]],Im[es[[2,2]]]},PlotLegends->{"Dominant","Re(Subdominant)","Im(Subdominant)"}];
Export["~/Desktop/eigenvectors.png",p];


es[[2,2]]


SelectClosest[NestAssociation[fps][1.,1.],2.]


sm = StabilityMatrix[integrandsListIsotropic[[1;;2]]/.n->1,fp,ConstRep[2., 5.1288,1.]];


es = EigenSys[sm];


ev2[\[Phi]_] := Abs[es[[2,2]]] Cos[Arg[es[[2,2]]]+\[Phi]]


p = Plot[es[[1,2]].ev2[\[Phi]]/Norm[ev2[\[Phi]]],{\[Phi],0,2Pi}];
Export["~/Desktop/overlap.png",p];


Norm[es[[2,2]]]


es[[1;;3,1]]


Keys[allExps]


a = 1.;
exp = \[Eta];
e1s= KeySelect[allExps,MatchQ[{exp,a,n_,d_}]];
data = NestAssociation[e1s][exp][a];


exp = e1;
n=1.;
e1s= KeySelect[allExps,MatchQ[{exp,aa_,n,dd_}]];
pmsdata = TransposeAssociation[TransposeAssociation[NestAssociation[e1s][exp]][n]];


ds= {1.35,1.375,1.4,1.45,1.5};
ListPlot[Map[SelectClosest[pmsdata,#]&,ds],PlotLegends->ds,PlotRange->{{1,3},{0.6,.9}}]


ds= {1.35,1.375,1.4,1.45,1.5};
ListPlot[Map[SelectClosest[pmsdata,#]&,ds],PlotLegends->ds]


pmsdata[1.4]


ListPlot[Values[data],PlotLegends->Keys[data],ImageSize->Large,PlotRange->{{0.95,3.05},{0,1.}}]


ListPlot[Values[data],PlotLegends->Keys[data],ImageSize->Large,PlotRange->{{0.95,3.05},{0,1.7}}]


MatchQ[{e1,1.`,1.`,3.`},{e1,b_,1.,d_}]


Keys[fps]


exps


exps[[2]]


exps = KeySelect[allExps,MatchQ[#,{a_,2.,c_,e1}]&];


newExps = NestAssociation[Association[Re[KeyMap[{#[[3]],#[[1]]}&,exps]]]];


ListPlot[Values[newExps],PlotLegends->Keys[newExps]]


newExps


newExps


nested = NestAssociation[exps];


Keys[nested[[1]]]


Length[exps]


Keys[allExps]


KeySelect


exps = NestAssociation[allExps];


Block[{ii=60},
Print[Keys[allExps][[ii]]];
Print[Keys[allFPS][[ii]]];
fp = allFPS[[ii]];
{dim,\[Rho]M} = Keys[allFPS][[ii,3;;4]];
Print[allExps[[ii]]]];


data = TransposeAssociation[nested[e1]][1.][[1;;4]];


PlotFixedPoint[fp,\[Rho]M]


{fps,exps} = Import["/home/andrzej/FunctionalRenormalization/ON_functional_smooth/fixed_points1.wdx"];


exps


PlotFixedPoint[NestAssociation[fps][1.,1.][[-1,1]],Keys[NestAssociation[fps][1.,1.][[-1]]][[1]]]


nested = NestAssociation[exps];





es = EigenSys[sm];


ListPlot[es[[1;;3,2]],PlotRange->All]


es[[1;;3,1]]


g[y_]=D[NumericDerivatives[integrandsListIsotropic[[1,3]]/.regulatorReplacement/.\[Rho]->5 eps],d]/.fp/.ConstRep[dim,\[Rho]M,1.]/.n->1;


eqs = IntegrateEquations[integrandsListIsotropic/.n->1];


der = D[eqs,d]/.fp/.ConstRep[dim,\[Rho]M,1];


ListPlot[der]


integrands = integrandsListIsotropic[[1;;2]]/.n->1;


eqs  = IntegrateEquations[integrands];


fp0 = Join[fp[[1;;2 gridSize+1]],fp[[-1;;-1]]];


NewFP[old_, \[Rho]Max_, dim_, dd_] := Block[{der, sm, dG, new, new\[Rho]},
	der = D[eqs,d]/.old/.ConstRep[dim,\[Rho]Max,1];
	der = Drop[der,{gridSize+2}];
	sm = StabilityMatrix[integrandsListIsotropic[[1;;2]]/.n->1,old,ConstRep[dim,\[Rho]Max,1.]];
	Print[EigenSys[sm][[1;;3,1]]];
	dG = -Inverse[sm].(der/.fp/.ConstRep[dim,\[Rho]M,1]);
	AppendTo[dG,0.]; (* Append no change to \[Eta] *)
	
	new = old;
	new[[All,2]] = old[[All,2]] - dd dG;
	
	{new,new\[Rho]} = Rescale\[Rho][new,\[Rho]Max];
	
	new = FindFixedPoint[eqs,new,ConstRep[dim-dd,new\[Rho],1.],4];
	Return[{new, new\[Rho], dim-dd}];
];


fplist2 = {};


fp0


Length[eqs]


Length[fp0]


fplist2 = Drop[fplist2,{-1}];


AppendTo[fplist2, NewFP[fp0,\[Rho]M, dim,0.01]];


AppendTo[fplist2, NewFP[fplist2[[-1,1]],fplist2[[-1,2]],fplist2[[-1,3]],0.0075]];
PlotFixedPoint[fplist2[[-1,1]],fplist2[[-1,2]]]
IsFixedPoint[eqs,fplist2[[-1,1]],ConstRep[fplist2[[-1,3]],fplist2[[-1,2]],1.]]


fplist = {};


fplist = Drop[fplist,{-1}];


AppendTo[fplist, NewFP[fp,\[Rho]M, dim,0.01]];


AppendTo[fplist, NewFP[fplist[[-1,1]],fplist[[-1,2]],fplist[[-1,3]],0.0025]];
PlotFixedPoint[fplist[[-1,1]],fplist[[-1,2]]]
Print[zp[60]/.fplist[[-1,1]]]


IsFixedPoint[eqs,fplist[[-1,1]],ConstRep[fplist[[-1,3]],fplist[[-1,2]],1.]]


h[y_] = NumericDerivatives[integrandsListIsotropic[[All,3]]/.n->1/.regulatorReplacement/.\[Rho]->1 eps]/.fplist[[-1,1]]/.ConstRep[fplist[[-2,3]],fplist[[-2,2]],1.];


GetIntegral[h]


NIntegrate[h[y],{y,0.001,25}]


h[z_] := FindRoot[D[z y^2 + r[y^2]/.regulatorReplacement/.\[Alpha]->1.,y]==0,{y,1.5}]


{fp3,\[Rho]M3} = Rescale\[Rho][fp2,\[Rho]M];


PlotFixedPoint[fp3,\[Rho]M3]


NIntegrate[f[y],{y,0.000001,7}]


GetIntegral[f]


nested = NestAssociation[exps];


ListPlot[nested[\[Nu],1.,1.]]


zint = Interpolation[Table[{i 8/gridSize, zs[i]/.fp[[1]]/.zs[i_]->1.},{i,0,gridSize}],InterpolationOrder->6];


Plot[{zint[\[Rho]],zint'[\[Rho]]},{\[Rho],0,8},PlotRange->All]


SelectRegulator["exponential"]


IntegralConfigurations={{0.0,0.1`,3},{0.1`,1.`,5},{1.`,2.`,5},{2.`,3.`,4},{3.`,5.`,3}};


integrands = integrandsListIsotropic[[All,3]] /. \[Rho]-> 57 eps /.n->1 /.regulatorReplacement // NumericDerivatives;
integrands = integrands/.ConstRep[dim,\[Rho]M,1.]/.fp;
i[q_] = integrands/.y->q;
i2[q_] = Normal[Series[integrands,{y,0,12}]]/.y->q;
conf = {1.,2.,5};
Norm[GetIntegral[i]-NIntegrate[i[q],{q,0.001,5}]]


NIntegrate[i2[q],{q,conf[[1]],conf[[2]]}]


GLIntegral[i,conf[[1]],conf[[2]],conf[[3]]]


simp[q_] = Simplify[i[q]];


integrands[[1;;2]]


Simplify[integrands[[1;;2]]]


Plot[integrands,{y,0.0001,2}]


Plot[simp/.ConstRep[1.4,8,1.6]/.fp[[1]],{y,0,0.01}]


PlotFixedPoint[nested[[-8]]]


nested = NestAssociation[exps];


data = nested[e1][1.85];
ListLinePlot[Values[data],PlotLegends->Keys[data],PlotRange->{{1,3},{0,1.7}}]


Keys[nested[e1]]


nested = DeleteCases[nested,Missing[KeyAbsent,x_]];


GetDim[nn_,d_,e_:e1] := Block[{assoc},
AddExp[alpha_] := Block[{},
	If[MemberQ[Keys[nested[e,alpha]],nn] && Min[Abs[Keys[nested[e,alpha,nn]]-d]]<10^-4,
		Return[Re[SelectClosest[nested[e,alpha,nn],d]]],
		Return[Null];
	];
];
assoc = Map[#->AddExp[#]&,Keys[nested[e]]];
assoc = DeleteCases[assoc,x_->Null];
Return[Association[assoc]];
];


as = Table[1.+0.075 * 4 i,{i,0,3}];


ListPlot[Map[nested[e1,#,1.]&,as],PlotLegends->as]


pmsdata = Map[#->GetDim[1.1,#,e1]&,{3.,2.,1.75,1.65,1.625,1.6125}];
ListPlot[Values[pmsdata],PlotLegends->Keys[pmsdata]]


nested[e1]


data = nested[e1][1.85];
ListLinePlot[Values[data],PlotLegends->Keys[data],PlotRange->{{1,3},{0,1.7}}]


ConstRep[a,b,c]


data = nested[\[Eta]][2.];
ListLinePlot[Values[data],PlotLegends->Keys[data],PlotRange->{{1,3},{0,1.7}}]


Keys[nested[e1]]


Keys[nested]


nested = NestAssociation[fps];


nested = KeySelect[nested,#>=2&];


PlotKappa[fps_] := Block[{len,ds,lim=4,\[Rho]Maxs,fixedpoints,vs,\[Pi]s,\[Sigma]s,plot,dir,NThread,data,stars,cardy,cardypoints},
AtKappa[nn_,dim_]:=Block[{fp, \[Rho]max, vfunc, \[Sigma]func,\[Pi]func, \[Rho]0},
	fp = fps[nn,dim][[1]];
	\[Rho]max = Keys[fps[nn,dim]][[1]];
	vfunc = Interpolation[Table[{i \[Rho]max/gridSize, v[i]/.fp},{i,0,gridSize}]];
	\[Sigma]func = Interpolation[Table[{i \[Rho]max/gridSize, z[i]/.fp/.z[0]->1},{i,0,gridSize}]];
	\[Pi]func = Interpolation[Table[{i \[Rho]max/gridSize, (z[i]-2 i \[Rho]max/gridSize yy[i])/.fp/.z[0]->1},{i,0,gridSize}]];
	\[Rho]0 = FindRoot[vfunc[\[Rho]],{\[Rho],\[Rho]max}][[1,2]];
	
	Return[{vfunc'[\[Rho]0],\[Sigma]func[\[Rho]0],\[Pi]func[\[Rho]0],vfunc[0],vfunc'[0]}]];
NThread[nn_] := AssociationMap[AtKappa[nn,#]&,Keys[fps[nn]]];
data = AssociationMap[NThread,Keys[fps]];

l1 = .025;
l2 = .001;
s = l1/2;
dash[1] = {};
dash[2] = {l1,s};
dash[3] = {l1,s,l2,s};
dash[4] = {l1, s, l2, s, l2, s};
dash[5] = {l2,s};
curves = Table[Directive[AbsoluteThickness[6],ColorData[97, "ColorList"][[i]]],{i,1,6}];
markers = Tuples[{{"\[FilledCircle]", "\[FilledSquare]", "\[FilledDiamond]", "\[FilledUpTriangle]","\[FilledDownTriangle]"},{30}}];
plot[i_,title_,label_]:=Block[{p,funcs,range=Association[]},
funcs = data[[All,All,i]];
range["sigma"] = {1.,2.5};
range["pi"] = {0.6,1.1};
range["u"] = {0.,0.7};
range["v"] = {-2,0.01};
range["u2"] = {0,0.08};

(*cardy = Association[Map[2+(#-2)4/Pi^2->SelectClosest[funcs[#],2+(#-2)4/Pi^2]&,Keys[funcs]]];*)
cardypoints = {{2,2},{2.05`,2.18325`},{2.1`,2.3649999999999998`},{2.15`,2.5549999999999997`},{2.2`,2.7375`},{2.25`,2.9625000000000004`},{2.3`,3.1625000000000005`},{2.4`,3.6125`},{2.5`,4.012499999999999`},{2.75`,4.8375`}};
cardypoints = Interpolation[Transpose[Reverse[Transpose[cardypoints]]]];
cardy = Association[Map[cardypoints[#]->SelectClosest[funcs[#],cardypoints[#]]&,Keys[funcs]]];
stars = ListPlot[cardy, PlotMarkers->{"\[FivePointedStar]",50}, PlotStyle->Black];
p = ListPlot[Values[funcs], PlotTheme->"Monochrome",
	PlotLegends->Placed[LineLegend[Map["N="<>ToString[PaddedForm[#,{3,2}]]&,Keys[funcs]], LegendFunction->"Frame", LegendLayout->"Row",
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 40, Black], LegendMarkerSize->{20,40}],{Center,Below}], 
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 40],
	AxesLabel->{Style["d",FontSize->40,Black,Bold],Style[title,FontSize->40,Black,Bold]}, 
	ImageSize->Full, PlotRange->{{1.97,3.03},range[label]}, PlotStyle->curves,
	Ticks->Automatic,GridLines->Automatic,PlotMarkers->markers];
Export[dir<>label<>".png",Show[p,stars]];
];
dir = "~/Dropbox/Paper_with_Andrzej_2020/fixed_points/";
If[!DirectoryQ[dir],CreateDirectory[dir]];
plot[1,U''[Subscript[\[Rho],0]],"u"];
plot[2,Subscript[Z,\[Sigma]][Subscript[\[Rho],0]],"sigma"];
plot[3,Subscript[Z,\[Pi]][Subscript[\[Rho],0]],"pi"];
plot[4,U'[0],"v"];
plot[5,U''[0],"u2"];
];


PlotKappa[KeySelect[nested,MemberQ[{4.,3.,2.5,2.25,2.},#]&]]


Keys[nested[2.5]]


2./Pi^2


PlotPotentials[nn_]:=Block[{len,ds,lim=4,\[Rho]Maxs,fixedpoints,vs,\[Pi]s,\[Sigma]s,plot,dir},
len = Length[nested[nn]]-6;
ds = Table[Keys[nested[nn]][[Floor[i*len/(lim+1)]+1]],{i,0,lim+1}];
ds = Drop[ds,{2}];
If[nn==2.5,ds = {3.`,2.5000000000000018`,2.200000000000003`,2.0500000000000025`}; lim=3;];
\[Rho]Maxs = Map[Keys[nested[nn,#]]&,ds];
fixedpoints = Map[nested[nn,#][[1]]&,ds];
vs = Table[Interpolation[Table[{i \[Rho]Maxs[[j]]/gridSize, v[i]/.fixedpoints[[j]]},{i,0,gridSize}]],{j,1,lim+1}];
\[Sigma]s = Table[Interpolation[Table[{i \[Rho]Maxs[[j]]/gridSize, z[i]/.fixedpoints[[j]]/.z[0]->1},{i,0,gridSize}]],{j,1,lim+1}];
\[Pi]s = Table[Interpolation[Table[{i \[Rho]Maxs[[j]]/gridSize, (z[i]-2 i \[Rho]Maxs[[j]]/gridSize yy[i])/.fixedpoints[[j]]/.z[0]->1},{i,0,gridSize}]],{j,1,lim+1}];
l1 = .025;
l2 = .001;
s = l1/2;
dash[1] = {};
dash[2] = {l1,s};
dash[3] = {l1,s,l2,s};
dash[4] = {l1, s, l2, s, l2, s};
dash[5] = {l2,s};
curves = Table[Directive[AbsoluteThickness[6],ColorData[97, "ColorList"][[i]]],{i,1,5}];
markers = Tuples[{{"\[FilledCircle]", "\[FilledSquare]", "\[FilledDiamond]", "\[FilledUpTriangle]",""},{30}}];
plot[fs_,title_,label_]:=Block[{p,funcs,range=Association[]},
funcs = MapThread[#2[t 0.845 #1]&,{\[Rho]Maxs,fs}];
range[Subscript[Z,\[Sigma]]] = {0.75,2.5};
range[Subscript[Z,\[Pi]]] = {0.6,1.1};
range["U'[\[Rho]]"] = {-2.05,2.05};

p = Plot[funcs,{t,0,1.15}, PlotTheme->"Monochrome",
	PlotLegends->Placed[LineLegend[Map["d="<>ToString[PaddedForm[#,{3,2}]]&,ds], LegendFunction->"Frame", LegendLayout->"Row",
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 40, Black], LegendMarkerSize->{20,40}],{Center,Below}], 
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 40],
	AxesLabel->{Style["\[Rho]",FontSize->40,Black,Bold],Style[title,FontSize->40,Black,Bold]}, 
	ImageSize->Full, PlotRange->{{0,1.15},range[title]}, PlotStyle->curves,
	Ticks->Automatic,GridLines->Automatic];
Export[dir<>label<>".png",p];
];
dir = "~/Dropbox/Paper_with_Andrzej_2020/fixed_points/"<>ToString[ToString[PaddedForm[nn,{3,2}]]]<>"/";
If[!DirectoryQ[dir],CreateDirectory[dir]];
plot[vs,"U'[\[Rho]]","v"];
plot[\[Sigma]s,Subscript[Z,\[Sigma]],"sigma"];
plot[\[Pi]s,Subscript[Z,\[Pi]],"pi"];
];


PlotPotentials[2.5]


Map[PlotPotentials,Keys[nested]];


integrands = integrandsListIsotropic;


AppendTo[integrands[[1]],integrands[[1,3]]/.\[Rho]->0];
AppendTo[integrands[[2]],integrands[[2,3]]/.\[Rho]->0];
AppendTo[integrands[[3]],integrands[[4,3]]];
integrands = Drop[integrands,{4}];


eqs = IntegrateEquations[integrands/.n->2.5];


sm3 = StabilityMatrix[integrands/.n->2.5, fp3,ConstRep[d3,\[Rho]Max3]];


sm2 = StabilityMatrix[integrands/.n->2.5, fp2,ConstRep[d2,\[Rho]Max2]];


PlotFixedPoint[fp2,\[Rho]Max2]


es3 = EigenSys[sm3];
es2 = EigenSys[sm2];


es2[[1,1]]


SerializeEigenvector[ev_,\[Rho]Max_]:=Block[{serialized=Association[]},
serialized[v] = Map[ev[[#]]&,Range[1,gridSize+1]];
serialized[zs] = Join[{0},Map[ev[[gridSize+1+#]]&,Range[1,gridSize]]];
serialized[zp] = Join[{0},Map[ev[[gridSize+1+#]]-2 \[Rho]Max #/gridSize ev[[#+2*gridSize+2]] &,Range[1,gridSize]]];
Return[serialized];
]


serial = SerializeEigenvector[Re[es2[[1,2]]],\[Rho]Max2];
serial2 = SerializeEigenvector[-Re[es3[[1,2]]],\[Rho]Max2];


ListPlot[{serial2[v],serial[v]}]


ListPlot[serial]


ListPlot[{-es3[[1,2]],es2[[1,2]]},PlotRange->All]


nested = KeySelect[nested, #>=2&];


high = nested[2.5][[1]];
mid = nested[2.5][[19]];
low = nested[2.5][[-1]];


{PlotFixedPoint[high],PlotFixedPoint[mid],PlotFixedPoint[low]}


MaxDer[fp_] := Block[{params = fp[[1]],\[Rho]Max=Keys[fp][[1]],\[Rho]Points,func,der},
\[Rho]Points= Range[0,gridSize]*\[Rho]Max/gridSize;
func = Map[(f[#]->z[#]+\[Rho]Points[[#+1]]*2*yy[#]&),Range[0,gridSize]]/.params/.ConstRep[d,\[Rho]Max];
der = 1/eps NumericDerivatives[Map[ f',Range[0,gridSize]]]/.func/.ConstRep[d,\[Rho]Max];
Return[Max[der]];
];


vders = Map[Map[(1/eps NumericDerivatives[v'[0]] /. #[[1]]/.ConstRep[d,Keys[#][[1]]]&) ,#]&,nested];
vdersSel = Map[Keys[Select[#,(#<0.019)&]][[1]]&,vders];
cardyv = Transpose[{Values[vdersSel],Keys[vdersSel]}];


ListPlot[Values[vders],PlotLegends->Keys[vders]]


zders = Map[Map[MaxDer,#]&,nested];
zdersSel =Map[Keys[Select[#,(#>0.16)&]][[1]]&,zders];
cardyz = Transpose[{Values[zdersSel],Keys[zdersSel]}];


ListPlot[Values[zders],PlotLegends->Keys[zders], PlotRange->{0,0.2}]


ListPlot[{zdersSel,vdersSel}]


Show[Plot[2+(x-2) Pi^2/4, {x, 2,2.4}],ListPlot[{cardyv,cardyz}, PlotMarkers->"OpenMarkers"],PlotRange->All]


Map[PlotFixedPoint[#[[1]],Keys[#][[1]]]&,{high,mid,low}]


low


V= Interpolation[Transpose[{Range[0,120]*50./120,Table[v[i],{i,0,120}]/.low[[1]]}]];


Vhigh= Interpolation[Transpose[{Range[0,120]*7.6/120,Table[v[i],{i,0,120}]/.high[[1]]}]];


U[x_] := NIntegrate[V[y],{y,0,x}]


Uh[x_] := NIntegrate[Vhigh[y],{y,0,x}]


Plot[Uh[x^2/2],{x,0,6}]


Plot[U[x^2/2],{x,0,10.2}]


{fps,exps} = Import["/home/andrzej/FunctionalRenormalization/ON_functional_smooth/regulator_dependence.wdx"];


nested = NestAssociation[exps];


Keys[nested[\[Nu],2.]]


range := Table[#[[-i]],{i,1,4}]&;
ListPlot[range[-nested[e2,2.]], PlotLegends->range[Keys[nested[e2,2.]]]]


range := Table[#[[-i]],{i,1,5}]&;
ListPlot[range[nested[e1,2.]], PlotLegends->range[Keys[nested[e1,2.]]]]


Keys[exps]


{fps,exps} = Import["/home/andrzej/FunctionalRenormalization/fixed_points1.wdx"];


packed = SelectClosest[NestAssociation[fps][1.],2.];
\[Rho]Max = Keys[packed][[1]];
fp=packed[[1]];


functions = {v,z,yy};
fp2 = Flatten[Map[Table[#[i]->(#[2*i]/.fp),{i,If[BooleanQ[#==z],1,0],gridSize}]&,functions]];
fp2 = Join[fp2,{\[Eta]->(\[Eta]/.fp)}];


AppendTo[integrandsListIsotropic[[1]],integrandsListIsotropic[[1,3]]/.\[Rho]->0];
AppendTo[integrandsListIsotropic[[2]],integrandsListIsotropic[[2,3]]/.\[Rho]->0];
AppendTo[integrandsListIsotropic[[3]],integrandsListIsotropic[[4,3]]];
integrandsListIsotropic = Drop[integrandsListIsotropic,{4}];


eqlist = IntegrateEquations[integrandsListIsotropic];


ntab = Join[Table[1+0.1 i,{i,0,5}],Table[1.5+0.05 i,{i,1,8}],Table[1.9+0.01 i,{i,1,10}]];


Length[fp2]


Length[eqlist]


eqlist/.n->1/.ConstRep[2.,\[Rho]Max, 2.]/.fp2


Get["FixedPoint.wl"]


fpdict = FindFixedPointDict[eqlist, fp2, ntab, {"N"}, <|"dim"->2.,"alpha"->2.,"\[Rho]Max"->\[Rho]Max|>];


exps = FindExponents[integrandsListIsotropic, {}, fpdict, {"N","\[Rho]Max"},  <|"dim"->2.,"alpha"->2.|>];


myExps = KeySelect[NestAssociation[UnnestAssociation[exps]],MemberQ[{e2,e3,e4,e5},#]&];


ListPlot[Im[myExps], AxesLabel->{"N"}]


Export["~/Desktop/eigenvalues.png", ListPlot[Re[myExps], PlotLabel->"Eigenvalues in d=2", AxesLabel->{"N"},ImageSize->Large]]


Export["~/Desktop/d=2O(N).wdx",myExps];


myExps = Import["~/Desktop/d=2O(N).wdx"];


{fps,allExponents} = Import[cacheDirectory<>REGULATORDEPENDENCEFILE];
FindPMS[exps_,e_,returnVals_]:=Block[{pms=Association[], scanned, SinglePMS, j, vals, keys, noPMS},
	scanned = TransposeAssociation[exps[e]];
	SinglePMS[assoc_]:= Block[{interpolation, pos, root, found=False},
		der[f_,x_] := (f[[2;;-1]]-f[[1;;-2]])/(x[[2;;-1]]-x[[1;;-2]]);
		deri = der[Keys[assoc],Values[assoc]];
		For[i=1,i<Length[assoc]-1,i++,
			If[Sign[deri[[i]]]!=Sign[deri[[i+1]]],
				found = True;
				pos = i+1;
				Break[];
			];
		];
		If[!found, Return[{}]];
		Return[{Keys[assoc][[pos]],assoc[[pos]]}];
	];
	If[returnVals, j=2,j=1];
	keys = Table[Keys[scanned][[i]],{i,1,Length[scanned]}];
	vals = Map[SinglePMS[KeySort[#]]&,Values[scanned]];
	noPMS = Position[vals,{}];
	keys = Delete[keys, noPMS];
	vals = Delete[vals, noPMS];
	Return[AssociationThread[keys, vals[[All,j]]]];
];


{fps,allExponents} = Import["/home/andrzej/Documents/Uczelnia/Anizotropie/FunctionalRenormalization/lowd_functional_smooth backup/fixed_points.wdx"];


nested = NestAssociation[allExponents];


nested


Keys[nested[e2,2.]]


Keys[nested[e2,2.]]


ListPlot[SelectClosest[SelectClosest[nested[e2],2.],2.]]



