(* ::Package:: *)

(* ::Subsubsection::Closed:: *)
(*Initialization*)


(* ::Text:: *)
(*Imports ModuleNaming module and sets the working directory*)


$HistoryLength = 10;
gridSize = 120;


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


(* ::Text:: *)
(*Serializes the command line arguments*)


If[debug,
	(* Debug run configuration *)
	task = "Z4"; (* ON/ Z4 / tricritical/ regulator /lowd*)
	configuration = 2;
	ON = False;
	regulatorName = "smooth";
	LPA = True,
	
	(* Release run configuration *)
	LPA = False;
	task = Null;
	anisotropySymmetry = 4;
	regulatorName = Null;
	configuration = Null;
	Block[{i},
	For[i = 1, i<=Length[arguments], i++,
		If[arguments[[i]] == "-task", task = arguments[[i+1]]];
		If[arguments[[i]] == "-regulator", regulatorName = arguments[[i+1]]];
		If[arguments[[i]] == "-conf", configuration = ToExpression[arguments[[i+1]]]];
		If[arguments[[i]] == "LPA", LPA = True];
		];
	];
];
ON = MemberQ[ONTASKS, task];

(* Check if all required options are specified *)
If[task == Null || regulatorName == Null || (anisotropySymmetry == "ON" && nn == Null),
	Print["Specify all required options: task, regulator"];
	Quit[];
	];


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
If[LPA,
	zNormalization = Floor[5*gridSize/6],
	zNormalization = 0];

Get["FixedPoint.wl"];
Get["Plots.wl"];

(* Select regulator function *)
SelectRegulator[regulatorName];
cacheDirectory = GetCacheDirectory[task, 0, regulatorName, LPA];


(* ::Subsubsection::Closed:: *)
(*Importing equations*)


(* ::Text:: *)
(*Imports cached flow equations*)


isotropicEquationsFile = GetEquationsFileName[task, 0, 0];
If[FileExistsQ[isotropicEquationsFile], 
	integrandsListIsotropic = Import[isotropicEquationsFile];
	If[LPA, integrandsListIsotropic= MakeLPA[integrandsListIsotropic];],
	Print["No isotropic equations found"]; Quit[]];
If[ON,
integrandsListAnisotropic = {},
anisotropicEquationsFile = GetEquationsFileName[task, 0, 1];
If[FileExistsQ[anisotropicEquationsFile], 
	integrandsListAnisotropic = Import[anisotropicEquationsFile],
	Print["No anisotropic equations found"]; Quit[]];
If[LPA, integrandsListAnisotropic= MakeLPA[integrandsListAnisotropic];]
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
		eqlistIntermediate = IntegrateEquations[integrandsListIsotropic /. n -> 4.];
		{guess, \[Rho]Max} = \[Rho]MaxScan[eqlistIntermediate, guess, dim, \[Rho]Max, 2., 0.05*\[Rho]Max, 40];
	];
	integrandsListIsotropic = (integrandsListIsotropic /. n -> constants["N"]);
	integrandsListAnisotropic = (integrandsListAnisotropic /. n -> constants["N"])];

eqlist = IntegrateEquations[integrandsListIsotropic];

scan = \[Rho]MaxScan[eqlist, guess, dim, \[Rho]Max, 2., 0.05*\[Rho]Max, 40];
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
(*Z4-anisotropy dimension dependence*)


If[ task == "Z4",
{guess, dim, \[Rho]Max} = Import[guessfile];
guess = AdjustNumberOfPoints[guess, gridSize,zNormalization];
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
	If[alpha>2, FindInitialGuess[2., alpha, +1.],
		scan = \[Rho]MaxScan[eqlist, guess, dim, \[Rho]Max, alpha, 0.05*\[Rho]Max, 40];
		CheckScan[scan];
		];
	];

fixedPointDict = FixedPointDimScan[eqlist, guess, constants, dim, \[Rho]Max];

exponents = FindExponents[integrandsListIsotropic, integrandsListAnisotropic, fixedPointDict, {"dim","\[Rho]Max"}, constants];

fixedPointDict = UnnestAssociation[TransposeAssociation[Association[alpha->fixedPointDict]]];
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


(* ::Subsubsection:: *)
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



(* ::Subsubsection::Closed:: *)
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
guess = AdjustNumberOfPoints[guess, gridSize,zNormalization];

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


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
(*Plotting critical exponents dimension dependence *)


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


If[ task == "ON" || task == "tricritical",
{fps,allExponents} = Import["ON_functional_smooth backup/"<>FIXEDPOINTSFILE];
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

data = Map[KeySelect[#,(#>=2)&]&, KeySort[allExponents[exp]]];

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
Map[PlotDimDependence, Keys[allExponents]];
];


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
(*Other*)


If[!debug, Quit[]];
(* Anything below this line won't be executed in command line run *)


{fps,exps} = Import[cacheDirectory<>FIXEDPOINTSFILES[3]];
nestedfps = NestAssociation[fps];
nestedexps = NestAssociation[exps];


{guess, dim, \[Rho]Max} = Import[guessfile];
guess = AdjustNumberOfPoints[guess, gridSize,zNormalization];

configurations = Import[runconf];

If[!Element[configuration,Integers] || configuration<1 || configuration > Length[configurations],
	Print["Choose proper run cofiguration id"];
	Quit[]];
constants = configurations[[configuration]];

eqlist = IntegrateEquations[integrandsListIsotropic];


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
