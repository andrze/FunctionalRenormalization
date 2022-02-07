(* ::Package:: *)

(* ::Chapter:: *)
(*Initialization*)


(* ::Subsection::Closed:: *)
(*Setting the project path and importing the naming module*)


(* ::Text:: *)
(*Imports ModuleNaming module and sets the working directory*)


$HistoryLength = 10;


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
	task = "Z4"; (* ON/ Z4 / tricritical/ regulator /lowd / ising / ising2 *)
	configuration = 24;
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


(* Set desired gridsize for potential minimum *)
gridSizeDict = <|"ON" -> 200, "Z4" -> 100, "regulator" -> 120, "lowd" -> 120, "ising" -> 300, "ising2" -> 300, "tricritical" -> 100|>;
gridSize = Lookup[gridSizeDict,task,120];

If[debug,
	gridSize = 30];


(* Set desired position for potential minimum *)
potentialMinimumDict = <|"ON" -> Floor[85/100*gridSize], "Z4" -> Floor[1/2*gridSize], 
	"regulator" -> Floor[50/60*gridSize], "lowd" -> Floor[55/60*gridSize], 
	"ising" -> Floor[7/20*gridSize], "ising2" -> Floor[7/20*gridSize], "tricritical" -> Floor[gridSize/2]|>;
potentialMinimum = Lookup[potentialMinimumDict,task,Floor[1/2*gridSize]];

(* Set z-normalization point *)

If[LPA,
	zNormalization = Floor[5*gridSize/6],
	zNormalization = 0];

largeNormalizationTasks = {"ising", "ising2"(*,"Z4","ON"*)};
If[MemberQ[largeNormalizationTasks,task],
	zNormalization = potentialMinimum-Floor[gridSize/20]];

Print["Grid size: "<> ToString[gridSize]<>" Potential minimum: "<>
	ToString[potentialMinimum]<>" Z normalization point: "<> ToString[zNormalization]];


(* ::Subsection::Closed:: *)
(*Choosing the configurations files and importing modules*)


(* Select run configuration file *)
runconf = Null;
If[task=="ON", runconf = ONRUNCONF];
If[task=="Z4", runconf = Z4RUNCONF];
If[task=="tricritical", runconf = TRICRITICALRUNCONF];
If[task=="ising", runconf = ISINGRUNCONF];
If[task=="ising2", runconf = ISING2RUNCONF];
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


(* ::Chapter:: *)
(*O(N) models*)


(* ::Subsection:: *)
(*2D Ising*)


(* ::Subsubsection::Closed:: *)
(*Ising*)


If[task=="ising" && ! debug,
	integrands = (integrandsListIsotropic/.n->1)[[1;;2]];
	If[zNormalization>0,
		integrands = (integrands /. zp[0] ->zs[0])];
	
	eqlist = IntegrateEquations[integrands];

	{guess, guessDim, \[Rho]Max} = Import[guessfile];
	guess = AdjustNumberOfPoints[guess, gridSize,zNormalization,\[Rho]Max];
	guess = Normal[KeySelect[Association[guess],(MatchQ[#,f_[i_]/;MemberQ[{v,zs},f]]||MatchQ[#,\[Eta]])&]];

	{newGuess, \[Rho]Max2} = \[Rho]MaxScan[eqlist, guess, guessDim, \[Rho]Max, 2, \[Rho]Max/20, 50];
	
	dict = FindFixedPointDict[eqlist, newGuess, 3-Range[10] 1/10, "dim", <|"N"->1,"alpha"->2,"\[Rho]Max"->\[Rho]Max2|>];
	fp = dict[[-1]];
	\[Rho]m = Keys[dict][[-1,2]];
	
	lowAlpha = FindFixedPointDict[eqlist, fp, 2-Range[0,16]/10, "alpha", <|"N"->1,"\[Rho]Max"->\[Rho]m, "dim"->2|>,True];
	highAlpha = FindFixedPointDict[eqlist, fp, 2+ Range[35]/5, "alpha", <|"N"->1,"\[Rho]Max"->\[Rho]m, "dim"->2|>,True];
	fixedPointDict = Join[lowAlpha,highAlpha];
	
	exps = FindExponents[integrands, {}, fixedPointDict, {"alpha","\[Rho]Max"}, <|"N"->1, "dim"->2|>];
	
	fixedPointDict = UnnestAssociation[fixedPointDict];
	exponents = UnnestAssociation[exps];
	
	AppendToCache[cacheDirectory, FIXEDPOINTSFILE, {fixedPointDict, exponents}]
]


(* ::Subsubsection::Closed:: *)
(*Ising2*)


If[task=="ising2" && ! debug,
	integrands = (integrandsListIsotropic/.n->1)[[1;;2]];
	If[zNormalization>0,
		integrands = (integrands /. zp[0] ->zs[0])];
	
	eqlist = IntegrateEquations[integrands];

	{guess, guessDim, \[Rho]Max} = Import[guessfile];
	guess = AdjustNumberOfPoints[guess, gridSize,zNormalization,\[Rho]Max];
	guess = Normal[KeySelect[Association[guess],(MatchQ[#,f_[i_]/;MemberQ[{v,zs},f]]||MatchQ[#,\[Eta]])&]];
	
	configurations = Import[runconf];
	startDim = 3;
	endDim = configurations[[configuration]]["dim"];
	steps = startDim + Sign[(endDim-startDim)]/10 Range[0,Abs[(endDim-startDim)]*10];
	steps = DeleteDuplicates[Append[steps, endDim]];
	
	dict = FindFixedPointDict[eqlist, guess, steps, "dim", <|"N"->1,"alpha"->2,"\[Rho]Max"->\[Rho]Max|>];
	fp = dict[[-1]];
	\[Rho]m = Keys[dict][[-1,2]];
	
	lowAlpha = FindFixedPointDict[eqlist, fp, 2-Range[0,16]/10, "alpha", <|"N"->1, "\[Rho]Max"->\[Rho]m, "dim"->endDim|>, True];
	highAlpha = FindFixedPointDict[eqlist, fp, 2+ Range[35]/5, "alpha", <|"N"->1, "\[Rho]Max"->\[Rho]m, "dim"->endDim|>, True];
	fixedPointDict = Join[lowAlpha, highAlpha];
	
	exps = FindExponents[integrands, {}, fixedPointDict, {"alpha","\[Rho]Max"}, <|"N"->1, "dim" -> endDim|>];
	
	fixedPointDict = UnnestAssociation[Association[endDim->fixedPointDict]];
	exponents = UnnestAssociation[Association[endDim->exps]];
	
	AppendToCache[cacheDirectory, FIXEDPOINTSFILE, {fixedPointDict, exponents}]
]


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
	If[constants["N"] > 2,
		eqlistIntermediate = IntegrateEquations[integrandsListIsotropic];
		nsteps = 2 + Range[0, 4*(constants["N"] - 2)]/4;
		dict = FindFixedPointDict[eqlistIntermediate, guess, nsteps, "N", <|"dim"->3,"alpha"->2,"\[Rho]Max"->\[Rho]Max|>];
		
		guess = dict[[-1]];
		\[Rho]Max = Keys[dict][[-1,2]];
	];
	integrands = (integrandsListIsotropic /. n -> constants["N"]);
	If[constants["N"]==1, 
		integrands = integrands[[1;;2]];
		guess = Join[guess[[1;;2*gridSize+1]],guess[[-1;;-1]]]],
	Print["Value of N for the calculation not specified"];
	Quit[];
	];
	
ChangeVariable[integrands_] := Block[{newIntegrands = integrands, VarChange},
	VarChange[i_] := (4 t^3 i /.{y->t^4})/.t->y;
	newIntegrands[[3]] = VarChange[newIntegrands[[3]]];
	If[Length[integrands] == 4,
		newIntegrands[[4]] = VarChange[newIntegrands[[4]]];
	];
	Return[newIntegrands];
];
integrands = Map[ChangeVariable, integrands];

Unprotect[IntegralConfigurations];
IntegralConfigurations = {{0,1/16,30},{1/16,9/10,40},{9/10,3/2,30}};
Protect[IntegralConfigurations];

eqlist = IntegrateEquations[integrands];
params = Flatten[Table[integrands[[i,1]]/.\[Rho]->eps j,{i,1,Length[integrands]}, {j, 4-Length[integrands[[i]]],gridSize}]];
params = DeleteCases[params, zs[zNormalization]];
AppendTo[params,\[Eta]];
jacobian = Grad[eqlist,params];

If[constants["alpha"]!=2,
	initialAlphas = 2+Range[0,constants["alpha"]-2];
	initialAlphas = DeleteDuplicates[Append[initialAlphas,constants["alpha"]]];
	alphaScan = FindFixedPointDict[eqlist, guess, initialAlphas, "alpha", <|"dim"->dim, "N"->constants["N"], "\[Rho]Max"->\[Rho]Max|>];
	If[Length[alphaScan]==0, Print["Initial alpha scan failed"]; Quit[]];
	guess = alphaScan[[-1]];
	\[Rho]Max = Keys[alphaScan][[-1,2]];];

(*{fixedPointDict,exponents} = FullDScan[eqlist, newGuess, dim, 1/20, \[Rho]Max, constants["alpha"]];*)
fixedPointDict = FixedPointDimScan[eqlist, jacobian, guess, constants, dim, \[Rho]Max];
If[Length[fixedPointDict]==0, 
	Print["No fixed points found"];
	Quit[],
	Print["Found "<>ToString[Length[fixedPointDict]]<>" fixed points."]];
exponents = FindExponents[integrands, {}, fixedPointDict, {"dim","\[Rho]Max"}, constants];

flatFixedPointDict = UnnestAssociation[<|constants["alpha"]-><|constants["N"]->fixedPointDict|>|>];

nestedExps = TransposeAssociation[<|constants["alpha"]->TransposeAssociation[<|constants["N"]->exponents|>]|>];
flatExps = UnnestAssociation[nestedExps];

AppendToCache[cacheDirectory, FIXEDPOINTSFILE, {flatFixedPointDict,flatExps}];
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


If[False && task == "ON",
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


(*{fps,allExponents} = Import["ON_functional_smooth backup2/"<>FIXEDPOINTSFILE];
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
curves = Table[Directive[AbsoluteThickness[6],ColorData[97, "ColorList"][[i]]],{i,1,5}];
markers = Tuples[{{"\[FilledCircle]", "\[FilledSquare]", "\[FilledDiamond]", "\[FilledUpTriangle]",""},{40}}];

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
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 50,Black],LegendMarkerSize->{40,40}],Bottom], LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 40],
	AxesLabel->{Style["d",FontSize\[Rule]50,Black,Bold],Style[label,FontSize\[Rule]50,Black,Bold]}, 
	ImageSize->Full, PlotRange->{{2,3},{min,max}},PlotMarkers->markers,
	PlotStyle->curves,Ticks->ticks[exp],GridLines->ticks[exp]];
Export["~/Desktop/"<>ToString[exp]<>".png",p];];
Map[PlotDimDependence, Keys[allExponents]];*)


If[ False && task == "ON",
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


(* ::Chapter:: *)
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
finalDimension = constants["dim"];
dimSteps = 3-Range[0,20*(3-finalDimension)]/20;
If[!MemberQ[dimSteps,finalDimension],
	AppendTo[dimSteps, finalDimension]];

eqlist = IntegrateEquations[integrandsListIsotropic];

dict = FindFixedPointDict[eqlist, guess, dimSteps, "dim", <|"N"->2,"alpha"->2,"\[Rho]Max"->\[Rho]Max|>];
fp = dict[[-1]];
\[Rho]m = Keys[dict][[-1,2]];

finalAlpha = constants["alpha"];
If[finalAlpha>2,
	aStep = 2/10;
	alphaSteps = 2+Range[1,(finalAlpha-2)/aStep]*aStep,
	
	aStep = 1/40;
	alphaSteps = 2-Range[1,(2-finalAlpha)/aStep]*aStep];
alphaSteps = DeleteDuplicates[Append[alphaSteps,finalAlpha]];

fixedPointDict = FindFixedPointDict[eqlist, fp, alphaSteps, "alpha", <|"N"->2,"\[Rho]Max"->\[Rho]m, "dim"->finalDimension|>];
If[Length[fixedPointDict] < Length[alphaSteps],
	Print["Alpha scan failed"];
	Quit[]];
fixedPointDict = fixedPointDict[[-1;;-1]];

exponents = FindExponents[integrandsListIsotropic, integrandsListAnisotropic, fixedPointDict, {"alpha","\[Rho]Max"}, <|"N"->2, "dim"->finalDimension|>];

fixedPointDict = UnnestAssociation[Association[finalDimension->fixedPointDict]];
exponents = UnnestAssociation[TransposeAssociation[Association[finalDimension->exponents]]];

AppendToCache[cacheDirectory,FIXEDPOINTSFILE,{fixedPointDict,exponents}];
];


If[ task == "Z4" && False,
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


If[False&& task == "Z4",
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

curves = Table[Directive[AbsoluteThickness[8],ColorData[97, "ColorList"][[i]]],{i,1,5}];
markers = Tuples[{{"\[FilledCircle]", "\[FilledSquare]", "\[FilledDiamond]", "\[FilledUpTriangle]",""},{40}}];

data = Association[Map[#->Quiet[FindPMS[allExponents,#,False]]&,exps[[1;;2]]]];
data2 = Association[Map[#->Quiet[FindPMS[allExponents,#,True]]&,exps]];
data3 = allExponents2[\[Eta],2.];
BigPlot[data_,legends_,axesLabels_,plotRange_:Automatic,ticks_:Automatic, PlotFunc_:ListPlot]:=Block[{p},
p = PlotFunc[data,PlotTheme->"Monochrome", PlotLegends->Placed[PointLegend[legends, LegendFunction->"Frame", LegendLayout->"Row",
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 50, Black],LegendMarkerSize->{50,50}],Bottom], 
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 50, Black],
	AxesLabel->{Style[axesLabels[[1]],FontSize->50,Black,Bold],Style[axesLabels[[2]],FontSize->50,Black,Bold]}, 
	ImageSize->Full, PlotRange->plotRange,PlotMarkers->markers,PlotStyle->curves,
	Ticks->ticks,GridLines->ticks];
	Return[p];
];
plotDir = "~/Desktop/";
Export[plotDir<>"PMS.png",
	BigPlot[Values[data],Map[("PMS("<>ToString[#]<>")")&,Keys[data]],{"d",Subscript[\[Alpha],PMS]}, {{1.975,3.025},{1.4,2.25}},Automatic,ListLinePlot]];
Export[plotDir<>"etaPMS.png",
	BigPlot[data2[\[Eta]],{},{"d",Subscript[\[Eta],PMS]},
		{{2.,3.01},{-0.005,0.305}},{{2.,2.2,2.4,2.6,2.8,3.},{0.,0.1,0.2,0.3}},ListLinePlot]];
Export[plotDir<>"etaPMS2.png",
	Show[BigPlot[{data2[\[Eta]]},{Subscript[\[Eta],PMS]},{"d",\[Eta]},
		{{2.,3.01},{-0.005,0.305}},{{2.,2.2,2.4,2.6,2.8,3.},{0.,0.1,0.2,0.3}},ListLinePlot],
		ListPlot[data3,PlotLegends->Placed[PointLegend[{"\[Eta](\[Alpha]=2)"}, LegendFunction->"Frame", LegendLayout->"Row",
			LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 50, Black], LegendMarkerSize->{40,20}],Bottom],
			PlotMarkers->{"\[FilledSquare]",56},PlotStyle-> ColorData[97, "ColorList"][[2]]]
		]];
Export[plotDir<>"nuPMS.png",
	BigPlot[data2[\[Nu]],{},{"d",Subscript[\[Nu],PMS]}]];
Export[plotDir<>"e1PMS.png",
	BigPlot[data2[e1],{},{"d",1/\[Nu]}]];
Export[plotDir<>"nu1.png",
	BigPlot[Values[allExponents[e1][[1,-5;;-3]]],Keys[allExponents[e1][[1]]][[-5;;-3]],{"\[Alpha]",1/\[Nu]}]];
Export[plotDir<>"eta1.png",
	BigPlot[Values[allExponents[\[Eta]][[1,-3;;-1]]],Map["d="<>ToString[#]&,Keys[allExponents[\[Eta]][[1]]][[-3;;-1]]],{"\[Alpha]",\[Eta]},
		{{0.95,3.05},{0.25,0.322}},{{1.0,1.5,2.0,2.5,3.0},{0.26,0.28,0.30,0.32}}]];

];


(*cacheDirectory = "ON_functional_smooth backup/";
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
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 50, Black],LegendMarkerSize->{40,20}],Bottom], 
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 50, Black],
	AxesLabel->{Style[axesLabels[[1]],FontSize\[Rule]50,Black,Bold],Style[axesLabels[[2]],FontSize\[Rule]50,Black,Bold]}, 
	ImageSize->Full, PlotRange->plotRange,PlotMarkers->{Automatic, 20},PlotStyle->AbsoluteThickness[4],
	Ticks->ticks,GridLines->ticks];
	Return[p];
];
cacheDirectory = "~/Desktop/";
Export[cacheDirectory<>"PMS.png",
	BigPlot[Values[data],Map[("PMS("<>ToString[#]<>")")&,Keys[data]],{"d",Subscript[\[Alpha],PMS]}, {{1.975,3.025},{1.4,2.25}}]];
Export[cacheDirectory<>"etaPMS.png",
	BigPlot[data2[\[Eta]],{},{"d",Subscript[\[Eta],PMS]},
		{{2.,3.01},{-0.005,0.305}},{{2.,2.2,2.4,2.6,2.8,3.},{0.,0.1,0.2,0.3}},ListLinePlot]];
Export[cacheDirectory<>"etaPMS2.png",
	Show[BigPlot[{data2[\[Eta]]},{Subscript[\[Eta],PMS]},{"d",\[Eta]},
		{{2.,3.01},{-0.005,0.305}},{{2.,2.2,2.4,2.6,2.8,3.},{0.,0.1,0.2,0.3}},ListLinePlot],
		ListPlot[data3,PlotLegends->Placed[PointLegend[{"\[Eta](\[Alpha]=2)"}, LegendFunction->"Frame", LegendLayout->"Row",
			LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 50, Black], LegendMarkerSize->{40,20}],Bottom],
			PlotMarkers->{"\[FilledSquare]",46},PlotStyle-> ColorData[97, "ColorList"][[2]]]
		]];
Export[cacheDirectory<>"nuPMS.png",
	BigPlot[data2[\[Nu]],{},{"d",Subscript[\[Nu],PMS]}]];
Export[cacheDirectory<>"e1PMS.png",
	BigPlot[data2[e1],{},{"d",1/\[Nu]}]];
Export[cacheDirectory<>"nu1.png",
	BigPlot[Values[allExponents[e1][[1,-5;;-3]]],Keys[allExponents[e1][[1]]][[-5;;-3]],{"\[Alpha]",1/\[Nu]}]];
Export[cacheDirectory<>"eta1.png",
	BigPlot[Values[allExponents[\[Eta]][[1,-3;;-1]]],Map["d="<>ToString[#]&,Keys[allExponents[\[Eta]][[1]]][[-3;;-1]]],{"\[Alpha]",\[Eta]},
		{{0.95,3.05},{0.25,0.322}},{{1.0,1.5,2.0,2.5,3.0},{0.26,0.28,0.30,0.32}}]];*)


(* ::Chapter:: *)
(*Other*)


If[!debug, Quit[]];
(* Anything below this line won't be executed in command line run *)


(* ::Section:: *)
(*Z4-symmetric anisotropy*)


(* ::Subsection::Closed:: *)
(*Check for grid size convergence*)


(*{fps,allExponents} = Import["Z4_functional_smooth 100/"<>FIXEDPOINTSFILE];
{fps2,allExponents2} = Import["Z4_functional_smooth 200/"<>FIXEDPOINTSFILE];
selected1 = KeySort[KeySelect[allExponents, MemberQ[Keys[allExponents2],#]&]];
selected2 = KeySort[KeySelect[allExponents2, MemberQ[Keys[allExponents],#]&]];
Max[Abs[selected1-selected2]/Abs[selected1]]*)


(* ::Subsection:: *)
(*Checks for differences between regulators*)


(* ::Subsubsection::Closed:: *)
(*Initialization*)


schemes = {"smooth", "smooth_strict","exponential", "exponential_strict","litim", "litim_strict"};
allExponents = <| |>;


ReadData[scheme_] := Block[{fps, exps, nested},
{fps,exps} = Import["Z4_functional_"<>scheme<>"/"<>FIXEDPOINTSFILE];
allExponents[scheme] = NestAssociation[exps];
Return[Null]];


Map[ReadData,schemes];


FindYReal[scheme_]:=Block[{newData = allExponents[scheme],isotropic1,isotropic2, anisotropic, yrealdata, threshold=10^-5},
	isotropic1 = UnnestAssociation[newData[e1]];
	isotropic2 = UnnestAssociation[newData[e2]];
	anisotropic = Map[UnnestAssociation[newData[#]]&,{y1,y2,y3,y4,y5}];
	
	yrealdata = <| |>;
	
	GetYReal[key_]:=Block[{diff1,diff2,drops1,drops2, positions},
		diff1 = Map[#[key]&,anisotropic] - isotropic1[key];
		diff2 = Map[#[key]&,anisotropic] - isotropic2[key];
		drops1 = Map[Abs[Im[#]]>0 || Abs[#]<threshold&,diff1];
		drops2 = Map[Abs[Im[#]]>0 && Abs[#]<threshold&,diff2];
		positions = Position[Map[drops1[[#]]||drops2[[#]]&,Range[Length[drops1]]],False];
		If[Length[positions]>0,
			yrealdata[key] = anisotropic[[positions[[1,1]]]][key];];
		Return[Null];];
	Map[GetYReal,Keys[isotropic1]];
	newData[yRe] = NestAssociation[yrealdata];
	allExponents[scheme] = newData;
Return[Null]];


FindYComplex[scheme_]:=Block[{newData = allExponents[scheme],isotropic1,isotropic2, anisotropic, yrealdata, threshold=10^-5},
	isotropic1 = UnnestAssociation[newData[e1]];
	isotropic2 = UnnestAssociation[newData[e2]];
	anisotropic = Map[UnnestAssociation[newData[#]]&,{y1,y2,y3,y4,y5}];
	
	yrealdata = <| |>;
	
	GetYReal[key_]:=Block[{diff1,diff2,drops1,drops2, positions},
		diff1 = Map[#[key]&,anisotropic] - isotropic1[key];
		diff2 = Map[#[key]&,anisotropic] - isotropic2[key];
		drops1 = Map[Abs[Im[#]]==0 || Abs[#]<threshold&,diff1];
		drops2 = Map[Abs[Im[#]]==0 && Abs[#]<threshold&,diff2];
		positions = Position[Map[drops1[[#]]||drops2[[#]]&,Range[Length[drops1]]],False];
		If[Length[positions]>0,
			yrealdata[key] = Re[anisotropic[[positions[[1,1]]]][key]]];
		Return[Null];];
	Map[GetYReal,Keys[isotropic1]];
	newData[yComplex] = NestAssociation[yrealdata];
	allExponents[scheme] = newData;
Return[Null]];


Map[FindYReal,schemes];


Map[FindYComplex,schemes];


GetLegend[scheme_]:=Block[{label,regulator,underscore,labelDict},
labelDict = <|"smooth" -> "Wetterich", "exponential" -> "Exponential", "litim" -> "Litim2"|>;
underscore = StringPosition[scheme,"_"];
If[Length[underscore]>0,
	regulator = StringTake[scheme,1;;underscore[[1,1]]-1];,
	regulator = scheme;];
label = labelDict[regulator];
If[Length[StringPosition[scheme,"_"]]>0,
	label = label<>" strict",
	label = label<>" anzatz"];
Return[label]];


(* ::Subsubsection::Closed:: *)
(*All regulators in one plot*)


PMS[dat_]:= Block[{norm,inter,pms},
norm = Map[{#[[1]],#[[2]]}&,Normal[dat]];
inter = Interpolation[norm];
pms = \[Alpha]/.FindRoot[inter'[\[Alpha]],{\[Alpha], Max[1, Min[norm[[All,1]]]+0.2], Min[norm[[All,1]]], Max[norm[[All,1]]]}, DampingFactor->2];
Return[{pms,inter[pms]}]];


PMSmin[dat_]:= Block[{norm,inter,pms},
norm = Map[{#[[1]],#[[2]]}&,Normal[dat]];
inter = Interpolation[norm];
pms = FindMinimum[inter[\[Alpha]],{\[Alpha], Max[1, Min[norm[[All,1]]]+0.2], Min[norm[[All,1]]], Max[norm[[All,1]]]}];
Return[pms[[1]]]];


PMSmax[dat_]:= Block[{norm,inter,pms},
norm = Map[{#[[1]],#[[2]]}&,Normal[dat]];
inter = Interpolation[norm];
pms = FindMaximum[inter[\[Alpha]],{\[Alpha], Max[1, Min[norm[[All,1]]]+0.2], Min[norm[[All,1]]], Max[norm[[All,1]]]}];
Return[pms[[1]]]];


PlotRegulatorComparison[dim_,exp_,customTitle_:"", customRange_:{{0,0},{0,0}}]:= Block[{title=customTitle, data, pms, p, norm, 
	pairs, inter, pms2, range=Automatic, axisLabel=exp, plotLabel, pmsvals},

data = TransposeAssociation[TransposeAssociation[allExponents][exp]][dim];

If[MemberQ[{e1,\[Eta]}, exp],
	pmsvals = Map[PMSmin,data],
	pmsvals = Map[PMSmax,data]];
plotLabel = "PMS values\n" <> Map[ToString[GetLegend[Keys[pmsvals][[#]]]]<>": "<>ToString[pmsvals[[#]]]<>"\n"&,Range[Length[pmsvals]]];

If[customRange != {{0,0},{0,0}},
	range = customRange,
	If[exp == yRe,
		range = {{0,5},{-1.1,0.1}}];
	If[exp == yComplex,
		range = {{0,5},{-1.2,0.2}}]];

If[TrueQ[exp==yComplex] && TrueQ[exp==yRe],
	axisLabel = Subscript[y,2];
];


p = ListPlot[Values[data], AxesLabel->{"\[Alpha]",axisLabel},ImageSize->Large, PlotLegends->PointLegend[Map[GetLegend,Keys[data]]],
	Ticks->Automatic,GridLines->Automatic, PlotRange->range, PlotLabel->plotLabel];

If[customTitle == "",
	title = ToString[exp]<>"(d="<>ToString[NumberForm[N[dim],{3,2}]]<>").png"];
Export["~/Desktop/Comparisons/"<>title,p];
Return[p];
]



Quiet[Map[PlotRegulatorComparison[#[[1]],#[[2]]]&, Tuples[{{3,275/100,25/10,24/10,23/10,22/10,210/100,205/100,201/100,2},{e1,\[Eta]}}]]];


Quiet[Map[PlotRegulatorComparison[#,yRe]&,Keys[allExponents["smooth",yRe]]]];
Quiet[Map[PlotRegulatorComparison[#,yComplex]&,Select[Keys[allExponents["smooth",yComplex]],#<2.35&]]];


Quiet[PlotRegulatorComparison[3, yRe, "Zoom In(d=3).png",{{0,5.1},{-0.125,-0.141}}]];
Quiet[PlotRegulatorComparison[2, yRe, "Zoom In(d=2).png",{{0,5.1},{0.1,-0.55}}]];


(* ::Subsubsection::Closed:: *)
(*Separate plots for regulators*)


PlotSelected[scheme_, exps_, dims_, title_, range_:{}]:= Block[{expons={}, reals, complexes,SelectStyle, p, p1, p2, legend, 
	yrange, min, max, curves, markersReal, markersComplex, RealQ, ComplexQ},
SelectStyle[length_]:=Block[{},
	curves = Join[
	Table[Directive[AbsoluteThickness[0.01],ColorData[97, "ColorList"][[i]]],{i,1,length}],
	Table[Directive[AbsoluteThickness[2],ColorData[97, "ColorList"][[i]]],{i,1,length}]];
	markersReal = Join[Tuples[{Table["\[FilledCircle]",{i,1,length}],{12}}],Tuples[{Table["\[Cross]",{i,1,length}],{1}}]];
	markersComplex = Join[Tuples[{Table["\[FivePointedStar]",{i,1,length}],{16}}],Tuples[{Table["\[Cross]",{i,1,length}],{1}}]];];

expons = Map[Values[KeySelect[allExponents[scheme,#],(Min[Abs[#-dims]]<0.001)&]]&,exps];

If[Length[exps]==2,
	SelectStyle[Length[dims]];
	legend = Map["d="<>ToString[N[#]]&,Sort[dims,#1>#2&]],
	If[Length[dims]<=2,
		expons = Transpose[expons];
		SelectStyle[Length[exps]];
		legend = Map[ToString[#]&,exps];,
		
		Print["Either exps or dims should have length no larger than 2"];
		Return[];
	];
];
expons = Flatten[expons];

If[Length[expons]==0,Return[]];

If[Length[range]==0,
min = Min[Re[expons]];
max = Max[Re[expons]];
yrange = {Min[min-(max-min)*0.05,-0.025],Max[max+(max-min)*0.05,0.025]},

yrange=range;
];

RealQ[x_] := Im[x]==0;
ComplexQ[x_] := Im[x]!=0;

reals = Map[Select[#,RealQ]&,expons];
complexes = Re[Map[Select[#,ComplexQ]&,expons]];
Map[If[Length[reals[[#]]]==0, reals[[#]] = <|-1->Mean[yrange]|>;]&,Range[Length[reals]]];
Map[If[Length[complexes[[#]]]==0, complexes[[#]] = <|-1->Mean[yrange]|>;]&,Range[Length[complexes]]];
p1 = ListLinePlot[reals, PlotStyle->curves, PlotMarkers->markersReal, PlotLegends->PointLegend[legend,LegendMarkerSize->{30,20}, 
	LabelStyle->Directive[Bold, 18]],PlotRange->{{0,5.1},yrange}];
p2 = ListLinePlot[complexes, PlotStyle->curves, PlotMarkers->markersComplex];
p = Show[p1,p2, LabelStyle->Directive[Bold, 18], 
	 PlotLabel->Style[ToString[title],FontSize->20],
	ImageSize->Large, AxesLabel->{"\[Alpha]"},PlotRange->{{0,5.1},yrange},GridLines->Automatic];
	
Export["~/Desktop/"<>GetLegend[scheme]<>"/"<>ToString[title]<>".png",p];
];


Keys[allExponents]


Map[CreateDirectory["~/Desktop/"<>ToString[#]]&,Keys[allExponents]] 


dims = {3,11/4,5/2,12/5,23/10,11/5,21/10,41/20,201/100,2};


dims = Keys[allExponents["smooth",y1]];


schemes = Keys[allExponents];


Map[PlotSelected[#,{y1,e1}, dims, "First EV"]&, schemes];


Map[PlotSelected[#, {y1,y2,y3,y4,y5,e1,e2},{2.},"Zoom In (d=2)",{1.025,-1.025}]&, schemes];


PlotDims[scheme_]:=Map[PlotSelected[scheme, {y1,y2,y3,y4,y5,y6,y7,y8,e1,e2},{#},"EVs (d="<>ToString[PaddedForm[N[#],{4,3}]]<>")",{2.025,-4}]&,dims];


PlotDims["smooth"]


Map[PlotDims, schemes];


Map[PlotSelected[#, {\[Eta],\[Eta]},{5/2,23/10,11/5,21/10,41/20,201/100,2},"eta"]&, schemes];


(* ::Subsubsection::Closed:: *)
(*Other*)


{fps,allExponents} = Import["Z4_functional_smooth/"<>FIXEDPOINTSFILE];


PMS[dat_]:= Block[{norm,inter,pms},
norm = Map[{#[[1]],#[[2]]}&,Normal[dat]];
inter = Interpolation[norm];
pms = \[Alpha]/.FindRoot[inter'[\[Alpha]],{\[Alpha],2}];
Return[{pms,inter[pms]}]];


PMS2[dat_]:= Block[{norm,inter,pms},
norm = Map[{#[[1]],#[[2]]}&,Normal[dat]];
Return[Max[norm[[All,2]]]]];


nested = NestAssociation[allExponents];


yData = allExponents["smooth"][yRe];
yData2 = allExponents["exponential"][yRe];


pmss = Map[PMS2,yData];
pmss2 = Map[PMS2,yData2];


p = ListPlot[{pmss,pmss2}, PlotStyle->Thick, PlotMarkers->{{"\[FilledSmallCircle]",25},{"\[FilledSquare]",20}}, AxesLabel->{d, Subscript[y,4]}, GridLines->Automatic, LabelStyle->Directive[Bold, Large], 
ImageSize->Large, PlotLegends->Placed[PointLegend[{"Wetterich", "Exponential"}, LegendLabel->"Cutoff function", LegendMarkerSize->20, LegendLayout->"Row"],Bottom]];


Export["~/Desktop/y4d.png",p];


PlotRegulatorDependence[data_,label_,title_]:= Block[{pms, p, plotLabel, norm, pairs, inter, pms2,MakeLabel},
pms = PMS[data];

plotLabel="PMS \[Alpha]="<>ToString[NumberForm[pms[[1]],{3,2}]]<>" "<>label<>"="<>ToString[NumberForm[pms[[2]],{6,5}]];

p = ListPlot[data, AxesLabel->{"\[Alpha]",label},ImageSize->Large, 
		PlotLabel->plotLabel, Ticks->Automatic,GridLines->Automatic];

Export["~/Desktop/"<>title<>".png",p];
]


PlotRegulatorDependence[exps[y2,3],"y2", "y2(d=3)"]


allExponents = NestAssociation[allExponents];
exps = {\[Nu],\[Eta],e1,y};
(*ticks = Association[\[Eta]-> {Table[2+0.25i,{i,0,4}],Table[0.1i,{i,0,4}]},e1-> {Table[2+0.25i,{i,0,4}],Table[0.5i,{i,0,4}]}];*)


(* ::Section::Closed:: *)
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



(* ::Subsubsection::Closed:: *)
(*2D Ising model*)


{fps, exps}=Import["/home/andrzej/FunctionalRenormalization/ising_functional_smooth/fixed_points.wdx"];
myExps = NestAssociation[exps];


PMS[exp_,func_:Re, start_:1] := Block[{norm=Normal[myExps[exp]], pairs, inter,pms},
norm = Select[norm, Abs[Im[#[[2]]]]<1&];
pairs =  Map[{norm[[#,1]],func[norm[[#,2]]]}&,Range[Length[norm]]];

inter = Interpolation[pairs];
pms = FindRoot[inter'[\[Alpha]]==0,{\[Alpha],start}];

Return[{\[Alpha], inter[\[Alpha]]}/.pms]];


PlotRegulatorDependence[exp_, f_: Identity]:= Block[{label = ToString[exp], pms, p, plotLabel, norm, pairs, inter, data, pms2,MakeLabel},
If[exp == e2, label = "\[CapitalOmega]"];
If[MemberQ[{Re,Im,Abs},f], label = ToString[f]<>"("<>label<>")";];
pms = PMS[exp, f];

plotLabel="PMS \[Alpha]="<>ToString[NumberForm[pms[[1]],{3,2}]]<>" "<>label<>"="<>ToString[NumberForm[pms[[2]],{6,5}]];
If[exp == e2,
	norm = Normal[myExps[e2]];
	norm = Select[norm, Abs[Im[#[[2]]]]<1&];
	MakeLabel[pms_] := "PMS("<>ToString[f]<>") \[Alpha]="<>ToString[NumberForm[pms[[1]],{3,2}]]<>" \[CapitalOmega]="<>ToString[NumberForm[inter[pms[[1]]],{6,5}]];
	pairs =  Map[{norm[[#,1]],Abs[Re[norm[[#,2]]]]+I Abs[Im[norm[[#,2]]]]}&,Range[Length[norm]]];
	inter = Interpolation[pairs];
	plotLabel=MakeLabel[pms];
	If[f==Re, 
		pms2 = PMS[exp, f, 7];
		plotLabel = plotLabel<>"\n"<>MakeLabel[pms2]];
];

data=f[Select[myExps[exp],Abs[Im[#]]<1&]];
If[exp==e2, data = Abs[f[Select[myExps[exp],Abs[Im[#]]<1&]]]];

If[exp==e2,
	If[f==Im,
	data = {data,{{1.2,0.2397},{1.3,0.2584},{1.46,0.2865},{1.56,0.3028},{1.7,0.3238},{1.8,0.3375},{1.9,0.3503}}}];
	If[f==Re,
	data = {data,{{0.6,1.6247},{1,1.6302},{1.1,1.6350},{1.2,1.6402},{1.3,1.6455},{1.46,1.6541},{1.56,1.6592},{1.7,1.6662},{1.8,1.6709},{1.9,1.6753},{2.2,1.6871},{2.6,1.6996},{2.8,1.7046},{3,1.7089},{3.2,1.7125},{3.8,1.7204},{4.4,1.7247},{5,1.7265},{5.6,1.7266},{6.2,1.7255},{6.8,1.7236}}}
	];

	p = ListPlot[data, AxesLabel->{"\[Alpha]",label},ImageSize->Large, PlotLegends->{"Mine","Bertrand's"}, 
		PlotLabel->plotLabel, Ticks->Automatic,GridLines->Automatic]];

If[MemberQ[{e1,\[Eta],\[Nu]},exp],
	p = ListPlot[data, AxesLabel->{"\[Alpha]",label},ImageSize->Large, 
		PlotLabel->plotLabel, Ticks->Automatic,GridLines->Automatic];];

Export["~/Desktop/Ising/"<>label<>".png",p];
]


PlotRegulatorDependence[e1]
PlotRegulatorDependence[\[Eta]]
PlotRegulatorDependence[\[Nu]]
PlotRegulatorDependence[e2,Re]
PlotRegulatorDependence[e2,Im]


(* ::Subsubsection::Closed:: *)
(*Integral configurations optimization*)


(* ::Text:: *)
(*Importing fixed point data for optimization*)


{fps, exps}=Import["/home/andrzej/FunctionalRenormalization/ising_functional_smooth/fixed_points.wdx"];


(* ::Text:: *)
(*Selecting fixed point and constants (N=1, d=2, \[Alpha]=1.9)*)


nested = NestAssociation[exps];


sel = KeySelect[nested, MemberQ[{e1,e2,e3,e4,e5},#]&];


ListPlot[Re[Values[sel]]]


Unprotect[{gridSize,zNormalization}];
gridSize = 150;
zNormalization = 93;


test = {v[0]->-1.87271757402827896562698322136918954822`24.,v[1]->-1.87152280137413683508977557074641652839`24.,v[2]->-1.87029681826318405422786271424941724907`24.,v[3]->-1.86903830427663563266806241890441631341`24.,v[4]->-1.86774586984118569430271210036373019018`24.,v[5]->-1.86641804585809172236082715956129091816`24.,v[6]->-1.86505327776708459238097041518495613062`24.,v[7]->-1.86364991943199526712260837272448850123`24.,v[8]->-1.86220622659546324099440285067139235628`24.,v[9]->-1.86072034982840925556329935677602829763`24.,v[10]->-1.85919032692554247674575821710755933193`24.,v[11]->-1.85761407470153959331856493917084813069`24.,v[12]->-1.85598938014177882630690276848893655231`24.,v[13]->-1.85431389085997603198966177826448647247`24.,v[14]->-1.85258510481361114417433222089261101799`24.,v[15]->-1.85080035922701645324710539062183783199`24.,v[16]->-1.84895681867169919530481310167643889805`24.,v[17]->-1.84705146225417952884266441117685740696`24.,v[18]->-1.84508106986368280805379691599312543049`24.,v[19]->-1.84304220743584596295829592609063195912`24.,v[20]->-1.84093121119468406303576649980779643905`24.,v[21]->-1.83874417084402059204327401652919995675`24.,v[22]->-1.83647691169213725751753504944688771239`24.,v[23]->-1.83412497571040051579096044438702049976`24.,v[24]->-1.83168360154906653148267489405299465049`24.,v[25]->-1.8291477035624912593189973392735452276`24.,v[26]->-1.82651184993285128884287269734193833584`24.,v[27]->-1.82377024002760160934393198825454658721`24.,v[28]->-1.82091668118271530069932398594696274868`24.,v[29]->-1.81794456517271684860765768458541983255`24.,v[30]->-1.81484684471095855300684741907532647271`24.,v[31]->-1.81161601042052108865275104165325028226`24.,v[32]->-1.80824406882802421735995379279335427334`24.,v[33]->-1.80472252205912477375857212789975408407`24.,v[34]->-1.80104235005389009955232574513167197951`24.,v[35]->-1.79719399626912017608402231241790528851`24.,v[36]->-1.7931673579872526572134284365685440531`24.,v[37]->-1.78895178249896712509367221751556724737`24.,v[38]->-1.78453607055673499464327523484410366869`24.,v[39]->-1.77990848859314598944594708362354666001`24.,v[40]->-1.77505679124069197454935144177509502649`24.,v[41]->-1.76996825565508523024269934340383583414`24.,v[42]->-1.76462972900613105962653671146393149741`24.,v[43]->-1.75902769023264724957772601056609238826`24.,v[44]->-1.75314832673832502253311695882743461548`24.,v[45]->-1.74697762611915666199194448408794483467`24.,v[46]->-1.74050148225894402080695733128508635996`24.,v[47]->-1.73370581422531389168928154416982717549`24.,v[48]->-1.72657669538618290689502011278962963177`24.,v[49]->-1.71910048911337539235179895702878027183`24.,v[50]->-1.71126398643827167766537351185761570599`24.,v[51]->-1.70305454018404645964151585914624267759`24.,v[52]->-1.69446018953619252832345309401091728117`24.,v[53]->-1.68546976883276212750191470745492929703`24.,v[54]->-1.67607299463362617114297680782012225978`24.,v[55]->-1.66626052589265797072480679474928567409`24.,v[56]->-1.65602399327821108134472273693012901785`24.,v[57]->-1.64535599527498470712931194869167958892`24.,v[58]->-1.63425006051194269296955290254452513074`24.,v[59]->-1.62270057762129815152304381306064634928`24.,v[60]->-1.61070269566081859567983783681925268779`24.,v[61]->-1.59825219956435808867073797254044417018`24.,v[62]->-1.58534536610527630689138793024330136702`24.,v[63]->-1.57197880640305904435233002325979550127`24.,v[64]->-1.5581493010739104291107663564924230019`24.,v[65]->-1.5438536337730668356624550981147530906`24.,v[66]->-1.52908842818994730141970449843872059247`24.,v[67]->-1.51384999264690197493374372948318184877`24.,v[68]->-1.49813417543058633330394879656904483054`24.,v[69]->-1.48193623295251531516461440376035593398`24.,v[70]->-1.46525071187194279434619886583986408778`24.,v[71]->-1.44807134547487182612217585707227281658`24.,v[72]->-1.43039096391831598732234884177750018917`24.,v[73]->-1.41220141742845613617006382704417583368`24.,v[74]->-1.39349351117855481813862432291896446764`24.,v[75]->-1.3742569503501500333346936576831427704`24.,v[76]->-1.35448029377601899028308637199887971785`24.,v[77]->-1.33415091455058201915567793416000120969`24.,v[78]->-1.31325496604855360278477321229832451465`24.,v[79]->-1.29177735189410175944087798128183463009`24.,v[80]->-1.26970169855242438454742274231064863898`24.,v[81]->-1.24701032935913887095328417079029168156`24.,v[82]->-1.22368423894945391182061771120820128493`24.,v[83]->-1.19970306719116774659934949839012337265`24.,v[84]->-1.17504507185817018071730981366444240766`24.,v[85]->-1.14968709940145395211463742422811216623`24.,v[86]->-1.12360455328136953802952635469877545506`24.,v[87]->-1.0967713594178139853588181290736239438`24.,v[88]->-1.06915992839481470850804428553647012175`24.,v[89]->-1.04074111412361251730332687747535082435`24.,v[90]->-1.01148416872518287196990864050013603189`24.,v[91]->-0.98135669344058452912776364208814438125`24.,v[92]->-0.95032458541701427607474930897740708156`24.,v[93]->-0.91835198025033066728520499693445501697`24.,v[94]->-0.88540119019233426278142560548086536945`24.,v[95]->-0.85143263795437150851495706414274319425`24.,v[96]->-0.81640478605883986395028583765760711716`24.,v[97]->-0.78027406170775077553930754840849417891`24.,v[98]->-0.74299477715336045044418524905113592264`24.,v[99]->-0.70451904557058920771020365863407342253`24.,v[100]->-0.66479669244498932491012245271513811796`24.,v[101]->-0.62377516250375813523248174251620304968`24.,v[102]->-0.58139942223100604993831211626706363818`24.,v[103]->-0.53761185802237543886827878017158273809`24.,v[104]->-0.49235217004829111786319194817537580736`24.,v[105]->-0.44555726190966805800243973914842094656`24.,v[106]->-0.39716112618481214639880398020522561644`24.,v[107]->-0.34709472598148131518789953535494678591`24.,v[108]->-0.29528587262353981484688544807651867808`24.,v[109]->-0.24165909961721286515557783011698375044`24.,v[110]->-0.18613553305747472682183186144386093752`24.,v[111]->-0.12863275865039462236854477042622119501`24.,v[112]->-0.06906468554211215565475344454591176394`24.,v[113]->-0.00734140715928714737029241148986878687`24.,v[114]->0.05663094072087790385507803484062338673`24.,v[115]->0.12295032444086958918016759289720159839`24.,v[116]->0.19171895923254647766729845175739212911`24.,v[117]->0.26304345803655533480206647821337502231`24.,v[118]->0.33703498211483236040777298257396596732`24.,v[119]->0.41380939319485641458831004970768599116`24.,v[120]->0.49348740688213013418176384898686890041`24.,v[121]->0.57619474707789424121340441547324415755`24.,v[122]->0.66206230114252813259771061484337037982`24.,v[123]->0.751226275551616832126102954611499502`24.,v[124]->0.84382835180138599247351968898061015199`24.,v[125]->0.94001584233318486820446053512255519646`24.,v[126]->1.03994184626293338975853499393682568461`24.,v[127]->1.14376540472087826546412291524557822947`24.,v[128]->1.25165165562948786538621111980283400086`24.,v[129]->1.36377198777264638782223379097962428586`24.,v[130]->1.4803041940372008298008640550152353001`24.,v[131]->1.60143262373801503979799058109877096087`24.,v[132]->1.72734833396957358933211592721969840494`24.,v[133]->1.85824923996037693720011464525570126347`24.,v[134]->1.99434026444035424480734771461625858996`24.,v[135]->2.1358334860657327906924983595940489811`24.,v[136]->2.28294828697967098677297791264648378386`24.,v[137]->2.43591149961985731063555818422836496976`24.,v[138]->2.59495755291599963779316720952379073179`24.,v[139]->2.76032861804777689891211505165587212752`24.,v[140]->2.93227475396980577760489550685940616142`24.,v[141]->3.11105405289334267995787914037949339934`24.,v[142]->3.29693278609410942000273692099988221449`24.,v[143]->3.49018554988293769980992747587912185605`24.,v[144]->3.69109541328873452880519997001180577911`24.,v[145]->3.89995406425853892410968716914865349062`24.,v[146]->4.11706196331930249818302282496428000961`24.,v[147]->4.34272848511082781417891380269433792223`24.,v[148]->4.57727209027847925912568496723526079807`24.,v[149]->4.8210204430489517016256899543735011035`24.,v[150]->5.07431062738771659015826803456583921913`24.,zs[0]->1.25910828872524484446388447342560910896`24.,zs[1]->1.26054055472692133788355399136786387382`24.,zs[2]->1.26196858751936513248266366663742070997`24.,zs[3]->1.26339054430425640654853285803711758454`24.,zs[4]->1.26480439380006129432922653882166952997`24.,zs[5]->1.26620789831182079130768510683000818843`24.,zs[6]->1.26759859420819949310792636154815005748`24.,zs[7]->1.2689737705945992741026994151269240949`24.,zs[8]->1.27033044605542521121420530213023752295`24.,zs[9]->1.27166534332122236046604470567650495608`24.,zs[10]->1.27297486171500348684967649304943050358`24.,zs[11]->1.27425504723438377574911289076756343153`24.,zs[12]->1.27550156013274197050152571468333273531`24.,zs[13]->1.27670963987500556323041336635061186361`24.,zs[14]->1.27787406736358925331086489670608022162`24.,zs[15]->1.27898912435966044600310533205945412852`24.,zs[16]->1.28004855006688584156684507468786278335`24.,zs[17]->1.28104549490229291980827814246681167923`24.,zs[18]->1.28197247155563895043595492948631576399`24.,zs[19]->1.28282130353918100154983361633685335771`24.,zs[20]->1.2835830715591804946658310703371941337`24.,zs[21]->1.28424805820480644824441434245239832881`24.,zs[22]->1.28480569165599565046149096458716398754`24.,zs[23]->1.28524448936656653693000227963266158994`24.,zs[24]->1.28555200299012085346364388640365844862`24.,zs[25]->1.28571476619162432227506768892175406616`24.,zs[26]->1.28571824743397890046624784201772367166`24.,zs[27]->1.28554681035169463479072859658768068962`24.,zs[28]->1.28518368492525293546734586697144501845`24.,zs[29]->1.28461095334738817204658909812204430857`24.,zs[30]->1.28380955521647837479979734476683654036`24.,zs[31]->1.28275931748235921682490666356996017827`24.,zs[32]->1.28143901537197166015554019833789043991`24.,zs[33]->1.27982647128397768967396384848148575847`24.,zs[34]->1.27789869928804328383247376634435333307`24.,zs[35]->1.27563210329483504443735482197439002165`24.,zs[36]->1.27300273704715794876066197506878760434`24.,zs[37]->1.26998663366300281419093405886344517425`24.,zs[38]->1.26656021135613211470225733047298190942`24.,zs[39]->1.26270075997608119173238155418810295411`24.,zs[40]->1.25838700996282442006672199207856517089`24.,zs[41]->1.25359978105751909393875264127218325853`24.,zs[42]->1.24832270258711649657821579102841557734`24.,zs[43]->1.24254299041675018480882009081786596814`24.,zs[44]->1.23625225799429786723086077876770078263`24.,zs[45]->1.22944733078130756216819481761545528897`24.,zs[46]->1.22213102551153846315441446574312283662`24.,zs[47]->1.21431284911640472559968844135797075209`24.,zs[48]->1.20600956794173461220506037528406973192`24.,zs[49]->1.19724559721266709777367023935622386102`24.,zs[50]->1.18805316456968469205355075774249755617`24.,zs[51]->1.17847221047948813554170446348703749109`24.,zs[52]->1.16855000237459145145194573226216862934`24.,zs[53]->1.15834045767621601254043836660538632134`24.,zs[54]->1.14790319179667336506295952747518404221`24.,zs[55]->1.13730232855175684643392096662887343666`24.,zs[56]->1.12660512957200236476554400588434756611`24.,zs[57]->1.11588051382635307707774822622801175105`24.,zs[58]->1.10519754636712274346269036350091761161`24.,zs[59]->1.09462397589743384147017237347827039317`24.,zs[60]->1.08422489387572232483796506506432422049`24.,zs[61]->1.07406157477705364867615324713593811862`24.,zs[62]->1.06419053978471711347332636391835986429`24.,zs[63]->1.05466286693636676011844523121911208375`24.,zs[64]->1.04552375191093886722756435600921872771`24.,zs[65]->1.03681230713581761997896718149418302859`24.,zs[66]->1.02856157401751503813386497117259637022`24.,zs[67]->1.02079871445925371577837259371698953215`24.,zs[68]->1.0135453433999326040470879772906470101`24.,zs[69]->1.00681796338680424138261729470930584531`24.,zs[70]->1.00062846438525261454098792190059894038`24.,zs[71]->0.99498465623200348109790629853118304194`24.,zs[72]->0.98989080648784945580961488829009779435`24.,zs[73]->0.98534816220783619004877368576729738354`24.,zs[74]->0.98135543976206284633123900312040573392`24.,zs[75]->0.97790927193265276452224372466987674959`24.,zs[76]->0.97500460586932507181245967868916266112`24.,zs[77]->0.97263504902530829842310507202276420718`24.,zs[78]->0.97079316292916908945868071845244288249`24.,zs[79]->0.96947070664871203581308734480692140146`24.,zs[80]->0.968658833174820450508834620003933805`24.,zs[81]->0.96834824281160567241607377551952581701`24.,zs[82]->0.96852929811734085480287625163375141087`24.,zs[83]->0.96919210509995147444624989031135889388`24.,zs[84]->0.9703265653175374266500220950497545025`24.,zs[85]->0.97192240333826097277824117590364370977`24.,zs[86]->0.97396917372901829527778493853756809284`24.,zs[87]->0.97645625140883622205258893478110522519`24.,zs[88]->0.97937280884954599474912056007008119064`24.,zs[89]->0.98270778325239970502968600275841565298`24.,zs[90]->0.98644983648727075592225373369820938983`24.,zs[91]->0.99058731025810194520194860622623936925`24.,zs[92]->0.99510817865789978934497614091973040068`24.,zs[94]->1.0052498695592705748159893720999032856`24.,zs[95]->1.01084437462632089287902235619768099464`24.,zs[96]->1.01676955306827459471064942646876789544`24.,zs[97]->1.02301085639983751860619877039066929415`24.,zs[98]->1.02955311819696835632042171075573671931`24.,zs[99]->1.03638052853133162948414113432340359943`24.,zs[100]->1.04347661496592725900604972553451703816`24.,zs[101]->1.05082423052998571212123422127075234832`24.,zs[102]->1.05840554898353252869056548015969856722`24.,zs[103]->1.06620206758801197160605691537122885169`24.,zs[104]->1.07419461751789733220567449563183611244`24.,zs[105]->1.08236338197795608320264327004182430443`24.,zs[106]->1.09068792203017318243875780473076136723`24.,zs[107]->1.09914721008140143920694165443704710726`24.,zs[108]->1.10771967093552276482926868713477455204`24.,zs[109]->1.11638323027003999040940371046108351609`24.,zs[110]->1.1251153703542968855648340701782349997`24.,zs[111]->1.13389319278273049168135600401186889261`24.,zs[112]->1.14269348794967280373814108369955845536`24.,zs[113]->1.15149281094053462040629822651939519677`24.,zs[114]->1.16026756345645643793952702375361697666`24.,zs[115]->1.16899408132497115143209013382558362387`24.,zs[116]->1.17764872707777720509685967935598420851`24.,zs[117]->1.18620798699891539741230894948101585851`24.,zs[118]->1.1946485719636966692871974694054135288`24.,zs[119]->1.20294752130250851835773238241887427711`24.,zs[120]->1.21108230883658842519393832551131729479`24.,zs[121]->1.21903095014794200325809949123621771435`24.,zs[122]->1.22677211006612033845873700031449487962`24.,zs[123]->1.23428520928409737054027776674887902494`24.,zs[124]->1.24155052895760499362230273888672238891`24.,zs[125]->1.24854931210047562130980737953645532431`24.,zs[126]->1.25526386056600181200653412761156923186`24.,zs[127]->1.26167762640378099096574643915463651581`24.,zs[128]->1.26777529640508512746194086335102309344`24.,zs[129]->1.27354286869884801984255096534927968189`24.,zs[130]->1.27896772033541602432897479908247724522`24.,zs[131]->1.284038664895867978511028589312227556`24.,zs[132]->1.2887459992896454561191208513211248369`24.,zs[133]->1.29308153905019472339565846780201666207`24.,zs[134]->1.29703864160419743660237763998899510239`24.,zs[135]->1.30061221717090201416615238061586833155`24.,zs[136]->1.30379872713948242184265474481427561154`24.,zs[137]->1.30659616997015208151435417881335783315`24.,zs[138]->1.30900405485800111985259281753036466638`24.,zs[139]->1.31102336361450617427333377110338475239`24.,zs[140]->1.3126565013007993543068650110482707826`24.,zs[141]->1.31390723671213208130212464064889394885`24.,zs[142]->1.31478063264232699325282898325172685471`24.,zs[143]->1.31528296984912819983397298367672256073`24.,zs[144]->1.31542165866565266202264590508334972237`24.,zs[145]->1.3152051558927526602114269740673325748`24.,zs[146]->1.31464285456125168429633098037042880502`24.,zs[147]->1.31374501087395415900593851944524795224`24.,zs[148]->1.31252260223629425531056903885455178276`24.,zs[149]->1.310987284499670791344878607171177528`24.,zs[150]->1.30915121752246950784505861632606559761`24.,zp[1]->1.25795782235583664939251378307539554786`24.,zp[2]->1.25673788775674013078321377422236941907`24.,zp[3]->1.25544622038522542591659362853882473767`24.,zp[4]->1.25408049044128985108574748066113907823`24.,zp[5]->1.25263830257276020246581296143136276102`24.,zp[6]->1.25111719545374485408426935181972732688`24.,zp[7]->1.24951464163873026929184255024705148369`24.,zp[8]->1.24782804770414532684420119132992502199`24.,zp[9]->1.24605475472058259573784630808265097325`24.,zp[10]->1.24419203910704306664744911123466302302`24.,zp[11]->1.24223711392514798829394575516118461477`24.,zp[12]->1.24018713067817463518594032669257748385`24.,zp[13]->1.23803918168723450673923492241605039249`24.,zp[14]->1.23579030312495077068688716075853685982`24.,zp[15]->1.23343747879559229935898030354032600825`24.,zp[16]->1.23097764475973852187584791769672963059`24.,zp[17]->1.22840769491110004204273526544497365207`24.,zp[18]->1.22572448762297607224486419112386596225`24.,zp[19]->1.22292485359180460020981655457805871684`24.,zp[20]->1.22000560501509537667888379725904128696`24.,zp[21]->1.21696354625037963082406834850810750604`24.,zp[22]->1.21379548611020410867939050616070275701`24.,zp[23]->1.21049825195504841332493822329958491641`24.,zp[24]->1.20706870575060507720590396193664929531`24.,zp[25]->1.20350376225720147440049128604130264108`24.,zp[26]->1.19980040951612683519044235973067337301`24.,zp[27]->1.19595573178889671896204195266868120532`24.,zp[28]->1.19196693508944490488504213827495554463`24.,zp[29]->1.18783137542404663833740395452883591945`24.,zp[30]->1.18354658981740146812620250613605532239`24.,zp[31]->1.17911033015353020150035106061987392688`24.,zp[32]->1.17452059979469090961334074925523847113`24.,zp[33]->1.16977569285819241257263041605822880539`24.,zp[34]->1.16487423592786928745404271824028078331`24.,zp[35]->1.15981523185275079064921201897448777456`24.,zp[36]->1.15459810513973627894306882130379315354`24.,zp[37]->1.14922274828092971093589987865784638416`24.,zp[38]->1.14368956817267963630496431373743188641`24.,zp[39]->1.13799953158781500268613478608828788306`24.,zp[40]->1.13215420846357895686825680025715611368`24.,zp[41]->1.12615581157726602921363279424107680924`24.,zp[42]->1.12000723101500959713244840930602517156`24.,zp[43]->1.11371206171521344655260627264725154487`24.,zp[44]->1.10727462230777781352523541290678325941`24.,zp[45]->1.10069996349528764485513650070722713026`24.,zp[46]->1.09399386435286243463406358485049016955`24.,zp[47]->1.08716281517495795723130861120541285033`24.,zp[48]->1.08021398587759377019836334315018971072`24.,zp[49]->1.07315517946946209451398617892235270758`24.,zp[50]->1.06599477071746257460915082785773700245`24.,zp[51]->1.05874163081873546050642218822402300456`24.,zp[52]->1.0514050396056556204006998575961019634`24.,zp[53]->1.04399458749560131587056271815859228763`24.,zp[54]->1.03652006999239635598990843140196703264`24.,zp[55]->1.02899137799329096731756695363534953727`24.,zp[56]->1.02141838740790408497540240418178231516`24.,zp[57]->1.01381085162607990350067323946177599784`24.,zp[58]->1.00617830017572515003017163354931083156`24.,zp[59]->0.99852994650919348068668369191253377186`24.,zp[60]->0.99087460728876481452185613364302190245`24.,zp[61]->0.98322063486443903496609485650421403801`24.,zp[62]->0.97557586391377517276838708170993844731`24.,zp[63]->0.96794757250541811646947146676544531316`24.,zp[64]->0.96034245720789819116345522933378529309`24.,zp[65]->0.95276662133194818798286925222144061372`24.,zp[66]->0.94522557499046094008331093315135335879`24.,zp[67]->0.93772424539204375609470910509105785235`24.,zp[68]->0.93026699564531966424959606690369148597`24.,zp[69]->0.92285765032537477035773217669709996274`24.,zp[70]->0.91549952611900170554987662583669549095`24.,zp[71]->0.90819546599747128332018641587629355261`24.,zp[72]->0.90094787554108161938325498718802564878`24.,zp[73]->0.89375876023803623364672564398853809618`24.,zp[74]->0.88662976278447078147811311840377107568`24.,zp[75]->0.87956219961003812067247633865780065636`24.,zp[76]->0.87255709603571453117203808567892512694`24.,zp[77]->0.86561521963221219893708997405833107783`24.,zp[78]->0.85873711148617000753247416804013030098`24.,zp[79]->0.851923115196828717649035289943344933`24.,zp[80]->0.84517340351930159326571734789685984418`24.,zp[81]->0.83848800264384813374268288551869704973`24.,zp[82]->0.83186681415623644674124198430494692589`24.,zp[83]->0.82530963476498293353355978459248253208`24.,zp[84]->0.81881617390956918698728092089502554269`24.,zp[85]->0.81238606938204206647509175988375034485`24.,zp[86]->0.80601890110482100896531822268805099756`24.,zp[87]->0.79971420321188516510392138451001066716`24.,zp[88]->0.7934714745803121216265777077688213406`24.,zp[89]->0.7872901879556310023968672740825339509`24.,zp[90]->0.78116979780862733865190068690835140499`24.,zp[91]->0.77510974705387139140829772155421809449`24.,zp[92]->0.76910947275193062171132366573936624017`24.,zp[93]->0.76316841090841730204531230349972579372`24.,zp[94]->0.75728600047404140675945346608367776224`24.,zp[95]->0.75146168664092114166484314330736554796`24.,zp[96]->0.74569492352171142203831720072311785955`24.,zp[97]->0.73998517628975333281775499945326723419`24.,zp[98]->0.73433192285049508171460876981302899439`24.,zp[99]->0.7287346551069295560624877104855891588`24.,zp[100]->0.72319287987475889667896271793338919658`24.,zp[101]->0.71770611949644386080378231794773006805`24.,zp[102]->0.71227391219722898492464362932965405413`24.,zp[103]->0.70689581222065313786978542513948652964`24.,zp[104]->0.70157138977595602235017971107388955772`24.,zp[105]->0.69630023082517017931250203744966384869`24.,zp[106]->0.69108193673353961258946449884730056855`24.,zp[107]->0.6859161238032235358307401097766661822`24.,zp[108]->0.68080242270701843205776797065873927296`24.,zp[109]->0.67574047783605268100031330516084547756`24.,zp[110]->0.6707299465730614528531270023413552547`24.,zp[111]->0.66577049850091773907133295102380577245`24.,zp[112]->0.66086181455455656039588635397239578561`24.,zp[113]->0.65600358612325752364391033895140184561`24.,zp[114]->0.65119551410941565892468942127088437234`24.,zp[115]->0.64643730794939747460065582166095200143`24.,zp[116]->0.64172868460181046064541696205752093208`24.,zp[117]->0.63706936750846899917608773481521810378`24.,zp[118]->0.63245908553347490369799987747260764913`24.,zp[119]->0.6278975718861026312277940271986070093`24.,zp[120]->0.62338456303354359474126308683752827961`24.,zp[121]->0.61891979760997800273763357099273669899`24.,zp[122]->0.61450301532886543892440903088118378697`24.,zp[123]->0.61013395590573925296628763914221054774`24.,zp[124]->0.6058123579991210596133722468432899457`24.,zp[125]->0.60153795817741131109618600456774620332`24.,zp[126]->0.59731048991973647002941958009115366013`24.,zp[127]->0.59312968265872502815517718104248234478`24.,zp[128]->0.58899526087303179534338743995345852069`24.,zp[129]->0.58490694323712691607456793699611486032`24.,zp[130]->0.58086444183541328956819124246906701593`24.,zp[131]->0.57686746144713941779000168246694000191`24.,zp[132]->0.57291569890784524180002867011956199689`24.,zp[133]->0.5690088425522317877853186606162631007`24.,zp[134]->0.56514657174240068767198529897622165845`24.,zp[135]->0.56132855648438899746562597688727243803`24.,zp[136]->0.55755445713485236495733271063267194534`24.,zp[137]->0.55382392419865030396903459887311872446`24.,zp[138]->0.55013659821698913535984041813895324781`24.,zp[139]->0.54649210974468956865045560595119002149`24.,zp[140]->0.54289007941417166883042550450088141581`24.,zp[141]->0.53933011808255490249173772322699233709`24.,zp[142]->0.53581182705823729464949255816529410038`24.,zp[143]->0.53233479839980909509470815975494020647`24.,zp[144]->0.52889861528748844785872092782196388425`24.,zp[145]->0.5255028524453548770280418518960106074`24.,zp[146]->0.52214707664485528510426284520012641123`24.,zp[147]->0.5188308471993315294651788288004388574`24.,zp[148]->0.51555371661634187870061163753824970353`24.,zp[149]->0.51231523105375050427732287754340297774`24.,zp[150]->0.50911493119447451666912000144677833317`24.,\[Eta]->0.23590572598007892234469142698930640866`24.};


dd = 200/100;
{aa,rr} = {2,144/10};
myfp = test;


myintegrands = integrandsListIsotropic[[All,3]];


myintegrands = myintegrands/.{y^2-> y^4, y^(-1+d)->2 y^(2d-1)};


Unprotect[IntegralConfigurations];
IntegralConfigurations = {{0,2/3,40},{2/3,12/10,40},{12/10,3,30}};
Protect[IntegralConfigurations];


(* ::Text:: *)
(*Sparsely probing \[Rho] grid *)


sparse = 10;
integrands1 = Flatten[Table[NumericDerivatives[myintegrands/.{n->10/10,\[Rho] -> (1+ sparse (i-1)) eps}/.regulatorReplacement]//.ConstRep[dd,rr,aa]/.myfp,{i,1,gridSize/sparse}]];
integrands[q_] = integrands1 /.y->q;


(* ::Text:: *)
(*Function determining maximal error between an exact integral and GLIntegral on specified integral with specified number of points*)


Diff[config_]:= Block[{exact,glint, int2},
	exact = NIntegrate[integrands[q],{q,config[[1]],config[[2]]}, WorkingPrecision->24];
	glint = GLIntegral[integrands,config[[1]],config[[2]],config[[3]]];
	Return[Max[Abs[exact-glint]]];
]


DiffExpand[config_]:= Block[{exact,glint, int2},
	If[config[[1]]==0.,
	int2[q_] = Simplify[(integrands[q]/.{Exp[x_] ->1+x+x^2/2+x^3/6+x^4/24+x^5/120+x^6/720+x^7/5040+x^8/40320+x^9/362880+x^10/3628800})];
	exact = NIntegrate[int2[q],{q,config[[1]],config[[2]]}, WorkingPrecision->30, MaxRecursion->50, AccuracyGoal->12, PrecisionGoal->12],
	exact = NIntegrate[integrands[q],{q,config[[1]],config[[2]]}, WorkingPrecision->30, MaxRecursion->50, AccuracyGoal->12, PrecisionGoal->12];
		];
	glint = GLIntegral[integrands,config[[1]],config[[2]],config[[3]]];
	Return[Max[Abs[exact-glint]]];
]


Quiet[Map[DiffExpand,{{0,1/20,30},{1/20,3,34},{3,6,12}}],NIntegrate::precw]


Quiet[Map[DiffExpand,IntegralConfigurations],NIntegrate::precw]


f = Simplify[integrands[q]/.{Exp[x_] ->1+x+x^2/2+x^3/6+x^4/24+x^5/120+x^6/720+x^7/5040+x^8/40320+x^9/362880+x^10/3628800}];


(* ::Text:: *)
(*Optimizing number of points for integrals on each specified interval not to exceed the threshold *)


Optimize[configs_] := Block[{current = configs, last = configs, diffs, threshold=10^-11, pos, i,Opti},
	Opti[conf_]:= Block[{diff, new=conf, growing},
		For[i=0,i<50,i++,
			diff = DiffExpand[new];
			If[i==0, If[diff>threshold, growing=False, growing=True]];
			
			If[diff>threshold,
				If[growing || new[[3]]+1>=MAXPOINTSLEGENDRE, 
					new[[3]]+=1; 
					Return[new]];
				new[[3]] += 1,
				
				If[!growing || new[[3]]==1, Return[new]];
				new[[3]] -= 1;
			];
		];
		Return[new];	 
	];
	
	Return[Map[Opti,configs]];
];


optimized = <| |>;


IsValid[c_]:= Block[{},
Return[AllTrue[Range[2,Length[c]],c[[#,1]]==c[[#-1,2]]&]]];


conf = {{0,0.1,4},{0.1, 1., 10},{1., 3., 10},{3., 6., 10},{6., 10., 2}};
optimized[conf] = Optimize[conf];
Print[Total[optimized[conf][[All,3]]]]


IntegralConfigurations


IntegralConfigurations


DiffExpand[{3,6,15}]


DiffExpand[{0,1/10,6}]


DiffExpand[IntegralConfigurations[[1]]]


conf = {{0,5.5,40}};


diffs = Map[Diff[{0,5.5,#}]&,Range[40]];


(1-q^2)^2 HeavisideTheta[1-q^2]


ListPlot[diffs, PlotRange->All]


Sum[DiffExpand[IntegralConfigurations[[i]]],{i,1,4}]


f = FullSimplify[integrands1[[3]]];


ders = Table[D[f,{y,2 n}],{n,1,10}];


maxs = Map[NMaxValue[{Abs[#],y<=4.2, y>=0.01},y]&,ders];


maxs2 = Map[NMaxValue[{Abs[#],y<=4.2, y>=1},y]&,ders];


maxs3 = Map[NMaxValue[{Abs[#],y<=2, y>=.1},y]&,ders];


Table[4.2^(2n+1) (n!)^4/((2n+1)*((2n)!)^3)*maxs[[n]],{n,1,10}]


Table[3.2^(2n+1) (n!)^4/((2n+1)*((2n)!)^3) maxs2[[n]],{n,1,10}]


Table[1.9^(2n+1) (n!)^4/((2n+1)*((2n)!)^3) maxs3[[n]],{n,1,10}]


ders = {f};
Map[AppendTo[ders, D[ders[[-1]],{y,2}]]&,Range[10]];


ListLogPlot[diffs,PlotRange->All]


(* ::Subsubsection::Closed:: *)
(*Data from O(N) calculations*)


{fps, exps} = Import["ising2_functional_smooth/"<>FIXEDPOINTSFILE];


Keys[exps]


data = OneLevelNest[KeySelect[exps, (MatchQ[#,{d_,e_,a_}/;(MemberQ[{e1,e2,e3,e4,e5},e])]&)]];


data2 = OneLevelNest[data];


PMS[data_,func_:Re] := Block[{norm=Normal[data], pairs, inter,pms},
pairs =  Map[{norm[[#,1]],func[norm[[#,2]]]}&,Range[Length[norm]]];

inter = Interpolation[pairs];
pms = FindRoot[inter'[\[Alpha]]==0,{\[Alpha],2}];

Return[inter[\[Alpha]]/.pms]];


pmss = TransposeAssociation[NestAssociation[Quiet[Map[PMS, data]]]];


p = ListPlot[Values[pmss], PlotLegends->Map["Re("<>ToString[#]<>")"&,Keys[pmss]], GridLines->Automatic, ImageSize->Large, AxesLabel->{"d"}, LabelStyle->Directive[Bold, Medium]]


Export["ising_eigenvalues.png", p]


(* ::Section:: *)
(*O(N)-models*)


(* ::Subsection::Closed:: *)
(*Momentum integrals*)


{fps,exps} = Import["ON_functional_smooth/fixed_points.wdx"];


SelectRegulator["smooth"]
Unprotect[zNormalization];
dat = KeySelect[fps, MatchQ[{20/10,1,15/10,r_}]];
{aa, nn, dd, rr} = Keys[dat][[1]];
test = dat[[1]];
zNormalization = 0;


threshold[d_,i_,m_] := ((r[#^2]-#^2 r[#^2]) #^(d-1) /(m+#^2+r[#^2])^i /.regulatorReplacement/.\[Alpha]->2)&
threshold2[d_,i_,m_] := ((r[#^(2/d)]-#^(2/d) r[#^(2/d)]) /(m+#^(2/d)+r[#^(2/d)])^i /.regulatorReplacement/.\[Alpha]->2)&
threshold3[d_,i_,m_] := 2 #^(2d-1)((r[#^(4)]-#^(4) r[#^(4)]) /(m+#^(4)+r[#^(4)])^i /.regulatorReplacement/.\[Alpha]->2)&
threshold4[d_,i_,m_] := 4 #^(4d-1)((r[#^(8)]-#^(8) r[#^(8)]) /(m+#^(8)+r[#^(8)])^i /.regulatorReplacement/.\[Alpha]->2)&


Unprotect[{gridSize,zNormalization}];
gridSize = 150;
zNormalization = 93;


test = {v[0]->-1.87271757402827896562698322136918954822`24.,v[1]->-1.87152280137413683508977557074641652839`24.,v[2]->-1.87029681826318405422786271424941724907`24.,v[3]->-1.86903830427663563266806241890441631341`24.,v[4]->-1.86774586984118569430271210036373019018`24.,v[5]->-1.86641804585809172236082715956129091816`24.,v[6]->-1.86505327776708459238097041518495613062`24.,v[7]->-1.86364991943199526712260837272448850123`24.,v[8]->-1.86220622659546324099440285067139235628`24.,v[9]->-1.86072034982840925556329935677602829763`24.,v[10]->-1.85919032692554247674575821710755933193`24.,v[11]->-1.85761407470153959331856493917084813069`24.,v[12]->-1.85598938014177882630690276848893655231`24.,v[13]->-1.85431389085997603198966177826448647247`24.,v[14]->-1.85258510481361114417433222089261101799`24.,v[15]->-1.85080035922701645324710539062183783199`24.,v[16]->-1.84895681867169919530481310167643889805`24.,v[17]->-1.84705146225417952884266441117685740696`24.,v[18]->-1.84508106986368280805379691599312543049`24.,v[19]->-1.84304220743584596295829592609063195912`24.,v[20]->-1.84093121119468406303576649980779643905`24.,v[21]->-1.83874417084402059204327401652919995675`24.,v[22]->-1.83647691169213725751753504944688771239`24.,v[23]->-1.83412497571040051579096044438702049976`24.,v[24]->-1.83168360154906653148267489405299465049`24.,v[25]->-1.8291477035624912593189973392735452276`24.,v[26]->-1.82651184993285128884287269734193833584`24.,v[27]->-1.82377024002760160934393198825454658721`24.,v[28]->-1.82091668118271530069932398594696274868`24.,v[29]->-1.81794456517271684860765768458541983255`24.,v[30]->-1.81484684471095855300684741907532647271`24.,v[31]->-1.81161601042052108865275104165325028226`24.,v[32]->-1.80824406882802421735995379279335427334`24.,v[33]->-1.80472252205912477375857212789975408407`24.,v[34]->-1.80104235005389009955232574513167197951`24.,v[35]->-1.79719399626912017608402231241790528851`24.,v[36]->-1.7931673579872526572134284365685440531`24.,v[37]->-1.78895178249896712509367221751556724737`24.,v[38]->-1.78453607055673499464327523484410366869`24.,v[39]->-1.77990848859314598944594708362354666001`24.,v[40]->-1.77505679124069197454935144177509502649`24.,v[41]->-1.76996825565508523024269934340383583414`24.,v[42]->-1.76462972900613105962653671146393149741`24.,v[43]->-1.75902769023264724957772601056609238826`24.,v[44]->-1.75314832673832502253311695882743461548`24.,v[45]->-1.74697762611915666199194448408794483467`24.,v[46]->-1.74050148225894402080695733128508635996`24.,v[47]->-1.73370581422531389168928154416982717549`24.,v[48]->-1.72657669538618290689502011278962963177`24.,v[49]->-1.71910048911337539235179895702878027183`24.,v[50]->-1.71126398643827167766537351185761570599`24.,v[51]->-1.70305454018404645964151585914624267759`24.,v[52]->-1.69446018953619252832345309401091728117`24.,v[53]->-1.68546976883276212750191470745492929703`24.,v[54]->-1.67607299463362617114297680782012225978`24.,v[55]->-1.66626052589265797072480679474928567409`24.,v[56]->-1.65602399327821108134472273693012901785`24.,v[57]->-1.64535599527498470712931194869167958892`24.,v[58]->-1.63425006051194269296955290254452513074`24.,v[59]->-1.62270057762129815152304381306064634928`24.,v[60]->-1.61070269566081859567983783681925268779`24.,v[61]->-1.59825219956435808867073797254044417018`24.,v[62]->-1.58534536610527630689138793024330136702`24.,v[63]->-1.57197880640305904435233002325979550127`24.,v[64]->-1.5581493010739104291107663564924230019`24.,v[65]->-1.5438536337730668356624550981147530906`24.,v[66]->-1.52908842818994730141970449843872059247`24.,v[67]->-1.51384999264690197493374372948318184877`24.,v[68]->-1.49813417543058633330394879656904483054`24.,v[69]->-1.48193623295251531516461440376035593398`24.,v[70]->-1.46525071187194279434619886583986408778`24.,v[71]->-1.44807134547487182612217585707227281658`24.,v[72]->-1.43039096391831598732234884177750018917`24.,v[73]->-1.41220141742845613617006382704417583368`24.,v[74]->-1.39349351117855481813862432291896446764`24.,v[75]->-1.3742569503501500333346936576831427704`24.,v[76]->-1.35448029377601899028308637199887971785`24.,v[77]->-1.33415091455058201915567793416000120969`24.,v[78]->-1.31325496604855360278477321229832451465`24.,v[79]->-1.29177735189410175944087798128183463009`24.,v[80]->-1.26970169855242438454742274231064863898`24.,v[81]->-1.24701032935913887095328417079029168156`24.,v[82]->-1.22368423894945391182061771120820128493`24.,v[83]->-1.19970306719116774659934949839012337265`24.,v[84]->-1.17504507185817018071730981366444240766`24.,v[85]->-1.14968709940145395211463742422811216623`24.,v[86]->-1.12360455328136953802952635469877545506`24.,v[87]->-1.0967713594178139853588181290736239438`24.,v[88]->-1.06915992839481470850804428553647012175`24.,v[89]->-1.04074111412361251730332687747535082435`24.,v[90]->-1.01148416872518287196990864050013603189`24.,v[91]->-0.98135669344058452912776364208814438125`24.,v[92]->-0.95032458541701427607474930897740708156`24.,v[93]->-0.91835198025033066728520499693445501697`24.,v[94]->-0.88540119019233426278142560548086536945`24.,v[95]->-0.85143263795437150851495706414274319425`24.,v[96]->-0.81640478605883986395028583765760711716`24.,v[97]->-0.78027406170775077553930754840849417891`24.,v[98]->-0.74299477715336045044418524905113592264`24.,v[99]->-0.70451904557058920771020365863407342253`24.,v[100]->-0.66479669244498932491012245271513811796`24.,v[101]->-0.62377516250375813523248174251620304968`24.,v[102]->-0.58139942223100604993831211626706363818`24.,v[103]->-0.53761185802237543886827878017158273809`24.,v[104]->-0.49235217004829111786319194817537580736`24.,v[105]->-0.44555726190966805800243973914842094656`24.,v[106]->-0.39716112618481214639880398020522561644`24.,v[107]->-0.34709472598148131518789953535494678591`24.,v[108]->-0.29528587262353981484688544807651867808`24.,v[109]->-0.24165909961721286515557783011698375044`24.,v[110]->-0.18613553305747472682183186144386093752`24.,v[111]->-0.12863275865039462236854477042622119501`24.,v[112]->-0.06906468554211215565475344454591176394`24.,v[113]->-0.00734140715928714737029241148986878687`24.,v[114]->0.05663094072087790385507803484062338673`24.,v[115]->0.12295032444086958918016759289720159839`24.,v[116]->0.19171895923254647766729845175739212911`24.,v[117]->0.26304345803655533480206647821337502231`24.,v[118]->0.33703498211483236040777298257396596732`24.,v[119]->0.41380939319485641458831004970768599116`24.,v[120]->0.49348740688213013418176384898686890041`24.,v[121]->0.57619474707789424121340441547324415755`24.,v[122]->0.66206230114252813259771061484337037982`24.,v[123]->0.751226275551616832126102954611499502`24.,v[124]->0.84382835180138599247351968898061015199`24.,v[125]->0.94001584233318486820446053512255519646`24.,v[126]->1.03994184626293338975853499393682568461`24.,v[127]->1.14376540472087826546412291524557822947`24.,v[128]->1.25165165562948786538621111980283400086`24.,v[129]->1.36377198777264638782223379097962428586`24.,v[130]->1.4803041940372008298008640550152353001`24.,v[131]->1.60143262373801503979799058109877096087`24.,v[132]->1.72734833396957358933211592721969840494`24.,v[133]->1.85824923996037693720011464525570126347`24.,v[134]->1.99434026444035424480734771461625858996`24.,v[135]->2.1358334860657327906924983595940489811`24.,v[136]->2.28294828697967098677297791264648378386`24.,v[137]->2.43591149961985731063555818422836496976`24.,v[138]->2.59495755291599963779316720952379073179`24.,v[139]->2.76032861804777689891211505165587212752`24.,v[140]->2.93227475396980577760489550685940616142`24.,v[141]->3.11105405289334267995787914037949339934`24.,v[142]->3.29693278609410942000273692099988221449`24.,v[143]->3.49018554988293769980992747587912185605`24.,v[144]->3.69109541328873452880519997001180577911`24.,v[145]->3.89995406425853892410968716914865349062`24.,v[146]->4.11706196331930249818302282496428000961`24.,v[147]->4.34272848511082781417891380269433792223`24.,v[148]->4.57727209027847925912568496723526079807`24.,v[149]->4.8210204430489517016256899543735011035`24.,v[150]->5.07431062738771659015826803456583921913`24.,zs[0]->1.25910828872524484446388447342560910896`24.,zs[1]->1.26054055472692133788355399136786387382`24.,zs[2]->1.26196858751936513248266366663742070997`24.,zs[3]->1.26339054430425640654853285803711758454`24.,zs[4]->1.26480439380006129432922653882166952997`24.,zs[5]->1.26620789831182079130768510683000818843`24.,zs[6]->1.26759859420819949310792636154815005748`24.,zs[7]->1.2689737705945992741026994151269240949`24.,zs[8]->1.27033044605542521121420530213023752295`24.,zs[9]->1.27166534332122236046604470567650495608`24.,zs[10]->1.27297486171500348684967649304943050358`24.,zs[11]->1.27425504723438377574911289076756343153`24.,zs[12]->1.27550156013274197050152571468333273531`24.,zs[13]->1.27670963987500556323041336635061186361`24.,zs[14]->1.27787406736358925331086489670608022162`24.,zs[15]->1.27898912435966044600310533205945412852`24.,zs[16]->1.28004855006688584156684507468786278335`24.,zs[17]->1.28104549490229291980827814246681167923`24.,zs[18]->1.28197247155563895043595492948631576399`24.,zs[19]->1.28282130353918100154983361633685335771`24.,zs[20]->1.2835830715591804946658310703371941337`24.,zs[21]->1.28424805820480644824441434245239832881`24.,zs[22]->1.28480569165599565046149096458716398754`24.,zs[23]->1.28524448936656653693000227963266158994`24.,zs[24]->1.28555200299012085346364388640365844862`24.,zs[25]->1.28571476619162432227506768892175406616`24.,zs[26]->1.28571824743397890046624784201772367166`24.,zs[27]->1.28554681035169463479072859658768068962`24.,zs[28]->1.28518368492525293546734586697144501845`24.,zs[29]->1.28461095334738817204658909812204430857`24.,zs[30]->1.28380955521647837479979734476683654036`24.,zs[31]->1.28275931748235921682490666356996017827`24.,zs[32]->1.28143901537197166015554019833789043991`24.,zs[33]->1.27982647128397768967396384848148575847`24.,zs[34]->1.27789869928804328383247376634435333307`24.,zs[35]->1.27563210329483504443735482197439002165`24.,zs[36]->1.27300273704715794876066197506878760434`24.,zs[37]->1.26998663366300281419093405886344517425`24.,zs[38]->1.26656021135613211470225733047298190942`24.,zs[39]->1.26270075997608119173238155418810295411`24.,zs[40]->1.25838700996282442006672199207856517089`24.,zs[41]->1.25359978105751909393875264127218325853`24.,zs[42]->1.24832270258711649657821579102841557734`24.,zs[43]->1.24254299041675018480882009081786596814`24.,zs[44]->1.23625225799429786723086077876770078263`24.,zs[45]->1.22944733078130756216819481761545528897`24.,zs[46]->1.22213102551153846315441446574312283662`24.,zs[47]->1.21431284911640472559968844135797075209`24.,zs[48]->1.20600956794173461220506037528406973192`24.,zs[49]->1.19724559721266709777367023935622386102`24.,zs[50]->1.18805316456968469205355075774249755617`24.,zs[51]->1.17847221047948813554170446348703749109`24.,zs[52]->1.16855000237459145145194573226216862934`24.,zs[53]->1.15834045767621601254043836660538632134`24.,zs[54]->1.14790319179667336506295952747518404221`24.,zs[55]->1.13730232855175684643392096662887343666`24.,zs[56]->1.12660512957200236476554400588434756611`24.,zs[57]->1.11588051382635307707774822622801175105`24.,zs[58]->1.10519754636712274346269036350091761161`24.,zs[59]->1.09462397589743384147017237347827039317`24.,zs[60]->1.08422489387572232483796506506432422049`24.,zs[61]->1.07406157477705364867615324713593811862`24.,zs[62]->1.06419053978471711347332636391835986429`24.,zs[63]->1.05466286693636676011844523121911208375`24.,zs[64]->1.04552375191093886722756435600921872771`24.,zs[65]->1.03681230713581761997896718149418302859`24.,zs[66]->1.02856157401751503813386497117259637022`24.,zs[67]->1.02079871445925371577837259371698953215`24.,zs[68]->1.0135453433999326040470879772906470101`24.,zs[69]->1.00681796338680424138261729470930584531`24.,zs[70]->1.00062846438525261454098792190059894038`24.,zs[71]->0.99498465623200348109790629853118304194`24.,zs[72]->0.98989080648784945580961488829009779435`24.,zs[73]->0.98534816220783619004877368576729738354`24.,zs[74]->0.98135543976206284633123900312040573392`24.,zs[75]->0.97790927193265276452224372466987674959`24.,zs[76]->0.97500460586932507181245967868916266112`24.,zs[77]->0.97263504902530829842310507202276420718`24.,zs[78]->0.97079316292916908945868071845244288249`24.,zs[79]->0.96947070664871203581308734480692140146`24.,zs[80]->0.968658833174820450508834620003933805`24.,zs[81]->0.96834824281160567241607377551952581701`24.,zs[82]->0.96852929811734085480287625163375141087`24.,zs[83]->0.96919210509995147444624989031135889388`24.,zs[84]->0.9703265653175374266500220950497545025`24.,zs[85]->0.97192240333826097277824117590364370977`24.,zs[86]->0.97396917372901829527778493853756809284`24.,zs[87]->0.97645625140883622205258893478110522519`24.,zs[88]->0.97937280884954599474912056007008119064`24.,zs[89]->0.98270778325239970502968600275841565298`24.,zs[90]->0.98644983648727075592225373369820938983`24.,zs[91]->0.99058731025810194520194860622623936925`24.,zs[92]->0.99510817865789978934497614091973040068`24.,zs[94]->1.0052498695592705748159893720999032856`24.,zs[95]->1.01084437462632089287902235619768099464`24.,zs[96]->1.01676955306827459471064942646876789544`24.,zs[97]->1.02301085639983751860619877039066929415`24.,zs[98]->1.02955311819696835632042171075573671931`24.,zs[99]->1.03638052853133162948414113432340359943`24.,zs[100]->1.04347661496592725900604972553451703816`24.,zs[101]->1.05082423052998571212123422127075234832`24.,zs[102]->1.05840554898353252869056548015969856722`24.,zs[103]->1.06620206758801197160605691537122885169`24.,zs[104]->1.07419461751789733220567449563183611244`24.,zs[105]->1.08236338197795608320264327004182430443`24.,zs[106]->1.09068792203017318243875780473076136723`24.,zs[107]->1.09914721008140143920694165443704710726`24.,zs[108]->1.10771967093552276482926868713477455204`24.,zs[109]->1.11638323027003999040940371046108351609`24.,zs[110]->1.1251153703542968855648340701782349997`24.,zs[111]->1.13389319278273049168135600401186889261`24.,zs[112]->1.14269348794967280373814108369955845536`24.,zs[113]->1.15149281094053462040629822651939519677`24.,zs[114]->1.16026756345645643793952702375361697666`24.,zs[115]->1.16899408132497115143209013382558362387`24.,zs[116]->1.17764872707777720509685967935598420851`24.,zs[117]->1.18620798699891539741230894948101585851`24.,zs[118]->1.1946485719636966692871974694054135288`24.,zs[119]->1.20294752130250851835773238241887427711`24.,zs[120]->1.21108230883658842519393832551131729479`24.,zs[121]->1.21903095014794200325809949123621771435`24.,zs[122]->1.22677211006612033845873700031449487962`24.,zs[123]->1.23428520928409737054027776674887902494`24.,zs[124]->1.24155052895760499362230273888672238891`24.,zs[125]->1.24854931210047562130980737953645532431`24.,zs[126]->1.25526386056600181200653412761156923186`24.,zs[127]->1.26167762640378099096574643915463651581`24.,zs[128]->1.26777529640508512746194086335102309344`24.,zs[129]->1.27354286869884801984255096534927968189`24.,zs[130]->1.27896772033541602432897479908247724522`24.,zs[131]->1.284038664895867978511028589312227556`24.,zs[132]->1.2887459992896454561191208513211248369`24.,zs[133]->1.29308153905019472339565846780201666207`24.,zs[134]->1.29703864160419743660237763998899510239`24.,zs[135]->1.30061221717090201416615238061586833155`24.,zs[136]->1.30379872713948242184265474481427561154`24.,zs[137]->1.30659616997015208151435417881335783315`24.,zs[138]->1.30900405485800111985259281753036466638`24.,zs[139]->1.31102336361450617427333377110338475239`24.,zs[140]->1.3126565013007993543068650110482707826`24.,zs[141]->1.31390723671213208130212464064889394885`24.,zs[142]->1.31478063264232699325282898325172685471`24.,zs[143]->1.31528296984912819983397298367672256073`24.,zs[144]->1.31542165866565266202264590508334972237`24.,zs[145]->1.3152051558927526602114269740673325748`24.,zs[146]->1.31464285456125168429633098037042880502`24.,zs[147]->1.31374501087395415900593851944524795224`24.,zs[148]->1.31252260223629425531056903885455178276`24.,zs[149]->1.310987284499670791344878607171177528`24.,zs[150]->1.30915121752246950784505861632606559761`24.,zp[1]->1.25795782235583664939251378307539554786`24.,zp[2]->1.25673788775674013078321377422236941907`24.,zp[3]->1.25544622038522542591659362853882473767`24.,zp[4]->1.25408049044128985108574748066113907823`24.,zp[5]->1.25263830257276020246581296143136276102`24.,zp[6]->1.25111719545374485408426935181972732688`24.,zp[7]->1.24951464163873026929184255024705148369`24.,zp[8]->1.24782804770414532684420119132992502199`24.,zp[9]->1.24605475472058259573784630808265097325`24.,zp[10]->1.24419203910704306664744911123466302302`24.,zp[11]->1.24223711392514798829394575516118461477`24.,zp[12]->1.24018713067817463518594032669257748385`24.,zp[13]->1.23803918168723450673923492241605039249`24.,zp[14]->1.23579030312495077068688716075853685982`24.,zp[15]->1.23343747879559229935898030354032600825`24.,zp[16]->1.23097764475973852187584791769672963059`24.,zp[17]->1.22840769491110004204273526544497365207`24.,zp[18]->1.22572448762297607224486419112386596225`24.,zp[19]->1.22292485359180460020981655457805871684`24.,zp[20]->1.22000560501509537667888379725904128696`24.,zp[21]->1.21696354625037963082406834850810750604`24.,zp[22]->1.21379548611020410867939050616070275701`24.,zp[23]->1.21049825195504841332493822329958491641`24.,zp[24]->1.20706870575060507720590396193664929531`24.,zp[25]->1.20350376225720147440049128604130264108`24.,zp[26]->1.19980040951612683519044235973067337301`24.,zp[27]->1.19595573178889671896204195266868120532`24.,zp[28]->1.19196693508944490488504213827495554463`24.,zp[29]->1.18783137542404663833740395452883591945`24.,zp[30]->1.18354658981740146812620250613605532239`24.,zp[31]->1.17911033015353020150035106061987392688`24.,zp[32]->1.17452059979469090961334074925523847113`24.,zp[33]->1.16977569285819241257263041605822880539`24.,zp[34]->1.16487423592786928745404271824028078331`24.,zp[35]->1.15981523185275079064921201897448777456`24.,zp[36]->1.15459810513973627894306882130379315354`24.,zp[37]->1.14922274828092971093589987865784638416`24.,zp[38]->1.14368956817267963630496431373743188641`24.,zp[39]->1.13799953158781500268613478608828788306`24.,zp[40]->1.13215420846357895686825680025715611368`24.,zp[41]->1.12615581157726602921363279424107680924`24.,zp[42]->1.12000723101500959713244840930602517156`24.,zp[43]->1.11371206171521344655260627264725154487`24.,zp[44]->1.10727462230777781352523541290678325941`24.,zp[45]->1.10069996349528764485513650070722713026`24.,zp[46]->1.09399386435286243463406358485049016955`24.,zp[47]->1.08716281517495795723130861120541285033`24.,zp[48]->1.08021398587759377019836334315018971072`24.,zp[49]->1.07315517946946209451398617892235270758`24.,zp[50]->1.06599477071746257460915082785773700245`24.,zp[51]->1.05874163081873546050642218822402300456`24.,zp[52]->1.0514050396056556204006998575961019634`24.,zp[53]->1.04399458749560131587056271815859228763`24.,zp[54]->1.03652006999239635598990843140196703264`24.,zp[55]->1.02899137799329096731756695363534953727`24.,zp[56]->1.02141838740790408497540240418178231516`24.,zp[57]->1.01381085162607990350067323946177599784`24.,zp[58]->1.00617830017572515003017163354931083156`24.,zp[59]->0.99852994650919348068668369191253377186`24.,zp[60]->0.99087460728876481452185613364302190245`24.,zp[61]->0.98322063486443903496609485650421403801`24.,zp[62]->0.97557586391377517276838708170993844731`24.,zp[63]->0.96794757250541811646947146676544531316`24.,zp[64]->0.96034245720789819116345522933378529309`24.,zp[65]->0.95276662133194818798286925222144061372`24.,zp[66]->0.94522557499046094008331093315135335879`24.,zp[67]->0.93772424539204375609470910509105785235`24.,zp[68]->0.93026699564531966424959606690369148597`24.,zp[69]->0.92285765032537477035773217669709996274`24.,zp[70]->0.91549952611900170554987662583669549095`24.,zp[71]->0.90819546599747128332018641587629355261`24.,zp[72]->0.90094787554108161938325498718802564878`24.,zp[73]->0.89375876023803623364672564398853809618`24.,zp[74]->0.88662976278447078147811311840377107568`24.,zp[75]->0.87956219961003812067247633865780065636`24.,zp[76]->0.87255709603571453117203808567892512694`24.,zp[77]->0.86561521963221219893708997405833107783`24.,zp[78]->0.85873711148617000753247416804013030098`24.,zp[79]->0.851923115196828717649035289943344933`24.,zp[80]->0.84517340351930159326571734789685984418`24.,zp[81]->0.83848800264384813374268288551869704973`24.,zp[82]->0.83186681415623644674124198430494692589`24.,zp[83]->0.82530963476498293353355978459248253208`24.,zp[84]->0.81881617390956918698728092089502554269`24.,zp[85]->0.81238606938204206647509175988375034485`24.,zp[86]->0.80601890110482100896531822268805099756`24.,zp[87]->0.79971420321188516510392138451001066716`24.,zp[88]->0.7934714745803121216265777077688213406`24.,zp[89]->0.7872901879556310023968672740825339509`24.,zp[90]->0.78116979780862733865190068690835140499`24.,zp[91]->0.77510974705387139140829772155421809449`24.,zp[92]->0.76910947275193062171132366573936624017`24.,zp[93]->0.76316841090841730204531230349972579372`24.,zp[94]->0.75728600047404140675945346608367776224`24.,zp[95]->0.75146168664092114166484314330736554796`24.,zp[96]->0.74569492352171142203831720072311785955`24.,zp[97]->0.73998517628975333281775499945326723419`24.,zp[98]->0.73433192285049508171460876981302899439`24.,zp[99]->0.7287346551069295560624877104855891588`24.,zp[100]->0.72319287987475889667896271793338919658`24.,zp[101]->0.71770611949644386080378231794773006805`24.,zp[102]->0.71227391219722898492464362932965405413`24.,zp[103]->0.70689581222065313786978542513948652964`24.,zp[104]->0.70157138977595602235017971107388955772`24.,zp[105]->0.69630023082517017931250203744966384869`24.,zp[106]->0.69108193673353961258946449884730056855`24.,zp[107]->0.6859161238032235358307401097766661822`24.,zp[108]->0.68080242270701843205776797065873927296`24.,zp[109]->0.67574047783605268100031330516084547756`24.,zp[110]->0.6707299465730614528531270023413552547`24.,zp[111]->0.66577049850091773907133295102380577245`24.,zp[112]->0.66086181455455656039588635397239578561`24.,zp[113]->0.65600358612325752364391033895140184561`24.,zp[114]->0.65119551410941565892468942127088437234`24.,zp[115]->0.64643730794939747460065582166095200143`24.,zp[116]->0.64172868460181046064541696205752093208`24.,zp[117]->0.63706936750846899917608773481521810378`24.,zp[118]->0.63245908553347490369799987747260764913`24.,zp[119]->0.6278975718861026312277940271986070093`24.,zp[120]->0.62338456303354359474126308683752827961`24.,zp[121]->0.61891979760997800273763357099273669899`24.,zp[122]->0.61450301532886543892440903088118378697`24.,zp[123]->0.61013395590573925296628763914221054774`24.,zp[124]->0.6058123579991210596133722468432899457`24.,zp[125]->0.60153795817741131109618600456774620332`24.,zp[126]->0.59731048991973647002941958009115366013`24.,zp[127]->0.59312968265872502815517718104248234478`24.,zp[128]->0.58899526087303179534338743995345852069`24.,zp[129]->0.58490694323712691607456793699611486032`24.,zp[130]->0.58086444183541328956819124246906701593`24.,zp[131]->0.57686746144713941779000168246694000191`24.,zp[132]->0.57291569890784524180002867011956199689`24.,zp[133]->0.5690088425522317877853186606162631007`24.,zp[134]->0.56514657174240068767198529897622165845`24.,zp[135]->0.56132855648438899746562597688727243803`24.,zp[136]->0.55755445713485236495733271063267194534`24.,zp[137]->0.55382392419865030396903459887311872446`24.,zp[138]->0.55013659821698913535984041813895324781`24.,zp[139]->0.54649210974468956865045560595119002149`24.,zp[140]->0.54289007941417166883042550450088141581`24.,zp[141]->0.53933011808255490249173772322699233709`24.,zp[142]->0.53581182705823729464949255816529410038`24.,zp[143]->0.53233479839980909509470815975494020647`24.,zp[144]->0.52889861528748844785872092782196388425`24.,zp[145]->0.5255028524453548770280418518960106074`24.,zp[146]->0.52214707664485528510426284520012641123`24.,zp[147]->0.5188308471993315294651788288004388574`24.,zp[148]->0.51555371661634187870061163753824970353`24.,zp[149]->0.51231523105375050427732287754340297774`24.,zp[150]->0.50911493119447451666912000144677833317`24.,\[Eta]->0.23590572598007892234469142698930640866`24.};


roots[n_] := Sort[Select[x/.Solve[0==HermiteH[n,x]],#>=0&]];
weights[n_] := Table[N[{roots[n][[i]],2^(n-1) n! Sqrt[Pi] /(n^2 HermiteH[n-1,roots[n][[i]]]^2)},60],{i,1,(n+1)/2}];


roots[2][[1]]


weights[4]


Map[(GHTable[#] = weights[#])&, Range[40]];


GHIntegral[f_, n_] := Sum[f[GHTable[n][[i,1]]]*GHTable[n][[i,2]],{i,1,(n+1)/2}]


GHIntegral[1/(1+#^2)&,40]


NIntegrate[Exp[-x^2] 1/(1+x^2), {x,0,Infinity}, AccuracyGoal->Infinity, WorkingPrecision->30]


Table[Total[GHTable[n][[All,2]]]-Mod[n,2]*GHTable[n][[1,2]]/2,{n,1,40}]


Show[ListLogPlot[{weights[10],weights[20],weights[30],weights[40]}],LogPlot[Exp[-x^2],{x,0,8}]]


dd = 200/100;
{aa,rr} = {2,144/10};
myfp = test;


myintegrands = integrandsListIsotropic[[1;;2,3]];


myintegrands = (4 t^3 myintegrands/.{y->t^4})/.t->y;
Unprotect[IntegralConfigurations];
IntegralConfigurations = {{0,1/16,30},{1/16,9/10,40},{9/10,3/2,30}};


sparse = 20;
integrands1 = Flatten[Table[NumericDerivatives[myintegrands/.{n->nn, \[Rho] -> (1+ sparse (i-1)) eps}/.regulatorReplacement]//.ConstRep[dd,rr,aa]/.test,{i,1,gridSize/sparse}]];
integrands[q_] = integrands1 /.y->q;


exact = NIntegrate[integrands[q],{q,0,4},WorkingPrecision->40, AccuracyGoal->16, PrecisionGoal->16, MaxRecursion->50];


Diff[f_] := Abs[GetIntegral[f]-NIntegrate[f[y],{y,0,3},WorkingPrecision->40, AccuracyGoal->16, PrecisionGoal->16, MaxRecursion->50]]


DiffExpand[f_, config_]:= Block[{exact,glint, int2},
	exact = NIntegrate[f[q],{q,config[[1]],config[[2]]},WorkingPrecision->40, AccuracyGoal->16, PrecisionGoal->16, MaxRecursion->50];
	glint = GLIntegral[f,config[[1]],config[[2]],config[[3]]];
	Return[Max[Abs[exact-glint]]];
]



IntegralConfigurations = {{0,1/4,25},{1/4,125/100,60},{125/100,3,20}};


ListLogPlot[Quiet[Diff[integrands], NIntegrate::precw]]


ListLogPlot[Quiet[Diff[integrands], NIntegrate::precw]]


err[d_]:= Total[Map[DiffExpand[threshold4[d,6,-195/100], #]&,ints]]


best = <|4/5->0.0000103561101300632411288218053145386096314391313445784708`28.431471862123985,9/10->9.588475029020712707064213358135158221841295957662651`27.486279512886448*^-7,1->2.51195806844213730184127922758660610583273508348`22.986872318893717*^-11,11/10->4.1386846073600888814229235702362096471935508001146`25.281179117700866*^-9,6/5->2.145353095046127199917830926067274863986441783785`24.06917959509147*^-10,13/10->4.29541348324581177026605515143137493630776753603`23.440552630726735*^-11,7/5->3.09099428319570055196777189990873141956775008976`23.364481068726867*^-11,3/2->2.9477213713302202588809811528280442992061937095`23.408044095939434*^-11,8/5->2.99365291739684294211375890033001839082245856591`23.47657457847854*^-11,17/10->3.0225157403675361596614795562020087427728995594`23.54044565911666*^-11,9/5->3.03960630032711346929116878599809784692932117481`23.60069350511252*^-11,19/10->3.04391131645738608686894186394678673469267267508`23.657376627554186*^-11,2->3.0358442798949954913457019572096805137668959311`23.71071041154204*^-11,21/10->3.01594598610105017380154505898418839721185137734`23.7608861896504*^-11,11/5->2.98481174146823161536706942837744154221475998911`23.80806727777779*^-1|>;
bestInts = {};
previous = <||>;


a1 = 2/32;
a2 = 90/100;
ints = {{0,a1,30},{a1,a2,40},{a2,15/10,30}};
previous = errdata;
errdata = Association[Map[(#->err[#])&,Range[20,60]/20]];
If[Mean[Values[errdata]]<Mean[Values[best]],
best = errdata;
bestInts = ints;
]


ListLogPlot[{errdata, best,previous}, PlotRange->All]


ListLogPlot[errdata]


IntegralConfigurations = {{0,1/10,20},{1/10,2,40},{2,4,20}};


Unprotect[IntegralConfigurations];
IntegralConfigurations = {{0,2/3,40},{2/3,12/10,40},{12/10,3,30}};


data = Association[Table[({d,m}/.{d->i/5, m->-1+(-j)/8}) -> (Diff[threshold[d,5,m]/.{d->i/5, m->-1+(-j)/8}]),{i,3,15},{j,0,7}]];


Get["FixedPoint.wl"]


data = Association[Table[({d,m}/.{d->i/5, m->-1+(-j)/8}) -> (Diff[threshold3[d,3,m]/.{d->i/5, m->-1+(-j)/8}]),{i,3,15},{j,0,7}]];


data = Association[Table[({d,m}/.{d->i/5, m->-1+(-j)/8}) -> (Diff[threshold4[d,3,m]/.{d->i/5, m->-1+(-j)/8}]),{i,3,15},{j,0,7}]];


nested = TransposeAssociation[NestAssociation[data]];


Show[ListLogPlot[Values[nested],PlotLegends->N[Keys[nested]],PlotRange->All], LogPlot[10^-8,{x,0,3}]]


Show[ListLogPlot[Values[nested],PlotLegends->N[Keys[nested]],PlotRange->All], LogPlot[10^-8,{x,0,3}]]


(* ::Subsection::Closed:: *)
(*afd*)


PMS[data_,func_:Re] := Block[{norm=Normal[data], pairs, inter,pms},
pairs =  Map[{norm[[#,1]],func[norm[[#,2]]]}&,Range[Length[norm]]];

inter = Interpolation[pairs];
pms = FindRoot[inter'[\[Alpha]]==0,{\[Alpha],2}];

Return[{\[Alpha], inter[\[Alpha]]}/.pms]];


SparsifyAssociation[association_,numPoints_] := Block[{ass,ratio, positions},
	ass = KeySort[association, #1>=#2&];
	If[numPoints >= Length[ass], 
		Return[ass];];
	
	ratio = Length[ass]/numPoints;
	positions = Floor[Range[0,numPoints-1]*ratio]+1;
	
	Return[Association[Map[Keys[ass][[#]]->ass[[#]]&,positions]]];
];



Remove[roots, weights]


{fps, exps}=Import["/home/andrzej/FunctionalRenormalization/ON_functional_smooth/fixed_points.wdx"];


partnest = OneLevelNest[fps];
vdata = NestAssociation[Association[Map[(Keys[partnest][[#]]->1+(v[0]/.partnest[[#,1]]) / Keys[partnest][[#,1]])&,Range[Length[partnest]]]]];



vdata = NestAssociation[Association[Map[(Keys[partnest][[#]]->v[150]/.partnest[[#,1]])&,Range[Length[partnest]]]]];


ready = vdata[15/10];


Export["~/Desktop/v0.png", ListLogPlot[Values[ready], LabelStyle->20, PlotLegends->PointLegend[N[Keys[ready]], LegendLabel->"N"], 
	GridLines->Automatic, AxesLabel->{ "1+V[0]/\[Alpha]","d"}, ImageSize->Large]]


Export["~/Desktop/vmax.png", ListLogPlot[Values[ready], LabelStyle->20, PlotLegends->PointLegend[N[Keys[ready]], LegendLabel->"N"], 
	GridLines->Automatic, AxesLabel->{"V[\[Rho]max]", "d"}, ImageSize->Large]]


sel = KeySelect[fps,MatchQ[{25/10,10/10,16/10, r_}]];
{aa, nn, dd, rr} = Keys[sel][[1]];
myfp = sel[[1]];
myfp[[1;;-2,2]] /= aa;


PlotFixedPoint[myfp,rr]


DimShift[ass_]:=Association[Map[Keys[ass][[#]]->ass[[#]]+Keys[ass][[#]]-2&,Range[Length[ass]]]];


shifted = Map[DimShift,ready];


test = KeySelect[fps, MatchQ[{2,19/10,d_,r_}]][[-11]];


Keys[KeySelect[fps, MatchQ[{2,19/10,d_,r_}]]][[-11]]


ListPlot[Values[shifted], PlotLegends->N[Keys[shifted]]]


data = NestAssociation[exps];


Map[Length,data[[1]]]


exp = e1;
aa = 2;
xaxis = d;
title = ToString[exp]<>"("<>ToString[xaxis]<>",a="<>ToString[N[aa]]<>").png";
sel = Re[data[exp,aa]];
p = ListPlot[Values[sel], LabelStyle->Directive[20], PlotLegends->PointLegend[N[Keys[sel]], LegendLabel->"N"], AxesLabel->{exp, xaxis}, GridLines->Automatic, ImageSize->Large];
Export["~/Desktop/"<>title,p];


Show[PlotAssociation[data[\[Eta],25/10]],Plot[2-d,{d,1,2}]]


PlotAssociation[N[TransposeAssociation[TransposeAssociation[data[e2]][2]][[1;;-1;;4]]]]


PlotAssociation[data[e1,20/10]]


PlotAssociation[data[\[Eta],10/10]]


ListPlot[Re[Values[sel]], PlotLegends->Keys[sel],PlotRange->All]


ListPlot[Re[Values[sel]], PlotLegends->Keys[sel],PlotRange->All]


{aa,nn,dd,\[Rho]m} = {2.,1.4,1.9,40};


Find\[Eta][integrands_,guess_,consts_]:=Block[{eq,dropped, int},
int[q_] =integrands[[2,4]]/.regulatorReplacement /. y->q;
eq = NumericDerivatives[(GetIntegral[int]+integrands[[2,2]])/.\[Rho]->0];
dropped = guess[[1;;-2]];
Return[Solve[0==eq/.dropped/.consts,{\[Eta]}]];
]


\[Rho]m=6;


gridSize = 60; zNormalization=0;


guessv = Table[v[i]->-1.9+10^-6x^10/.x->i*\[Rho]m/gridSize,{i,0,gridSize}];
guesszs = Table[zs[i]->1+x/(\[Rho]m)-x^2/\[Rho]m^2/2 + 7. Exp[-(x-\[Rho]m 7/12)^2/(1/6 \[Rho]m)]/.x->i*\[Rho]m/gridSize,{i,1,gridSize}];
guesszp = Table[zp[i]->1.-x/(2\[Rho]m) /.x->i*\[Rho]m/gridSize,{i,1,gridSize}];
guess = Join[guessv,guesszs,guesszp,{\[Eta]->0.2}];


reorganized = KeyMap[{#[[1]],#[[3]],#[[4]],#[[2]]}&,KeySelect[exps,MatchQ[{e_, a_, nn_,dd_}/;MemberQ[{e1,e2,\[Eta]},e]]]];
partNested = OneLevelNest[reorganized];


pmsrawdata = Quiet[Map[PMS[#][[2]]&,partNested]];


pmsdata = NestAssociation[pmsrawdata];


PlotAssociation[pmsdata[\[Eta]]]


PMS[partNested[[1]]]


partNested


Keys[reorganized]


e1Exps = NestAssociation[KeySelect[exps, MatchQ[{e1, a_, nn_,dd_}]]];
etaExps = NestAssociation[KeySelect[exps, MatchQ[{\[Eta], a_, nn_,dd_}]]];


exps = Join[e1Exps,etaExps];


PlotANDep[exp_,a_]:=Block[{range,data, p},
range = <|\[Eta]->{0,1},e1->{0,3}|>[exp];
data = exps[exp,a];
p = ListPlot[Values[data], PlotLegends->PointLegend[N[Keys[data]], LegendLabel->"N"], PlotRange->{{1,2.2},range}, AxesLabel->{d,exp}, GridLines->Automatic];
Export["~/Desktop/"<>ToString[exp]<>"(d,a="<>ToString[N[a]]<>").png",p];
];


exps[e1,2]


PlotANDep[\[Eta],4]


CritDimension[data_] := Block[{sel, nested, criticalDimensions, CheckDimension},
	sel = OneLevelNest[KeySelect[data, MatchQ[{e1, aa_, nn_,dd_}]]];
	nested = NestAssociation[Map[Min[Keys[#]]&,sel]][e1];
	criticalDimensions = <| |>;
	CheckDimension[newData_] := Block[{},
		Map[If[!MemberQ[Keys[criticalDimensions],#] || criticalDimensions[#]<newData[#], 
			criticalDimensions[#]=newData[#]]&, Keys[newData]];
	];
	Map[CheckDimension,Values[nested]];
	
	Return[criticalDimensions];
];


{fps, exps}=Import["/home/andrzej/FunctionalRenormalization/ON_functional_smooth/fixed_points.wdx"];


etadata = KeySelect[exps, #[[1]]==\[Eta]&];


nested = NestAssociation[etadata][\[Eta]];


DimensionShift[ass_]:= Block[{keys, vals, newvals},
keys = Keys[ass];
vals = Values[ass];
newvals = keys-2+vals;
Return[Association[Map[#[[1]]->#[[2]]&,Transpose[{keys,newvals}]]]]];


plotData = Map[DimensionShift,nested[2]];


ListLinePlot[Values[nested[2]], PlotLegends->N[Keys[nested[2]]]]


ListLinePlot[Values[plotData], PlotLegends->N[Keys[plotData]]]


Keys[nested[2]]


Length[etadata]


Keys[exps]


n2exps = NestAssociation[KeySelect[exps, MatchQ[{e_,a_,n_,d_}/; n>=2]]];


ListPlot[Values[n2exps[e1,1]],PlotLegends->Keys[n2exps[e1,1]]]


dimensions = CritDimension[exps];
myPoints = Transpose[{Values[dimensions],Keys[dimensions]}];
lauPoints = 1+{{0,0},{63,3},{105,13},{130,25},{159,50},{180,75},{197,100},{211,125},{222,150},{231,175},{239,200},{250,250}}/250.;


sel = NestAssociation[KeySelect[exps, MatchQ[{e1, 2, nn_, dd_}]]][e1,2];
dimensions2 = Normal[Map[Min[Keys[#]]&,sel]];
myPoints2 = Transpose[{Values[dimensions2],Keys[dimensions2]}];


p = Show[ListPlot[{myPoints, myPoints2,lauPoints}, PlotRange->{{0.95,2.05},{0.95,2.05}}, PlotLegends->{"My", "My(\[Alpha]=2)", "Lau"}],
Plot[2+ (d-2)*Pi^2/4,{d,0.95,2.05},PlotLegends->LineLegend[{"Cardy"},LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 20]],PlotStyle->Red],
GridLines->Automatic, AxesLabel->{"d","N"}]
Export["~/Desktop/Lauline.png",p];


PlotADep[exp_,nval_]:= Block[{selected, lowd, data1,data2,data,p},
selected = NestAssociation[KeySelect[exps, MatchQ[{exp, aa_, nval, dd_}/; 1.5<=dd<=3]]];
lowd = NestAssociation[KeySelect[exps, MatchQ[{exp, aa_, nval, dd_}/; dd<1.5]]];

data1 = SparsifyAssociation[UnnestAssociation[TransposeAssociation[selected[exp]][nval]],7*50];
data2 = SparsifyAssociation[UnnestAssociation[TransposeAssociation[lowd[exp]][nval]],7*50];
data = NestAssociation[Join[data1,data2]];
p = ListPlot[Values[data],PlotLegends->PointLegend[Keys[data], LegendLabel->"\[Alpha]"], AxesLabel-> {"d","e1"}];
Export["~/Desktop/"<>ToString[exp]<>"(a).png",p];
Return[p];
]


(* ::Subsection:: *)
(*Tricritical Fixed Points*)


{fps,exps} = Import["/home/andrzej/FunctionalRenormalization/tricritical_functional_smooth/fixed_points.wdx"];


nestedFps = NestAssociation[fps];
nestedExps = NestAssociation[exps];


nestedExps


PlotFixedPoint[nestedFps[4,2][[-1,1]],1]


data = TransposeAssociation[nestedExps[e2]][2];


p1 = ListPlot[Values[data],PlotLegends->PointLegend[N[Keys[data]],LegendLabel->"\[Alpha]"],AxesLabel->{d,e3},LabelStyle->Directive[Bold, Medium],ImageSize->Large, GridLines->Automatic];


data2 = TransposeAssociation[data][[1;;-6]];
p2 = ListPlot[Values[data2],PlotLegends->PointLegend[N[Keys[data2]],LegendLabel->"d"],AxesLabel->{\[Alpha],e3},LabelStyle->Directive[Bold, Medium],ImageSize->Large, GridLines->Automatic,PlotRange->All];


Export["~/Desktop/e3(d).png",p1];
Export["~/Desktop/e3(a).png",p2];


p2


PlotAssociation[TransposeAssociation[data]]


ListPlot[nestedExps[e3,4]]


Keys[nestedFps]


Keys[fps]


nestedFps



data = Import["~/Desktop/action.csv"];
data = data[[2;;122]];
funcs = {v, zs, zp};
Unprotect[gridSize, zNormalization];
gridSize = Length[data]-1;
guess = Flatten[Table[funcs[[i]][j-1]->data[[j,i]],{i,1,3},{j,1,gridSize+1}]];
znormdata = Select[guess, MatchQ[zs[zNormalization] -> v_],1];
normFactor = znormdata[[1,2]];
guess[[All,2]] /= normFactor;
guess = DeleteCases[guess, zs[zNormalization]->v_];
guess = DeleteCases[guess, zp[0]->v_];
guess = Append[guess, \[Eta]->0.2];


Get["FixedPoint.wl"]


integrands = integrandsListIsotropic/.n->14/10;


ChangeVariable[integrands_] := Block[{newIntegrands = integrands, VarChange},
	VarChange[i_] := (4 t^3 i /.{y->t^4})/.t->y;
	newIntegrands[[3]] = VarChange[newIntegrands[[3]]];
	If[Length[integrands] == 4,
		newIntegrands[[4]] = VarChange[newIntegrands[[4]]];
	];
	Return[newIntegrands];
];
integrands = Map[ChangeVariable, integrands];

Unprotect[IntegralConfigurations];
IntegralConfigurations = {{0,1/16,30},{1/16,9/10,40},{9/10,3/2,30}};
Protect[IntegralConfigurations];


eqs = IntegrateEquations[integrands];


fp = FindFixedPoint[eqs, guess, ConstRep[19/10, 7/10, 12/10]];


PlotFixedPoint[fp, 7/10]


fp


MatchQ[zp[20]->.5,zp[i_] -> v_ /; v<=1]


DeleteCases[guess, zs[zNormalization]->1]


Length[]


guess = 


(* ::Subsection:: *)
(*Archive*)


Quit[];


(* ::Subsubsection::Closed:: *)
(*PlotKappa*)


(* ::Text:: *)
(*Plots FP characteristics at potential minimum*)


{fps,exps} = Import["/home/andrzej/FunctionalRenormalization/ON_functional_smooth backup/fixed_points.wdx"];
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
dir = "~/Desktop/";
If[!DirectoryQ[dir],CreateDirectory[dir]];
plot[1,U''[Subscript[\[Rho],0]],"u"];
plot[2,Subscript[Z,\[Sigma]][Subscript[\[Rho],0]],"sigma"];
plot[3,Subscript[Z,\[Pi]][Subscript[\[Rho],0]],"pi"];
plot[4,U'[0],"v"];
plot[5,U''[0],"u2"];
];

PlotKappa[KeySelect[nested,MemberQ[{4.,3.,2.5,2.25,2.},#]&]]


(* ::Subsubsection:: *)
(*Plot potentials*)


(* ::Text:: *)
(*Plots full shapes of fixed points*)


{fps,exps} = Import["/home/andrzej/FunctionalRenormalization/ON_functional_smooth backup/fixed_points.wdx"];
nested = NestAssociation[fps];
nested = KeySelect[nested,#>=2&];

PlotPotentials[nn_]:=Block[{len,ds,lim=4,\[Rho]Maxs,fixedpoints,vs,\[Pi]s,\[Sigma]s,plot,dir, grid},
grid = Length[nested[[1,1,1]]]/3-1;
len = Length[nested[nn]]-6;
ds = Table[Keys[nested[nn]][[Floor[i*len/(lim+1)]+1]],{i,0,lim+1}];
ds = Drop[ds,{2}];
If[nn==2.5,ds = {3.`,2.5000000000000018`,2.200000000000003`,2.0500000000000025`}; lim=3;];
\[Rho]Maxs = Flatten[Map[Keys[nested[nn,#]]&,ds]];
fixedpoints = Map[nested[nn,#][[1]]&,ds];
vs = Table[Interpolation[Table[{i \[Rho]Maxs[[j]]/grid, v[i]/.fixedpoints[[j]]},{i,0,grid}]],{j,1,lim+1}];
\[Sigma]s = Table[Interpolation[Table[{i \[Rho]Maxs[[j]]/grid, z[i]/.fixedpoints[[j]]/.z[0]->1},{i,0,grid}]],{j,1,lim+1}];
\[Pi]s = Table[Interpolation[Table[{i \[Rho]Maxs[[j]]/grid, (z[i]-2 i \[Rho]Maxs[[j]]/grid yy[i])/.fixedpoints[[j]]/.z[0]->1},{i,0,grid}]],{j,1,lim+1}];

l1 = .025;
l2 = .001;
s = l1/2;
dash[1] = {};
dash[2] = {l1,s};
dash[3] = {l1,s,l2,s};
dash[4] = {l1, s, l2, s, l2, s};
dash[5] = {l2,s};
dashing = {Null, Dashed, DotDashed, Dotted,Null};
curves = Table[Directive[AbsoluteThickness[14],ColorData[97, "ColorList"][[i]]],{i,1,5}];
markers = Tuples[{{"\[FilledCircle]", "\[FilledSquare]", "\[FilledDiamond]", "\[FilledUpTriangle]",""},{30}}];
plot[fs_,title_,label_]:=Block[{p,funcs,range=Association[]},
funcs = MapThread[#2[t 0.845 #1]&,{\[Rho]Maxs,fs}];
range[Subscript[Z,\[Sigma]]] = {0.75,2.5};
range[Subscript[Z,\[Pi]]] = {0.6,1.1};
range["U'[\[Rho]]"] = {-2.05,2.05};

p = Plot[funcs,{t,0,1.15}, 
	PlotLegends->Placed[LineLegend[Map["d="<>ToString[PaddedForm[#,{3,2}]]&,ds], LegendFunction->"Frame", LegendLayout->"Row",
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 50, Black], LegendMarkerSize->{30,40}],Below],  
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 50],
	AxesLabel->{Style["\[Rho]",FontSize->40,Black,Bold],Style[title,FontSize->50,Black,Bold]}, 
	ImageSize->Full, PlotRange->{{0,1.15},range[title]}, PlotStyle->curves,PlotTheme->"Monochrome",
	Ticks->Automatic,GridLines->Automatic];
Export[dir<>label<>".png",p];
];
dir = "~/Desktop/"<>ToString[ToString[PaddedForm[nn,{3,2}]]]<>"/";
If[!DirectoryQ[dir],CreateDirectory[dir]];
plot[vs,"U'[\[Rho]]","v"];
plot[\[Sigma]s,Subscript[Z,\[Sigma]],"sigma"];
plot[\[Pi]s,Subscript[Z,\[Pi]],"pi"];
];


PlotPotentials[2.5]
