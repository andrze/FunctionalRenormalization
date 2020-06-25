(* ::Package:: *)

(* ::Text:: *)
(*Imports ModuleNaming module and FixedPoint library*)


gridSize = 120;


arguments = $CommandLine;
notebookPath = arguments[[3]];
slashes = StringPosition[notebookPath,"/"];
If[Length[slashes]!=0,
	projectPath = Strin \[AliasDelimiter]gTake[notebookPath, slashes[[-1,1]]];
	SetDirectory[projectPath];
];
Get["ModuleNaming.wl"];
Get["FixedPoint.wl"];
Get["Plots.wl"];


(* ::Text:: *)
(*An initializer for debugging process. To debug change mode to debug*)


mode = "run";
If[ mode == "debug",
SetDirectory["/home/andrzej/Studia/Anizotropie/FunctionalRenormalization"];
Get["ModuleNaming.wl"];
Get["FixedPoint.wl"];
Get["Plots.wl"];
task = "fixedPointDict";
anisotropySymmetry = 4;(*"ON";*)
regulatorName = "smooth";
Print["Debug mode is active!"]
]


(* ::Text:: *)
(*Serializes the command line arguments*)


task = "none";
anisotropySymmetry = 4;
regulatorName = Null;
configuration = Null;
Block[{i},
For[i = 1, i<=Length[arguments], i++,
	If[arguments[[i]] == "Zp", 
		If[arguments[[i+1]] == "ON",
		anisotropySymmetry = "ON";
		ON = True,
		anisotropySymmetry = ToExpression[arguments[[i+1]]]],
		ON = False;];
	If[arguments[[i]] == "task", task = arguments[[i+1]]];
	If[arguments[[i]] == "regulator", regulatorName = arguments[[i+1]]];
	If[arguments[[i]] == "N", nn = ToExpression[arguments[[i+1]]]];
	If[arguments[[i]] == "conf", configuration = ToExpression[arguments[[i+1]]]];
	];
];
If[anisotropySymmetry=="ON",
runconf = ONRUNCONF,
runconf = Z4RUNCONF];


If[anisotropySymmetry == Null || regulatorName == Null || (anisotropySymmetry == "ON" && nn == Null),
	Print["Specify all required options: task, regulator"];
	Quit[];
	];


SelectRegulator[regulatorName];
cacheDirectory = GetCacheDirectory[anisotropySymmetry, 0, regulatorName];


(* ::Text:: *)
(*Imports cached flow equations*)


isotropicEquationsFile = GetEquationsFileName[anisotropySymmetry, 0, 0];
If[FileExistsQ[isotropicEquationsFile], 
	integrandsListIsotropic = Import[isotropicEquationsFile],
	Print["No isotropic equations found"]; Quit[]];
If[ON,
integrandsListAnisotropic = {},
anisotropicEquationsFile = GetEquationsFileName[anisotropySymmetry, 0, 1];
If[FileExistsQ[anisotropicEquationsFile], 
	integrandsListAnisotropic = Import[anisotropicEquationsFile],
	Print["No anisotropic equations found"]; Quit[]];
];


guess3d = Import[INITIALGUESSFILE];
{guess3d, zNormalization} = AdjustNumberOfPoints[guess3d,gridSize];

{keyLabels, configurations} = Import[runconf];

If[!Element[configuration,Integers] || configuration<1 || configuration > Length[configurations],
	Print["Choose proper run cofiguration id"];
	Quit[]
	];
{keys, constants} = configurations[[configuration]];

If[MemberQ[Keys[constants], "N"], 
	integrandsListIsotropic = (integrandsListIsotropic /. n -> constants["N"]);
	integrandsListAnisotropic = (integrandsListAnisotropic /. n -> constants["N"])];

constants = Drop[constants, Position[Keys[constants],"N"][[1]]];

eqlist = IntegrateEquations[integrandsListIsotropic];

{newGuess, new\[Rho]Max} = \[Rho]MaxScan[eqlist, guess3d, 3., constants["\[Rho]Max"], 2., 0.1*constants["\[Rho]Max"], 20];

fixedPointDict = FindFixedPointDict[eqlist, newGuess, keys, keyLabels, constants];
	
AppendTo[keyLabels,"\[Rho]Max"];
constants = Drop[constants, Position[Keys[constants],"\[Rho]Max"][[1]]];
exponents = FindExponents[integrandsListIsotropic, integrandsListAnisotropic, fixedPointDict, keyLabels, constants];

If[ON,
	nestedFPD = NestAssociation[fixedPointDict];
	fixedPointDict = UnpackAssociation[Association[constants[["N"]]->nestedFPD]];]

If[ON,
	exponents = UnpackAssociation[TransposeAssociation[Association[constants[["N"]]->exponents]]],
	exponents = UnpackAssociation[exponents]]


Export[cacheDirectory<>FIXEDPOINTSFILES[configuration],{fixedPointDict,exponents}];
	
If[FileExistsQ[cacheDirectory<>FIXEDPOINTSFILE],
	{allFPs, allExponents} = Import[cacheDirectory<>FIXEDPOINTSFILE];
	allFPs = Join[allFPs, fixedPointDict];
	allExponents = Join[allExponents, exponents],
	
	allFPs = fixedPointDict;
	allExponents = exponents];

If[ON,
	plotExp = NestAssociation[FindPMSValues[allExponents]];
	For[i=1, i<=Length[plotExp], i++,
		PlotDataSeries[plotExp[[i]], {"d", Keys[plotExp][[i]]}, "N", cacheDirectory];
	];]



(*selector = Table[Divisible[i-1,4],{i,1,Length[exp[\[Nu]][2.]]}];
data = Pick[exp[\[Nu]][2.],selector];
ListPlot[Values[data], PlotLegends->PointLegend[Keys[data], LegendLabel\[Rule]"d", LegendFunction->"Panel"]]*)


(*zNormalization=40;
{newfp, zNormalization} = AdjustNumberOfPoints[fps[[1]],60];
eq = integrandsListIsotropic/.n\[Rule]2;
eqlist = IntegrateEquations[eq];
newfp = FindFixedPoint[eqlist, newfp, ConstRep[3.,0.12,2]];
PlotFixedPoint[newfp]
sm3 = StabilityMatrix[eq,newfp,ConstRep[3.,0.12,2]];*)


If[task == "plot",
	{keys, keyLabels, constants} = Import[ONRUNCONF];
	
	allExponents = Association[];
	For[configuration=1, configuration<=Length[constants], configuration++,
		name = cacheDirectory<>FIXEDPOINTSFILES[configuration];
		
		If[FileExistsQ[name],
			{fixedPoints, exponents, matrices} = Import[name];
			allExponents[constants[[configuration]]["N"]] = exponents;
		];
	];
	transposed = TransposeAssociation[allExponents];
	For[i=1, i<=Length[transposed], i++,
	
		PlotDataSeries[transposed[[i]], {"d", Keys[transposed][[i]]}, "N", cacheDirectory];
		];
];


(*allFPs = Association[];
allExponents = Association[];

For[conf=1, conf\[LessEqual] Length[FileNames[All,cacheDirectory]],conf++,
	If[FileExistsQ[cacheDirectory<>FIXEDPOINTSFILES[conf]],
		Quiet[Check[{fixedPointDict, exponents} = Import[cacheDirectory<>FIXEDPOINTSFILES[conf]],
		Print[conf];
		Continue[]]];
		
		exponents = UnpackAssociation[exponents];
		allFPs = Join[allFPs, fixedPointDict];
		allExponents = Join[allExponents, exponents]];
]
Export[cacheDirectory<>FIXEDPOINTSFILE, {allFPs, allExponents}];
	

allExponents = NestAssociation[allExponents];*)


(*allFPs = Association[];
allExponents = Association[];

For[conf=1, conf\[LessEqual] Length[FileNames[All,cacheDirectory]],conf++,
	If[FileExistsQ[cacheDirectory<>FIXEDPOINTSFILES[conf]],
		Quiet[Check[{fixedPointDict, exponents} = Import[cacheDirectory<>FIXEDPOINTSFILES[conf]],
		Print[conf];
		Continue[]]];
		
		allFPs = Join[allFPs, fixedPointDict];
		allExponents = Join[allExponents, exponents]];
]	

allExponents = NestAssociation[allExponents];*)


(*exp = \[Eta];
data = TransposeAssociation[allExponents[exp]][2.01];
ListPlot[Values[data],PlotLegends->PointLegend[Round[Keys[data],0.01], LegendLabel\[Rule]"N", LegendFunction->"Panel",
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 20]],LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 16],
	AxesLabel\[Rule]{\[Alpha],exp},PlotRange\[Rule]{{0.85,6.1},{0.,0.4}}, ImageSize\[Rule]Large]*)


(*exp = \[Nu];
selector = Table[Divisible[i-1,4],{i,1,Length[allExponents[exp]]}];
selector[[-1]]=False;
data = Pick[allExponents[exp],selector];
p=ListPlot[Values[data], PlotLegends->PointLegend[Round[Keys[data],0.01], LegendLabel\[Rule]"d", LegendFunction->"Panel",
	LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 20]],LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 16],
	AxesLabel\[Rule]{\[Alpha],exp},PlotRange\[Rule]{{0.85,6.1},{-0.15,0.4}}, ImageSize\[Rule]Large];
Export[cacheDirectory<>"PMS/"<>ToString[exp]<>".png",p];*)
