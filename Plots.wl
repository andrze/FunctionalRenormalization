(* ::Package:: *)

Plot\[Rho]Dependence[fixedPointDict_, cacheDirectory_]:= Block[
	{i, graphDirectory, dimensions, \[Rho]Dependence, \[Rho]Points, graphPoints, graph, dim},
	graphDirectory = cacheDirectory <> "rho_optimization_graphs/";
	
	If[!DirectoryQ[graphDirectory], CreateDirectory[graphDirectory]];
	
	dimensions = Keys[fixedPointDict];

	For[i=1, i<=Length[dimensions], i++,
		dim = dimensions[[i]];
		\[Rho]Dependence = fixedPointDict[dim];
		\[Rho]Points = Keys[\[Rho]Dependence];
		graphPoints = Table[{\[Rho]Points[[j]], \[Eta] /. \[Rho]Dependence[\[Rho]Points[[j]]]}, {j, 1, Length[\[Rho]Points]}];
		graph = ListPlot[graphPoints, AxesLabel->{\[Rho], \[Eta]}, GridLines->Automatic, PlotRange->{0.,0.3}];
		Export[graphDirectory <> ToString[NumberForm[dim,{5,4}]]<>".png", graph];
		];
];



PlotDimensionDependence[exponentsDict_, cacheDirectory_]:= Block[
	{graphDirectory, exponents, exp, graphDict, graph, i},
	graphDirectory = cacheDirectory <> "exponents/";
	
	If[!DirectoryQ[graphDirectory], CreateDirectory[graphDirectory]];
	
	exponents = Keys[exponentsDict];
	For[i=1, i<=Length[exponents], i++,
		exp = exponents[[i]];
		graphDict = exponentsDict[exp];
		graph = ListPlot[graphDict, AxesLabel->{"d", ToString[exp]}, GridLines->Automatic];
		Export[graphDirectory <> ToString[exp]<>".png", graph];
		];
	];


PlotRegulatorDependence[regulatorDependence_, cacheDirectory_] :=
	Block[{graphDirectory, dimensions, graph, data, exponent, exponentsLabels, dim, i, j},
	graphDirectory = cacheDirectory <> "regulator/";
	
	If[!DirectoryQ[graphDirectory], CreateDirectory[graphDirectory]];
	
	dimensions = Keys[regulatorDependence];
	
	For[j=1, j<=Length[dimensions], j++,
		dim = dimensions[[j]];
		exponentsLabels = Keys[regulatorDependence[dim]];
		For[i=1, i<=Length[exponentsLabels], i++,
			exponent = exponentsLabels[[i]];
			data = regulatorDependence[dim][exponent];
			graph = ListPlot[data, AxesLabel->{"\[Alpha]", ToString[exponent]}, GridLines->Automatic];
			Export[graphDirectory <> ToString[exponent] <> ToString[dim] <> ".png", graph];
		];
	];
];


PlotDataSeries[data_, axesLabels_, legendLabel_, cacheDirectory_]:= Block[
	{graphDirectory, exponents, exp, graphDict, graph, i},
	graphDirectory = cacheDirectory <> "exponents/";
	
	If[!DirectoryQ[graphDirectory], CreateDirectory[graphDirectory]];
	
	graph = ListPlot[Values[data], AxesLabel->axesLabels, GridLines->Automatic, LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 16],
		PlotLegends->PointLegend[Keys[data], LegendLabel->legendLabel, LegendFunction->"Panel", LabelStyle -> Directive[FontFamily -> "Arial", FontSize -> 20]],
		PlotRange->All, ImageSize->Large];
	Export[graphDirectory <> ToString[axesLabels[[2]]]<>".png", graph];
	];
