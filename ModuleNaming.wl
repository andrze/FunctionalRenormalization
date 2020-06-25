(* ::Package:: *)

EQUATIONSDIRECTORY = "flow_equations/";


GetEquationsFileName[anisotropySymmetry_, fieldOrder_, anisotropyOrder_] :=
	If[anisotropySymmetry == "ON",
	StringJoin[{EQUATIONSDIRECTORY,
	"ON_", If[fieldOrder==0, "functional", {"field", ToString[fieldOrder]}], ".wdx"}],
	StringJoin[{EQUATIONSDIRECTORY,
	If[anisotropyOrder!=0,{"Z", ToString[anisotropySymmetry], "_"},{}],
    If[fieldOrder==0, "functional", {"field", ToString[fieldOrder]}],
    If[anisotropyOrder!=0,{"_anisotropy", ToString[anisotropyOrder]},"_isotropic"],
    ".wdx"}]];


GetCacheDirectory[anisotropySymmetry_, fieldOrder_, regulatorName_] :=
Block[{directoryName},
If[ anisotropySymmetry == "ON",
directoryName = StringJoin[{"ON_",
	If[fieldOrder==0, "functional", {"field", ToString[fieldOrder]}],
	"_", regulatorName, "/"}],
directoryName = StringJoin[{"Z", ToString[anisotropySymmetry],"_",
	If[fieldOrder==0, "functional", {"field", ToString[fieldOrder]}],
	"_", regulatorName, "/"}]];
If[! DirectoryQ[directoryName], CreateDirectory[directoryName]];
Return[directoryName]; 
];


INITIALGUESSFILE = "initial_guess.wdx";
FIXEDPOINTSFILE = "fixed_points.wdx";
FIXEDPOINTSFILES[i_] := "fixed_points"<>ToString[i]<>".wdx";
DIMENSIONDEPENDENCEFILE = "dimension_dependence.wdx";
REGULATORDEPENDENCEFILE = "regulator_dependence.wdx";
ONRUNCONF = "on_run_configuration.wdx";
Z4RUNCONF = "z4_run_configuration.wdx";
