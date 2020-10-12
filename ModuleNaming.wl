(* ::Package:: *)

EQUATIONSDIRECTORY = "flow_equations/";


GetEquationsFileName[task_, fieldOrder_, anisotropyOrder_] :=
	If[MemberQ[ONTASKS,task],
		StringJoin[{EQUATIONSDIRECTORY,
			"ON_", If[fieldOrder==0, "functional", {"field", ToString[fieldOrder]}], ".wdx"}],
		StringJoin[{EQUATIONSDIRECTORY,
			If[anisotropyOrder!=0, {task, "_"}, {}],
			If[fieldOrder==0, "functional", {"field", ToString[fieldOrder]}],
			If[anisotropyOrder!=0,{"_anisotropy", ToString[anisotropyOrder]},"_isotropic"],
			".wdx"}]];


GetCacheDirectory[task_, fieldOrder_, regulatorName_] :=
Block[{directoryName,dirPrefix=task},
If[task=="regulator", dirPrefix="ON"];

directoryName = StringJoin[{dirPrefix,"_", 
	If[fieldOrder==0, "functional", {"field", ToString[fieldOrder]}],
	"_", regulatorName, "/"}];
If[! DirectoryQ[directoryName], CreateDirectory[directoryName]];
Return[directoryName]; 
];


INITIALGUESSFILE = "initial_guess.wdx";
LOWDGUESSFILE = "lowd_guess.wdx";
TRICRITICALGUESSFILE = "tricrit_guess.wdx";
FIXEDPOINTSFILE = "fixed_points.wdx";
FIXEDPOINTSFILES[i_] := "fixed_points"<>ToString[i]<>".wdx";
DIMENSIONDEPENDENCEFILE = "dimension_dependence.wdx";
REGULATORDEPENDENCEFILE = "regulator_dependence.wdx";
ONRUNCONF = "on_run_configuration.wdx";
Z4RUNCONF = "z4_run_configuration.wdx";
LOWDRUNCONF = "lowd_run_conf.wdx";
TRICRITICALRUNCONF = "tricrit_run_configuration.wdx";


ONTASKS = {"ON","tricritical","lowd","regulator"};
