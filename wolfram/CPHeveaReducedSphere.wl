BeginPackage["MyFun`"]

CPHeveaReducedSphere::usage =
  "CPHeveaReducedSphere[] 按 Bartzos、Borrelli、Denis、Lazarus、Rohmer 与 Thibert (2017) 的特征流凸积分公式，对小球内的短球面实施三个 corrugation，输出论文图 9 中 f1,3 型有限离散实验。默认采用已认证的 4000×20000 网格、可见褶皱数 {21,142,997}、BallRadius 0.52、Eta 0.5 与 TargetMetricFraction 0.237；每个振荡方向必须至少有 8 个采样点。OutputMode 默认 \"Preview\"，只写四个 sampled_ 全局文件以避免约 30 GiB 的重复 ASCII 阶段 VTK；需要完整阶段 VTK 时显式指定 \"Full\"。ProcessSucceeded 表示 native 进程完成，NumericallyValid 还要求 primitive coordinate、flow Jacobian、包围半径以及默认论文配置的指标带通过。论文的初始 profile 与具体过渡参数未公开；本实现使用作者公开 WRL 包络约束的十参数十三次 Hermite profile，故不是其未公开参数的逐位恢复。严格的 C1 等距嵌入是无限阶段极限 fInfinity。";

Begin["`Private`"]

ClearAll[CPHeveaReducedSphere, heveaReducedSphereLibraryPath,
  heveaReducedSphereFunctionCache, ensureHeveaReducedSphereFunction,
  parseHeveaReducedSphereNumber, parseHeveaReducedSphereValue];

Options[CPHeveaReducedSphere] = {"LongitudeSamples" -> 4000,
  "LatitudeSamples" -> 20000, "VisibleRidgeCounts" -> {21, 142, 997},
  "BallRadius" -> 0.52, "Eta" -> 0.5, "TargetMetricFraction" -> 0.237,
  "OutputMode" -> "Preview", "DiskSafetyMarginBytes" -> 1073741824,
  "OutputDirectory" -> Automatic, "TimeLimit" -> 3600};

CPHeveaReducedSphere::lib = "找不到 reduced-sphere LibraryLink 本地库：`1`。";
CPHeveaReducedSphere::arg = "选项不合法：LongitudeSamples 至少为第二方向可见褶皱数的 8 倍，LatitudeSamples 至少为第三方向可见褶皱数的 8 倍；VisibleRidgeCounts 应为三个正整数，BallRadius 应在 0.1 与 1 之间，Eta 应在 0 与 1.5 之间，TargetMetricFraction 应在 0 与 1 之间，OutputMode 应为 \"Preview\" 或 \"Full\"，DiskSafetyMarginBytes 应为非负整数，TimeLimit 应为 1 到 86400 秒的整数。";
CPHeveaReducedSphere::native = "Reduced-sphere LibraryLink 调用失败。";
CPHeveaReducedSphere::io = "无法创建本次 reduced-sphere 的唯一输出目录：`1`。";

$heveaReducedSpherePackageDirectory = DirectoryName[ExpandFileName[$InputFileName]];
$heveaReducedSphereProjectDirectory = DirectoryName[$heveaReducedSpherePackageDirectory];

heveaReducedSphereLibraryPath[] := Module[{override, candidates},
  override = Quiet@Environment["HEVEA_REDUCED_SPHERE_LIBRARY"];
  candidates = DeleteDuplicates@Join[
    If[StringQ[override] && StringLength[override] > 0, {ExpandFileName[override]}, {}],
    FileNameJoin[{$heveaReducedSphereProjectDirectory, #}] & /@
      {"libheveaReducedSphere_mma_static.so",
       FileNameJoin[{"build", "libheveaReducedSphere_mma_static.so"}],
       FileNameJoin[{"build", "bin", "libheveaReducedSphere_mma_static.so"}]}];
  SelectFirst[candidates, FileExistsQ, First[candidates]]];

heveaReducedSphereFunctionCache = <||>;

ensureHeveaReducedSphereFunction[] := Module[
  {key = "Run", lib = heveaReducedSphereLibraryPath[]},
  If[KeyExistsQ[heveaReducedSphereFunctionCache, key] &&
    Head[heveaReducedSphereFunctionCache[key]] === LibraryFunction,
    Return[heveaReducedSphereFunctionCache[key]]];
  If[!FileExistsQ[lib], Message[CPHeveaReducedSphere::lib, lib]; Return[$Failed]];
  heveaReducedSphereFunctionCache[key] = LibraryFunctionLoad[lib,
    "heveaReducedSphereRun_mma", {"UTF8String", Integer, Integer, Integer,
      Integer, Integer, Real, Real, Real, Integer, Integer, Integer}, "UTF8String"]];

parseHeveaReducedSphereNumber[value_String] := If[
  MemberQ[{"nan", "-nan", "inf", "+inf", "-inf"}, ToLowerCase[value]],
  Indeterminate, Quiet@Check[ToExpression@StringReplace[value,
    RegularExpression["([0-9.])e([+-]?[0-9]+)"] -> "$1*^$2"], Indeterminate]];

parseHeveaReducedSphereValue[value_String] := With[{parts = StringSplit[value, ","]},
  If[Length[parts] == 1, parseHeveaReducedSphereNumber[First[parts]],
    parseHeveaReducedSphereNumber /@ parts]];

CPHeveaReducedSphere[OptionsPattern[]] := Module[
  {nx = OptionValue["LongitudeSamples"], ny = OptionValue["LatitudeSamples"],
   ridges = OptionValue["VisibleRidgeCounts"], ballRadius = N@OptionValue["BallRadius"],
   eta = N@OptionValue["Eta"], fraction = N@OptionValue["TargetMetricFraction"],
   outputMode = OptionValue["OutputMode"],
   diskSafetyMargin = OptionValue["DiskSafetyMarginBytes"],
   directory = OptionValue["OutputDirectory"], timeLimit = OptionValue["TimeLimit"],
   fun, raw, result, metricReports, processSucceeded, numericallyValid,
   baseDirectory, runId, stageReports, paperConfiguration, paperMetricValid,
   parameters, files, manifestFile, manifestTemporary, manifestWritten,
   previewPoints, estimatedOutputBytes, requiredFreeBytes},
  If[!IntegerQ[nx] || nx < 32 || !IntegerQ[ny] || ny < 64 ||
      !MatchQ[ridges, {_Integer?Positive, _Integer?Positive, _Integer?Positive}] ||
      (MatchQ[ridges, {_Integer?Positive, _Integer?Positive, _Integer?Positive}] &&
        (nx < 8 ridges[[2]] || ny < 8 ridges[[3]])) ||
      !NumericQ[ballRadius] || !(0.1 < ballRadius < 1.) ||
      !NumericQ[eta] || !(0. < eta < 1.5) ||
      !NumericQ[fraction] || !Between[fraction, {0., 1.}] || fraction == 0. ||
      !MemberQ[{"Preview", "Full"}, outputMode] ||
      !IntegerQ[diskSafetyMargin] || diskSafetyMargin < 0 ||
      !IntegerQ[timeLimit] || !Between[timeLimit, {1, 86400}],
    Message[CPHeveaReducedSphere::arg]; Return[$Failed]];
  previewPoints = Min[nx, 800] Min[ny, 6000];
  estimatedOutputBytes = Ceiling[4 96 previewPoints +
    If[outputMode === "Full", 4 96 nx ny, 0]];
  requiredFreeBytes = estimatedOutputBytes + diskSafetyMargin;
  runId = CreateUUID["hevea-reduced-sphere-"];
  If[directory === Automatic,
    directory = Quiet@Check[CreateDirectory[FileNameJoin[{$TemporaryDirectory, runId}]], $Failed],
    baseDirectory = ExpandFileName[directory];
    If[!DirectoryQ[baseDirectory], baseDirectory = Quiet@Check[
      CreateDirectory[baseDirectory, CreateIntermediateDirectories -> True], $Failed]];
    directory = If[baseDirectory === $Failed, $Failed, Quiet@Check[
      CreateDirectory[FileNameJoin[{baseDirectory, runId}]], $Failed]]];
  If[directory === $Failed, Message[CPHeveaReducedSphere::io,
    OptionValue["OutputDirectory"]]; Return[$Failed]];
  fun = ensureHeveaReducedSphereFunction[]; If[fun === $Failed, Return[$Failed]];
  raw = Quiet@Check[fun[directory, nx, ny, Sequence @@ ridges, ballRadius, eta,
    fraction, timeLimit, Boole[outputMode === "Preview"], requiredFreeBytes], $Failed];
  If[raw === $Failed, Message[CPHeveaReducedSphere::native]; Return[$Failed]];
  result = Quiet@Check[ImportString[raw, "RawJSON"], $Failed];
  If[!AssociationQ[result], Message[CPHeveaReducedSphere::native]; Return[$Failed]];
  metricReports = Map[Association@Map[With[{pair = StringSplit[#, "="]},
      pair[[1]] -> parseHeveaReducedSphereValue[pair[[2]]]] &,
    StringSplit[StringTrim[#]][[2 ;;]]] &,
    StringCases[result["Log"], RegularExpression["METRIC [^\\n]+"]]];
  processSucceeded = result["ExitCode"] === 0;
  stageReports = If[Length[metricReports] == 4, Rest[metricReports], {}];
  paperConfiguration = nx === 4000 && ny === 20000 && ridges === {21, 142, 997} &&
    Abs[ballRadius - 0.52] < 10^-12 && Abs[eta - 0.5] < 10^-12 &&
    Abs[fraction - 0.237] < 10^-12;
  paperMetricValid = !paperConfiguration || (Length[stageReports] == 3 &&
    And @@ Thread[Abs[Lookup[stageReports, "RoundMean"] - {0.83, 0.73, 0.66}] <= {0.06, 0.06, 0.06}] &&
    And @@ Thread[Abs[Lookup[stageReports, "RoundMax"] - {1.03, 0.95, 0.94}] <= {0.12, 0.12, 0.12}] &&
    And @@ Thread[Abs[Lookup[stageReports, "TargetMean"] - {0.14, 0.07, 0.03}] <= {0.04, 0.03, 0.02}] &&
    And @@ Thread[Abs[Lookup[stageReports, "TargetMax"] - {0.24, 0.16, 0.18}] <= {0.08, 0.06, 0.07}]);
  numericallyValid = Length[stageReports] == 3 &&
    AllTrue[stageReports, NumericQ[Lookup[#, "MinPrimitive", Indeterminate]] &&
      Lookup[#, "MinPrimitive"] > 0 &&
      NumericQ[Lookup[#, "MinFlowJacobian", Indeterminate]] &&
      Lookup[#, "MinFlowJacobian"] > 0 &] &&
    AllTrue[metricReports, NumericQ[Lookup[#, "BoundingRadius", Indeterminate]] &&
      Lookup[#, "BoundingRadius"] <= ballRadius + 10^-6 &] && paperMetricValid;
  parameters = <|"LongitudeSamples" -> nx, "LatitudeSamples" -> ny,
    "VisibleRidgeCounts" -> ridges, "BallRadius" -> ballRadius, "Eta" -> eta,
    "TargetMetricFraction" -> fraction, "OutputMode" -> outputMode,
    "EstimatedOutputBytes" -> estimatedOutputBytes,
    "RequiredFreeBytes" -> requiredFreeBytes, "TimeLimit" -> timeLimit|>;
  files = Sort@Join[FileNames["reduced_sphere_*.vtk", directory],
    FileNames["sampled_reduced_sphere_*.vtk", directory]];
  manifestFile = FileNameJoin[{directory, "run-manifest.json"}];
  manifestTemporary = manifestFile <> ".tmp";
  manifestWritten = Quiet@Check[Export[manifestTemporary, <|"RunID" -> runId,
    "Parameters" -> parameters, "ExitCode" -> result["ExitCode"],
    "ProcessSucceeded" -> processSucceeded, "NumericallyValid" -> numericallyValid,
    "Succeeded" -> (processSucceeded && numericallyValid),
    "MetricReports" -> metricReports, "Files" -> files|>, "RawJSON"];
    RenameFile[manifestTemporary, manifestFile]; True, False];
  If[!manifestWritten && FileExistsQ[manifestTemporary], DeleteFile[manifestTemporary]];
  Join[result, <|"ProcessSucceeded" -> processSucceeded,
    "NumericallyValid" -> numericallyValid,
    "Succeeded" -> (processSucceeded && numericallyValid && manifestWritten),
    "RunID" -> runId, "Parameters" -> parameters, "MetricReports" -> metricReports,
    "ManifestFile" -> If[manifestWritten, manifestFile, Missing["NotWritten"]],
    "Files" -> files|>]];

CPHeveaReducedSphere[___] := (Message[CPHeveaReducedSphere::arg]; $Failed);

End[]
EndPackage[]
