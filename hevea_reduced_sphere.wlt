root = DirectoryName[ExpandFileName[$InputFileName]];
report = TestReport[{
  VerificationTest[
    Get[FileNameJoin[{root, "wolfram", "CPHeveaReducedSphere.wl"}]];
    MemberQ[Names["MyFun`CP*"], "CPHeveaReducedSphere"],
    True,
    TestID -> "public-symbol"
  ],
  VerificationTest[
    Options[MyFun`CPHeveaReducedSphere][[1 ;; 8]],
    {"LongitudeSamples" -> 4000, "LatitudeSamples" -> 20000,
      "VisibleRidgeCounts" -> {21, 142, 997}, "BallRadius" -> 0.52,
      "Eta" -> 0.5, "TargetMetricFraction" -> 0.237,
      "OutputMode" -> "Preview", "DiskSafetyMarginBytes" -> 1073741824},
    TestID -> "paper-inspired-defaults"
  ],
  VerificationTest[
    Module[{directory = CreateDirectory[], stale, result, manifest},
      stale = FileNameJoin[{directory, "sampled_reduced_sphere_stale.vtk"}];
      Export[stale, "stale", "Text"];
      SetEnvironment["HEVEA_CORRECTION_SCALES" -> "0,0,0"];
      result = MyFun`CPHeveaReducedSphere[
        "LongitudeSamples" -> 96, "LatitudeSamples" -> 256,
        "VisibleRidgeCounts" -> {2, 3, 8},
        "OutputMode" -> "Full", "OutputDirectory" -> directory, "TimeLimit" -> 120];
      SetEnvironment["HEVEA_CORRECTION_SCALES" -> None];
      manifest = Import[result["ManifestFile"], "RawJSON"];
      {result["ProcessSucceeded"], Length[result["MetricReports"]],
        result["NumericallyValid"], result["Succeeded"],
        Length[result["Files"]], And @@ (FileExistsQ /@ result["Files"]),
        AllTrue[result["Files"],
          FileNameSplit[DirectoryName[#]] ===
            FileNameSplit[result["OutputDirectory"]] &],
        !MemberQ[result["Files"], stale], manifest["Files"] === result["Files"]}],
    {True, 4, True, True, 8, True, True, True, True},
    TestID -> "native-smoke-unique-manifest"
  ],
  VerificationTest[
    Module[{result = MyFun`CPHeveaReducedSphere[
        "LongitudeSamples" -> 96, "LatitudeSamples" -> 256,
        "VisibleRidgeCounts" -> {2, 3, 8},
        "OutputDirectory" -> CreateDirectory[], "TimeLimit" -> 120]},
      {AssociationQ[result], result["ProcessSucceeded"],
        result["NumericallyValid"], result["Succeeded"],
        FileExistsQ[result["ManifestFile"]]}],
    {True, False, False, False, True},
    TestID -> "native-structured-failure-kernel-survives"
  ],
  VerificationTest[
    Module[{result = MyFun`CPHeveaReducedSphere[
        "LongitudeSamples" -> 96, "LatitudeSamples" -> 256,
        "VisibleRidgeCounts" -> {2, 3, 8},
        "DiskSafetyMarginBytes" -> 10^18,
        "OutputDirectory" -> CreateDirectory[], "TimeLimit" -> 120]},
      {AssociationQ[result], result["ExitCode"], result["ProcessSucceeded"],
        result["Succeeded"], Length[result["Files"]],
        StringContainsQ[result["Log"], "insufficient disk space"],
        FileExistsQ[result["ManifestFile"]]}],
    {True, 75, False, False, 0, True, True},
    TestID -> "disk-preflight-structured-failure"
  ],
  VerificationTest[
    Module[{result = MyFun`CPHeveaReducedSphere[
        "LongitudeSamples" -> 1136, "LatitudeSamples" -> 7976,
        "VisibleRidgeCounts" -> {21, 142, 997},
        "OutputDirectory" -> CreateDirectory[], "TimeLimit" -> 1]},
      {AssociationQ[result], result["ExitCode"] >= 128,
        result["ProcessSucceeded"], result["Succeeded"],
        FileExistsQ[result["ManifestFile"]]}],
    {True, True, False, False, True},
    TestID -> "timeout-structured-failure"
  ],
  VerificationTest[
    Quiet[MyFun`CPHeveaReducedSphere[
      "LongitudeSamples" -> 96, "LatitudeSamples" -> 256,
      "VisibleRidgeCounts" -> {2, 3, 8},
      "OutputDirectory" -> "/proc/hevea-reduced-sphere-test"],
      MyFun`CPHeveaReducedSphere::io],
    $Failed,
    TestID -> "unwritable-output-directory"
  ],
  VerificationTest[
    MyFun`CPHeveaReducedSphere["LongitudeSamples" -> 12],
    $Failed,
    {MyFun`CPHeveaReducedSphere::arg},
    TestID -> "invalid-grid"
  ],
  VerificationTest[
    MyFun`CPHeveaReducedSphere[
      "LongitudeSamples" -> 800, "LatitudeSamples" -> 4000,
      "VisibleRidgeCounts" -> {21, 142, 997}],
    $Failed,
    {MyFun`CPHeveaReducedSphere::arg},
    TestID -> "reject-undersampled-paper-frequencies"
  ],
  VerificationTest[
    Module[{directory = CreateDirectory[], result, file, points, bounds},
      result = MyFun`CPHeveaReducedSphere[
        "LongitudeSamples" -> 96, "LatitudeSamples" -> 256,
        "VisibleRidgeCounts" -> {2, 3, 8},
        "OutputDirectory" -> directory, "TimeLimit" -> 120];
      file = SelectFirst[result["Files"],
        StringContainsQ[#, "sampled_reduced_sphere_stage=0"] &];
      points = FirstCase[Import[file], GraphicsComplex[p_, ___] :> p, {}, Infinity];
      bounds = MinMax /@ Transpose[points];
      Max[Abs@Flatten[bounds]] <= 0.501],
    True,
    TestID -> "caps-remain-inside-small-ball"
  ]
}];
Print["sphere_wlt=", If[report["TestsFailedCount"] == 0, "pass", "fail"],
  " succeeded=", report["TestsSucceededCount"],
  " failed=", report["TestsFailedCount"]];
If[report["TestsFailedCount"] != 0,
  Print[ToString[report["TestsFailed"], InputForm]]; Exit[1]];
