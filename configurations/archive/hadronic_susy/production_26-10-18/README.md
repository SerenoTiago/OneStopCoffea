# Hadronic SUSY production configuration, archived 2026-10-19 (state of 2026-10-18)

Exact copy (`cp -p`, working tree, not git HEAD) of
`configurations/production/hadronic_susy/` immediately before the Run-3 MET
noise filters and jet veto map were added.  It includes the 2026-10-18
dataset update (34 backgrounds: 6 Z->nunu LO, 12 W+jets LO MLM, 10 QCD,
3 ttbar incl. TTto2L2Nu, 3 diboson) and the event selection used by the
2026-10-07 production and the 2026-10-15 freeze.  `__pycache__` not copied.

These files were not edited after copying.  Notes:
- `backgrounds_zjets_rerun.py` loads `Path(__file__).with_name("backgrounds.yaml")`,
  so the archived copy reproduces the archived configuration.  Its usage
  string, `legacy_cuts.yaml` and
  `docs/investigations/dask_exec_silent_data_loss_2026-10-15.md:51` still
  mention the old production path.
- Dataset and era definitions are NOT copied; the configs resolve them from
  `analyzer_resources/` at run time.  Their checksums at archival are below.
- `backgrounds.yaml.diff_vs_HEAD` is `git diff HEAD` of backgrounds.yaml
  (34-sample update + two comment lines removed outside these edits).
- The two comment lines above `location_priorities` that are missing relative
  to git HEAD were removed before archival by an edit not made in this work.

Git HEAD at archival: 876a26d7daf0fb2bafa5635202936263fb6c8ec7

## sha256 (archived files)
```
220d3bc04a010da55748484a5a731d716049dda0f1b7a4282ed2fe84eb6d29c0  backgrounds.yaml
86c254fdf6fd3e296fa482fa7fcc8abcc479d9ec3249f756abccc836ef755446  split_signals.yaml
03e7d9906e71f29c157643c3b4f4a68c1625e78a8fe6c24f9b7dc9cf8829cd3c  postprocessing.yaml
eaca5c1a491d2f74b1b7af6c95e7f975ffb4a73527d32f316ae3d55251822158  backgrounds_zjets_rerun.py
72b5014281f5925e878aaa7da3936ecfdfcad6cfdb0ad1d9c418643bdce2679c  legacy_cuts.yaml
9af16253dac1b1edf068ae5eb8602e7bc11a6771ac8e20feb64237d48bfba7b1  backgrounds.yaml.diff_vs_HEAD
```

## sha256 (resolved dataset and era files at archival)
```
7102b6198b6075fe7214dd12906c32b502da66ed4b1e48b01925f35258a049fa  analyzer_resources/datasets/2024/diboson_2024.yaml
16053fc55510814629f5c8be6d8770ca3ed734f9cdfd58f218725067db0f4970  analyzer_resources/datasets/2024/qcd_inclusive_2024.yaml
1ef4900096c3552aff261d7139f364f3e37b1a2193b80e54e9db10335d6854f7  analyzer_resources/datasets/2024/tt_dileptonic.yaml
bba541889a7fdbc6da4dcba319b56ccbba1338b8c4207e002b528b00e0308b13  analyzer_resources/datasets/2024/tt_hadronic.yaml
7b897e2d88582c52ee350db2f976163ab0bda32662440c58c74100d2dfaaa247  analyzer_resources/datasets/2024/tt_semileptonic.yaml
804e3dda3eb7f5c687c3f5900ed42df52b3062455b55f6dd515092b8a0ef9a56  analyzer_resources/datasets/2024/wjets_2024.yaml
5fd31799b24d93670504f99ab8a3d4214bb14aba3f9ac8ed4028c57f8a091281  analyzer_resources/datasets/2024/zjets_2024.yaml
4ac181e5091ea49239e54002d1ad74ecd8a67a4fff179247300a357381b7c438  analyzer_resources/datasets/2024/2stop_split/2stop_1000_2024_splitLSP.yaml
0ef03ac6419ea19dbbe7f4ee017a71e65c5bdefedb026a4226c8b2e08758f9dc  analyzer_resources/datasets/2024/2stop_split/2stop_1050_2024_splitLSP.yaml
bfad12211a26db5c80f4ab8f6270c0298694a9ac0723cfaa0ac713596a0fd5f1  analyzer_resources/datasets/2024/2stop_split/2stop_1100_2024_splitLSP.yaml
3b952f3131e883dd1e8df8b01bbd7576a91b77a40f21fd2306e27691edf40e4a  analyzer_resources/datasets/2024/2stop_split/2stop_1150_2024_splitLSP.yaml
380f12518ddfa31fb1f64476d980dd2a0c8382622513955814da1766b7922fc5  analyzer_resources/datasets/2024/2stop_split/2stop_1200_2024_splitLSP.yaml
471d34b2edb5dbf93e48e28b3476665530cc3171f7b609a0e00197aabbed5c7c  analyzer_resources/datasets/2024/2stop_split/2stop_1250_2024_splitLSP.yaml
a2661a860728bfb5c8572b0b45d0cd042f64d01681be1780880c6dc01d53dd32  analyzer_resources/datasets/2024/2stop_split/2stop_1300_2024_splitLSP.yaml
3f097a2a2af99115902782b37f3ff5f821f180a52862b63fc91cb4da8d49d42e  analyzer_resources/datasets/2024/2stop_split/2stop_1350_2024_splitLSP.yaml
33d4c74d0ab9266a45135134b60a6036edef51c91b64583e6c728c01c1a312f5  analyzer_resources/datasets/2024/2stop_split/2stop_1400_2024_splitLSP.yaml
1faca9835b91aba4844203c455e52d8eee50dba5d70b9ff419ed3f398dd63990  analyzer_resources/datasets/2024/2stop_split/2stop_1450_2024_splitLSP.yaml
40770c835d1ad56597e9fca18e3803fbfe1c974c309ed55ba4c4f6d8e493ec40  analyzer_resources/datasets/2024/2stop_split/2stop_1500_2024_splitLSP.yaml
50bf82d6c7598464d4668106a2c0f6c40257ba9aa7c1535568f3d1995fb91f31  analyzer_resources/datasets/2024/2stop_split/2stop_1600_2024_splitLSP.yaml
036f9d075d452c3573025aa68a4c5d62321f9905e16ccb39d36fae5b005f4bba  analyzer_resources/datasets/2024/2stop_split/2stop_1700_2024_splitLSP.yaml
942dd41439166691de87f97796225df667b050263c54622bb9918f9c7541cdc1  analyzer_resources/datasets/2024/2stop_split/2stop_1800_2024_splitLSP.yaml
7cff15c0e6911df0b9593b62a23082d1bb1f4fbf03a9ba0d6b72648080d3c38e  analyzer_resources/datasets/2024/2stop_split/2stop_1900_2024_splitLSP.yaml
3b17617d8564e5dc160ea3dd7bf18e9351e4245a0a071e3cc23c8eba2546f69b  analyzer_resources/datasets/2024/2stop_split/2stop_600_2024_splitLSP.yaml
418b6ec2a299c3536a209c77e35f298455e99a4db9f1d0b45545da7d2f2df84d  analyzer_resources/datasets/2024/2stop_split/2stop_700_2024_splitLSP.yaml
3ede9f3d0ea9228aa328eb230f4fda48df42cdb224869729be497c8b280652df  analyzer_resources/datasets/2024/2stop_split/2stop_900_2024_splitLSP.yaml
f49f0af05f3f646414cc9fa096223faf5d77c0fc1786618aaf23a78a7c1f69f4  analyzer_resources/datasets/2024/2stop_split/2stop_950_2024_splitLSP.yaml
54b56c80cf54f3b0beebc068864130382b011d16f424c8944a68081cf97893a0  analyzer_resources/eras/2024.yaml
```
