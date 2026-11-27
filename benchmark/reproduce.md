[Back to the benchmarks](index.md)

# How to reproduce

- **Data:** the benchmark's datasets, fold files and rejection files, Zenodo record 3357167.
- **Folds:** the A tasks as described under [Folds](methods.md#folds); the B tasks and the rejection benchmark with the
  benchmark's own files.
- **SCARLET:** SCARLET 0.6.3 with `blocklist none`, every other option at its default: `learn` + `assign` (A;
  abstention Parts 2–3) and `transfer` (B; abstention Part 1); quick mode with `--draws 1 --seeds S`, S = 1 to 10.
  The same inputs, options and seeds give the same calls, whatever `--ncpu`.
- **Benchmarking script and run records:** the other tools' wrappers and the per-run records are not part of this
  repository; they are available on request.

## The inputs and their SHA-256

Every file the formal benchmark runs of SCARLET 0.6.3 read, with its SHA-256. These were recorded after the fact: the runs did not hash their inputs when they read them. That they are the files the runs read is shown by the formal 0.6.3 rerun reproducing the earlier tables byte for byte, and by their modification times, before every formal run; from now on every run records its inputs' hashes. The data files are h5ad conversions of the benchmark's CSV files (Zenodo record 3357167); the split files are the benchmark's fold files, or, for the within-dataset tasks, the folds regenerated as described under [Folds](methods.md#folds).

| file | role | used by | bytes | SHA-256 |
| :--- | :--- | :--- | :--- | :--- |
| `within-dataset folds, regenerated: 10x_5cl/folds.json` | split | A/10x_5cl | 131,497 | `b37f853e09cb1547791ee3db991dbbc666ad24b96570ebd2b6231657f595b4e2` |
| `within-dataset folds, regenerated: amb/folds.json` | split | A/amb | 470,502 | `bd5456de5b36b181329568c6db52862e5c68e36b2f1bdffd0333875fa1d7e215` |
| `within-dataset folds, regenerated: amb16/folds.json` | split | A/amb16 | 469,845 | `e73def701d8bed1f70368113ae6daf3f2bb4ba58d111ec8dfa46cfce195ac7a1` |
| `within-dataset folds, regenerated: celseq2_5cl/folds.json` | split | A/celseq2_5cl | 17,254 | `0d21adb7a4e9c2f910e5c933403a21867235116e1a69b51cda86b2c34ece7251` |
| `within-dataset folds, regenerated: pancreas_baron_human/folds.json` | split | A/pancreas_baron_human | 302,870 | `019b60b1583bd1d96495034a4ffe81932a9add8259302023c8544d516a5d4bd6` |
| `within-dataset folds, regenerated: pancreas_baron_mouse/folds.json` | split | A/pancreas_baron_mouse | 61,586 | `ee2029cd5e58d3cf188e4ecf0b08ff957aae767532a49926baa5e76efae20912` |
| `within-dataset folds, regenerated: pancreas_muraro/folds.json` | split | A/pancreas_muraro | 70,894 | `3acd2ec369a5c3e614daa6b4c68bc23a356d83da7771de6822fc11b6fb187fd0` |
| `within-dataset folds, regenerated: pancreas_segerstolpe/folds.json` | split | A/pancreas_segerstolpe | 70,739 | `db2b22370abcda6d3519730c289c35772ba71dfd1644cfc40e5df3d425ef1163` |
| `within-dataset folds, regenerated: pancreas_xin/folds.json` | split | A/pancreas_xin | 46,753 | `d9dc9286241f3a15a9c1eabc14c58cc4629ee44c682d52be8e46086bec328e82` |
| `within-dataset folds, regenerated: tm/folds.json` | split | A/tm | 2,194,059 | `cbb765638b5f084b815bb958260d76241e54ff1585709750d804fa4b036973aa` |
| `within-dataset folds, regenerated: zheng_68k/folds.json` | split | A/zheng_68k | 2,648,257 | `40a471714c6de4a26a8a224ff988c7099caed491eb74113f0b0aaa2abbaafab9` |
| `within-dataset folds, regenerated: zheng_sorted/folds.json` | split | A/zheng_sorted | 764,594 | `b9e374ab6855120770ee500ab9d34e71be83ed54318cb746fa16c87034448579` |
| `Inter-dataset/Brain/MouseV1_MouseALM/MouseV1_MouseALM.h5ad` | data | B/Brain_MouseV1_MouseALM_Labels3, B/Brain_MouseV1_MouseALM_Labels34 | 2,654,930,952 | `4546f548b91e7ccea3f4bca6148ecf2174ea8644f62e2f10b7dd329a6d368c68` |
| `Inter-dataset/CellBench/Combined_10x_CelSeq2_5cl_data.h5ad` | data | B/CellBench | 217,378,114 | `9c69afe6762d4a1926dfa3742f87b2108a9e856fdf162f1500c56169af704027` |
| `Inter-dataset/PbmcBench/10Xv2/10Xv2_pbmc1.h5ad` | data | B/PbmcBench_10Xv2 | 248,010,212 | `b0b491574244f0f283041e9162859e5e4422a691bd641a18a2bd709fda76ca4a` |
| `Inter-dataset/PbmcBench/10Xv3/10Xv3_pbmc1.h5ad` | data | B/PbmcBench_10Xv3 | 200,866,868 | `30aa096e6ccab69fbb121f734b593a080f9f1dbcfc58b79bbc32ba1e0cd8adc3` |
| `Inter-dataset/PbmcBench/CEL-Seq/CL_pbmc1.h5ad` | data | B/PbmcBench_CEL-Seq | 203,747,668 | `f6b30fd3850eb5b658651832580dda3aadca159c2e3b34948b2ed62d43e65720` |
| `Inter-dataset/PbmcBench/Drop-Seq/DR_pbmc1.h5ad` | data | B/PbmcBench_Drop-Seq | 241,193,604 | `63aa70e4f27dd5e21862ca862307949a89530ade9515044847531dfcf8d47da6` |
| `Inter-dataset/PbmcBench/Seq-Well/SW_pbmc1.h5ad` | data | B/PbmcBench_Seq-Well | 184,966,076 | `f460ef60f7de6fdebfba216a7329cd95d94b45c0b64f3abd4c418001e27bc8ed` |
| `Inter-dataset/PbmcBench/Smart-Seq2/SM2_pbmc1.h5ad` | data | B/PbmcBench_Smart-Seq2 | 192,045,404 | `bbd60a21140cd2f101e37f0f9b6d4d7c318b264b4c27e19353b760a3ed734c5c` |
| `Inter-dataset/PbmcBench/inDrop/iD_pbmc1.h5ad` | data | B/PbmcBench_inDrop | 207,679,872 | `74a20a9df304b47af7505ee42eaca3081436c39121989a4db5a9da7911c998ae` |
| `Intra-dataset/AMB/Filtered_mouse_allen_brain_data.h5ad` | data | A/amb | 1,606,369,568 | `6f8bd7f3aaaf617e6cedfdb0165b0450a7774ec8f8320639d271231afa8b7065` |
| `Intra-dataset/CellBench/10x_5cl/10x_5cl_data.h5ad` | data | A/10x_5cl | 199,294,235 | `c85b69502a1c674b2c656e4d8370d4531c58a77908798bde9e4f65d8a1d16b14` |
| `Intra-dataset/CellBench/CelSeq2_5cl/CelSeq2_5cl_data.h5ad` | data | A/celseq2_5cl | 37,311,292 | `960155c731f6aade232617c1e70607c35d7d0fd054af16ce021aac3b444dfca5` |
| `Intra-dataset/Pancreatic_data/Baron Human/Filtered_Baron_HumanPancreas_data.h5ad` | data | A/pancreas_baron_human | 195,508,649 | `161d1396bc591e052f2603fb54355170a75052de879de987b25a3658f8fce81f` |
| `Intra-dataset/Pancreatic_data/Baron Mouse/Filtered_MousePancreas_data.h5ad` | data | A/pancreas_baron_mouse | 37,925,338 | `b5a7fd13647a12909ca11a8388c4fa6f7dc6dcc4ddb47b52f0717e08624ce8cd` |
| `Intra-dataset/Pancreatic_data/Muraro/Filtered_Muraro_HumanPancreas_data.h5ad` | data | A/pancreas_muraro | 132,419,046 | `3363c773a2b41ac37294ec39787ee17db29364491ad9019c41f48c4c82cacfc6` |
| `Intra-dataset/Pancreatic_data/Segerstolpe/Filtered_Segerstolpe_HumanPancreas_data.h5ad` | data | A/pancreas_segerstolpe | 164,311,629 | `53477fb717cbc5d139f8216f21440c35083509ba2dfc176b37a654ecbd9a9d57` |
| `Intra-dataset/Pancreatic_data/Xin/Filtered_Xin_HumanPancreas_data.h5ad` | data | A/pancreas_xin | 100,972,344 | `b7c37b34430567a578951a04cce3a6de08bbc25a723a199b22ed634faa7bfe73` |
| `Intra-dataset/TM/Filtered_TM_data.h5ad` | data | A/tm | 1,266,329,121 | `2cdb1412cbe8e771baae7aabe6c45c095de845388cb41ac0e4014ade61004873` |
| `Intra-dataset/Zheng 68K/Filtered_68K_PBMC_data.h5ad` | data | A/zheng_68k | 444,029,751 | `11a3e5aa222f51d9f82ef966ed1145b2cc6d96b41b8d6a5faf86f8d0d33b2b5e` |
| `Intra-dataset/Zheng sorted/Filtered_DownSampled_SortedPBMC_data.h5ad` | data | A/zheng_sorted | 154,495,424 | `5e42d67c8c1950be909fb43d5374301458e42093d87c8641e697819bdbb6190f` |
| `Rejection/Negative_controls/Human/Human_PBMC_Pancreas.h5ad` | data | REJ__part1_a, REJ__part1_b | 233,416,810 | `a9987e67fd2a8e217f61a10a8dacb847ae7468b47d4f727b2a7c40cb3d679f4f` |
| `Rejection/Negative_controls/Mouse/Mouse_Cortex_Pancreas.h5ad` | data | REJ__part1_c, REJ__part1_d | 1,267,165,172 | `f8fcbf00c8d6eeca036efe0ac73fe2455b28619f1edbe0915f157571e669e40e` |
| `Rejection/Unseen populations/Rejection_data.h5ad` | data | REJ_part2, REJ_part3 | 50,589,830 | `5ef542570db99c74a041be8a92cd2bb317c9fc5207c6b2183d6b129a02ad6404` |
| `agy_inter/Brain/MouseV1_MouseALM_Labels3/cv_folds.json` | split | B/Brain_MouseV1_MouseALM_Labels3 | 722,453 | `5752688c3709547bc4236e986d676d3400e081b13b45000932a61b04ef506f83` |
| `agy_inter/Brain/MouseV1_MouseALM_Labels34/cv_folds.json` | split | B/Brain_MouseV1_MouseALM_Labels34 | 722,454 | `544bba709aac87ddacc3e140eae6e2a017e1bce3e2ff2f4985fd1faf4a04c1ad` |
| `agy_inter/CellBench/cv_folds.json` | split | B/CellBench | 76,623 | `a36938ead14af4f4d5adaf5bfda08372557a714c73abaf0e19f550b690145996` |
| `agy_inter/PbmcBench/10Xv2/cv_folds.json` | split | B/PbmcBench_10Xv2 | 521,590 | `c21cd37ad4abe6970a1e91dfca110b98775a914cd443028e2c2551afb2bf972a` |
| `agy_inter/PbmcBench/10Xv3/cv_folds.json` | split | B/PbmcBench_10Xv3 | 441,918 | `1c27d93eb7d4705fd5aa436c58af63ad9777e59217ff983df1cf764721603e2c` |
| `agy_inter/PbmcBench/CEL-Seq/cv_folds.json` | split | B/PbmcBench_CEL-Seq | 443,390 | `66664ffe15f9b0c829f6f3c2f8508c0101af8026535cf6b60aaddf7e6a82cdb8` |
| `agy_inter/PbmcBench/Drop-Seq/cv_folds.json` | split | B/PbmcBench_Drop-Seq | 521,590 | `1e37728c2e9b272d4f55a763a9d37e72c29eb84244affedc85c4beb212fd4911` |
| `agy_inter/PbmcBench/Seq-Well/cv_folds.json` | split | B/PbmcBench_Seq-Well | 425,266 | `b90d9f74da8359811bfbc1a1a53e2f4e98c3ec5496167bb4ec43237de5358996` |
| `agy_inter/PbmcBench/Smart-Seq2/cv_folds.json` | split | B/PbmcBench_Smart-Seq2 | 423,426 | `30892425f664d170d512772f11ee926f5841096786c2dbbc18fe1d3691e1a93a` |
| `agy_inter/PbmcBench/inDrop/cv_folds.json` | split | B/PbmcBench_inDrop | 491,184 | `5a2acd3b8b055055d8d1179bee69574fcaf7da0e07ef5a459962dc3d49d8378f` |
| `rejection part2: Rejection_folds.json` | split | REJ_part2 | 172,758 | `f9aec58c84e642a0e96bbc55ea5835236b5ffb0972941e5c50f677dab83079dd` |
| `rejection part3: split.json` | split | REJ_part3 | 172,771 | `83d28cc2b103a312a941a5a9cafa688f9acdec9c370539730c477189a1f6f850` |
| `derived: Pancreatic_unaligned/Pancreatic_unaligned.h5ad` | data | NEW/Pancreatic_unaligned | 253,094,566 | `0578993642f48b146368b5c215717c37aa440b5eedffd1e17144edbcf76fa96b` |
| `derived: Pancreatic_unaligned/folds.json` | split | NEW/Pancreatic_unaligned | 300,796 | `c309cb8ae209103b3aaf6e4236e4d2606abf09140d6f5174f0a40e96477f5b8b` |
| `derived: inputs/amb16.h5ad` | data | A/amb16 | 1,604,186,913 | `731c74a9ebefef89bcfd71c3c2300192c9f7e22c3550f1f01732352bdf06483e` |
