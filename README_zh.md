# Campylotropis 物种复合体：分析代码与支持数据

本项目提供论文 **Genomic, ecological, and morphological insights into speciation with introgression in the *Campylotropis macrocarpa* complex** 的分析脚本、处理后数据及主要结果。

[English](README.md) · [分析指南](docs/ANALYSIS_GUIDE.md) · [方法说明](docs/METHODS.md) · [数据与代码获取](docs/DATA_AND_CODE_AVAILABILITY.md)

## 项目内容

分析涵盖群体结构、系统发育、种群历史、基因组分化与等位基因共享、形态关联、基因型–环境关联及空间遗传关系。各分析可独立运行；从原始测序数据开始的处理需要另行获取上游数据及软件。

- `scripts/`：分析及绘图脚本。
- `data/`：处理后的输入、注释、样本信息及模型文件。
- `results/`：统计汇总与主要分析结果。
- `docs/`：方法、依赖软件及运行说明。
- `vendor/`：附带署名与许可证的第三方代码。
- `CONTENTS.tsv`：文件说明、大小及 SHA-256 校验值。

## 分析设置

- 形态 RDA 与基因组 pRDA 使用同一套气候数据，预测变量为 bio3、bio7、bio10、bio13、bio15 和 srad_07。
- 形态 RDA 包含 28 个居群的 192 个完整个体。R² = 34.76%，adjusted R² = 32.64%；前两轴分别解释总变异的 20.98% 和 10.45%。
- 基因组 pRDA 控制遗传 PC1–PC3；候选 SNP 的载荷超过 3.5 个标准差，且基因型–环境相关系数绝对值大于 0.5。
- GEMMA 使用亲缘关系矩阵和遗传 PC1–PC3；BayPass 对过滤后的全基因组 SNP 进行环境关联检验。
- Dsuite Dinvestigate 使用 50 个 usable SNP 的窗口和 25 个 SNP 的步长。物种层面采用 MY1–CYU–CPO 与 CY32–CPO–CYU 两套配置，并在 10-kb 网格中汇总，同时提供参考配置敏感性分析。
- ADMIXTURE 使用 MAF 0.01、位点检出率 0.90、最小间距 1 kb 的基因间 SNP；对应图为 `maf01ms09_re.pdf`。

## 运行示例

在项目根目录执行：

```bash
Rscript scripts/morphology_rda_pca.R
python scripts/gea_overlap.py
python scripts/summarize_fsc.py
```

以上分析使用项目内数据。其他脚本的输入与软件要求见[分析指南](docs/ANALYSIS_GUIDE.md)。Fig. 3、Fig. S7、Table S9 和 Fig. 5d 对应的脚本见[基因流分析说明](docs/INTROGRESSION.md)。

## 数据获取

参考基因组及相关测序数据存于 CNGBdb，项目编号为 **CNP0005557**；全基因组重测序数据的项目编号为 **CNP0007152**。

代码与处理后数据仓库：[Fengyaa/campylotropis-species-complex](https://github.com/Fengyaa/campylotropis-species-complex)。第三方软件和数据库的署名与许可证说明见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
