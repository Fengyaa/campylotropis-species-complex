# Campylotropis 文章脚本与必要数据

本项目已精简为按分析分类的代码与数据档案，不要求整体复现全部文章。已去掉总运行器、持续集成、整套环境配置、重复探索版本及生成的work目录。旧项目在本地另行备份。

## 最终确认

- 形态RDA与基因组pRDA使用同一套修正后的气候数据：bio3、bio7、bio10、**bio13**、bio15、srad_07。最后确认实际使用bio13，未替换为bio12。
- 形态RDA为192个体/28居群，R²=34.76%，adjusted R²=32.64%；两轴占总变异20.98%和10.45%。
- GEMMA未使用环境因子，保留遗传PC和亲缘关系矩阵。
- BayPass最终关联预测使用全基因组过滤后的所有SNP。
- species和sympatric fd均记录为50 SNP窗口、25 SNP步长。
- ADMIXTURE只保留maf01ms09_re.pdf及其实际对应的maf01ms90thin1k Q矩阵；其他MAF版本不打包。

`scripts/`存放分析脚本和整理后的上游命令；`data/`存放必要小型输入、候选表和模型文件；`results/`存放关键汇总和指定ADMIXTURE图。`CONTENTS.tsv`记录逐文件来源和校验值。

新增PSMC、SVDquartets、ADMIXTURE命令来自本次提供的记录，清理了重复调用、文件名笔误及服务器绝对路径。整理后的命令与已实际执行的历史命令有所区别，其他上游命令明确标为参数模板。外部输入要求见[分析说明](docs/ANALYSIS_GUIDE.md)。

PSMC保留所提供的μ=8.17e−8、g=10，未据此重画旧图或更改FSC时间。原vcf2fq.py与原sample_list仍是外部依赖。旧species fd文件的100_50后缀与本次确认的50/25分开记录，保留结果数据值，不根据文件名擅改结果。

项目尚未上传；公开时补充真实仓库地址、文章引用及代码/数据许可证。参数与数据来源见[方法说明](docs/METHODS.md)。
