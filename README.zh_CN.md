# DimuonSpectrum2011 MC&Data

本项目是一个使用 CMS Open Data 计算双μ子不变质量谱的简单分析样例。

本项目的分析代码基于 CERN Open Data 门户网站上的 [Example code to produce the di-muon spectrum](http://opendata.web.cern.ch/record/5001) from [a CMS 2011 or 2012 primary dataset](http://opendata.web.cern.ch/record/14)(Geiser, Achim. Dutta, Irene. Hirvonsalo, Harri. Sheeran, Bridget. (2017). CERN Open Data Portal. DOI: 10.7483/OPENDATA.CMS.B8MR.C4A2). 同时，本项目还参考了另一个 GitHub 项目 [DimuonSpectrum2011](https://github.com/cms-opendata-analyses/DimuonSpectrum2011). 并加入了对蒙特卡洛样本的比较分析。分析流程参照 [Searching in CMS Open Data for Dimuon Resonances with Substantial Transverse Momentum](https://arxiv.org/abs/1902.04222) 一文，以下称为例文。

对原始代码的修改如下：

- 为了避免与工作区中任何现有的 `DemoAnalyzer` 插件冲突，将类名从 `DemoAnalyzer` 更改为`DimuonSpectrum2011MC`.
- 文件路径已被修改为相对于配置文件中的相对路径，即它们指向 `datasets` 目录，该目录位于将要运行程序的目录下。
- 增加了例文中的选择条件，具体请参考 `src/DimuonSpectrum2010.cc` 文件



## 运行方式

您可以在 [CMS Open Data VM](http://opendata.web.cern.ch/VM/CMS/2010) 上运行此代码。如果您尚未安装 CMSSW，请执行以下操作：

```
cmsrel CMSSW_5_3_32
```

您也可以在 [CMS Open Data Docker](http://opendata.cern.ch/docs/cms-guide-docker) 的容器上运行此代码。如执行以下操作：

```bash
docker run --name dimu2011mc -it cmsopendata/cmssw_5_3_32 bash
```

如果您已经安装了 CMSSW 或运行了 Docker 容器，请直接运行如下命令：

```bash
cd CMSSW_5_3_32/src
cmsenv
```

接下来，您需要创建一个工作目录，可以将其命名为 `WorkDir` 或任意其他名称。并将本项目下载到此目录下。

```bash
mkdir WorkDir
cd WorkDir
git clone git://github.com/PhyM73/DimuonSpectrum2011MC.git
```

再转到本项目的目录，使用 `scram b` 进行编译。

```bash
cd DimuonSpectrum2011MC
scram b
```

其中输入文件已经在配置文件 ` demoanalyzer_cfg.py` 中定义，并且放置在 `datasets` 目录下。接下来直接运行配置文件即可。

```bash
cmsRun demoanalyzer_cfg.py
```

---

为了方便对全部文件进行分析，您还可以在高通量计算平台上运行本项目。



其输出是一个包含多个直方图的 ROOT 文件，默认情况下命名为 DoubleMu.root, 包含有10000个输入 `event`。您可以使用 ROOT 查看这些内容。您还可以修改 `demoanalyzer_cfg.py` 中的相关部分以选择 `datasets` 中其他的输入文件，并重新运行和比较。或者使用 ROOT 软件将各个输出文件融合起来以获得更全面的分析样本，同时也请注意调节 `demoanalyzer_cfg.py` 中的 `process.maxEvent` 参数。现有的 DoubleMu2010.root 文件可以提供参考。

有关更多详细信息，请参阅 `src/DimuonSpectrum2010.cc` 中的注释。

