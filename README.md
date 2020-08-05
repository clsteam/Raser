# Raser
**A pipeline that automatically analyzes RNA-Seq data**
## Introduction
RNA-Seq is a new transcriptome research method, with high efficiency, high sensitivity, and full genome analysis (for any species without pre-designing probes) and other advantages. Currently, a variety of analysis tools have been developed for RNA-Seq data, including data preprocessing, sequence alignment, transcriptome assembly, gene expression estimation, and non-coding RNA detection. However, these analysis tools basically exist independently, lacking a relatively complete system to integrate different tools to complete most of the analysis.

Raser was born from this. He helps you realize most of the software installation-free configuration, parameter configuration, accurate management of multiple samples, complete log management for each sample, and some visualization tasks
## Installing Raser

Raser requires the following software and data resources to be installed. 
>Note, if you can use our [Docker](https://github.com/STAR-Fusion/STAR-Fusion/wiki#Docker)  images, then you'll have all the software pre-installed and can hit the ground running. 

###  1. Downloading from GitHub Clone
```
    $ git clone --recursive git@github.com:clsteam/RASER.git
```
The --recursive parameter is needed to integrate the required submodules.

###  2. Tools Required

 * Raser is developed based on Python 3.7 (if you have not installed it, you can go to the official Python website to download and install), you need to **enter the root directory of Raser** and run the following command to install the required Python dependency packages:

```
    $ pip3 install -r requirements.txt
```

* Raser packages most of the software in a separate folder, and users can use it after downloading, but some software needs to be manually installed and compiled:
    1. 


R
1:FastQC(>=0.11.5)

2:Trimmomatic(>=0.32)

3:Samtools(>=1.3.1)

4:Tophat2(>=2.1.0)

5:Cufflinks(>=2.2.1)

---

#### manager.py
主要提供外部接口参数
- culster：运行环境选在集群中
---

raser/
-
#### settings.py
- 参数配置文件位置
- 工具选择
---

core/
-
#### base.py
提供接口类
#### rdconf.py
提供基本的常量（从配置文件中读取）
- 配置文件句柄conf
- CMD_PBS
- EXPERIMENT_GROUP, CONTROL_GROUP
---

### INSTALL
1. need to export path

    
    RASER_TOOL_HOME = ~/tools
    export PATH=${RASER_TOOL_HOME}/blast

2. need to installed and export path
- samtools

3. python3
- RSeQC
- HTSeq

R
- DESeq2
- Ballgown


