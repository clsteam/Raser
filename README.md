# README
所有的流程文件中包含相应的软件，定义规则如下：
- raser.py决定流程顺序， 享有绝对决定权，其他文件什么的命名顺序只是常规的分析流程顺序。
- 以软件优先原则，继承Tools的子类可以增加相应的功能函数
（不能因为此功能在其他流程中而在此创建新Tools子类）
Installation
******
Note: all tools need to be in the environment variable.

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


