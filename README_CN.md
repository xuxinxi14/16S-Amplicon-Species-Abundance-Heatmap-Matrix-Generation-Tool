# 16S扩增子物种丰度热图矩阵生成工具 使用说明

***
下载generate_heatmap_matrix.py或者拷贝代码到本地

## 一、环境配置

需使用 Python 3.8 及以上版本，执行以下命令安装依赖包：

```bash
pip install pandas numpy scipy
```

也可使用 conda 安装：

```bash
conda install pandas numpy scipy
```

***

## 二、参数配置方法

所有配置项均在代码顶部的 `CONFIG` 字典中修改，常用配置修改示例如下：

1.  **切换分类水平（示例改为门水平）**
    ```python
    "tax_level": "phylum"
    ```
    可选值：`phylum`（门） | `class`（纲） | `order`（目） | `family`（科） | `genus`（属）

2.  **按分组输出（计算组内平均相对丰度）**
    ```python
    "aggregate_by": "group"
    ```
    可选值：`sample`（按样本输出） | `group`（按分组输出）

3.  **关闭标准化（输出原始相对丰度）**
    ```python
    "normalization": "none"
    ```
    可选值：`none`（无标准化） | `zscore`（Z-score标准化） | `minmax`（最小-最大归一化）

4.  **关闭低丰度过滤**
    ```python
    "abundance_threshold": 0
    ```

5.  **修改输入文件路径**
    ```python
    "metadata_path": "path/to/your/metadata.txt"
    ```

6.  **关闭上级分类注释拼接**
    ```python
    "add_higher_taxonomy": False
    ```

***

## 三、输入文件格式要求

所有输入文件均为**制表符（Tab）分隔**的文本文件（TSV格式），具体要求如下：

### 3.1 metadata.txt 样本信息表

必须包含 `sampleID`（样本ID）和 `group`（分组）两列，其余表型列可自由扩展，示例如下：

    sampleID    group       age     sex
    Sample1     Control     25      M
    Sample2     Control     30      F
    Sample3     Treatment   28      M

### 3.2 otutab.txt OTU/ASV丰度表

第一列为OTU/ASV ID（无列名、或列名为`#OTU ID`均可），其余列均为样本名，示例如下：

    #OTU ID    Sample1    Sample2    Sample3
    OTU001     100        200        50
    OTU002     0          300        120

### 3.3 taxonomy.txt 物种分类注释表

第一列为OTU/ASV ID，列名为7个标准分类层级（kingdom/phylum/class/order/family/genus/species），示例如下：

    #OTU ID    kingdom      phylum          class           ...
    OTU001     k__Bacteria  p__Firmicutes   c__Bacilli      ...

注意：列名大小写不敏感，程序会自动识别；若列名完全不匹配（如使用中文列名），需手动统一为标准英文分类层级名。

***

## 四、输出文件说明

输出文件为 **UTF-8编码** 的CSV格式文件，第一列为物种名/样本名，其余列为样本名/分组名。

*   标准化方法为`none`时：数值为**相对丰度**（0\~1之间的小数）
*   标准化方法为`zscore`时：数值为Z-score值（可为负数，适配热图配色逻辑）
*   标准化方法为`minmax`时：数值归一化至0\~1区间（另一种标准化方案）

输出文件完美兼容以下主流热图绘制工具：

*   ✅ TBtools（直接导入CSV即可绘制热图）
*   ✅ R pheatmap/ggplot2（`read.csv()`直接读取，转为矩阵后绘图）
*   ✅ GraphPad Prism（导入CSV后直接选择热图模板）
*   ✅ Origin（导入CSV后插入热图）
*   ✅ Python seaborn/matplotlib（`pd.read_csv()`直接读取绘图）

***

## 五、常见问题排查

**Q1：运行报错 "文件不存在"**
A：请检查文件路径是否正确，推荐使用绝对路径；Windows系统路径中的反斜杠`\`需写为双反斜杠`\\`，或直接改用正斜杠`/`。

**Q2：运行报错 "metadata 与 otutab 无任何共有样本"**
A：请检查otutab的样本列名，与metadata的sampleID列内容是否完全一致，**严格区分大小写、空格和特殊字符**。

**Q3：运行报错 "过滤后无任何物种保留"**
A：请降低`abundance_threshold`阈值，例如改为0.0001（即0.01%），或直接设为0关闭低丰度过滤功能。

**Q4：用Excel打开输出文件后中文乱码**
A：输出文件默认使用utf-8-sig编码，Excel可正常识别；若仍出现乱码，请使用Excel「数据→自文本/导入文本」功能，导入时选择UTF-8编码。

**Q5：分类注释表列名不标准（如Phylum而非phylum）**
A：程序会自动识别列名大小写，无需手动修改；仅当列名与标准分类层级名完全不符时（如使用中文列名），才需要手动统一列名。

**Q6：热图中出现大量NaN值**
A：该问题通常是Z-score标准化时，某物种在所有样本中丰度完全一致（方差为0）导致；程序已自动将这类行的数值置0，不影响后续热图绘制。

***

## 六、引用与合规说明

本工具的核心计算逻辑完全符合QIIME2、USEARCH的16S扩增子分析行业标准：

*   相对丰度计算：按样本总reads数归一化，与QIIME2 `feature-table relative-frequency` 结果完全一致
*   Z-score标准化：按物种（行）维度计算，使用样本标准差（ddof=1）
*   低丰度过滤：基于物种在所有样本中的平均相对丰度过滤，与QIIME2 `filter-features` 逻辑一致

