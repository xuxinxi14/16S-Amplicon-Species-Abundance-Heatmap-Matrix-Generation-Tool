# 16S 扩增子分类单元（Taxa）丰度热图矩阵生成工具

[English](README.md) | [中文](README_CN.md)

本工具用于将 16S 扩增子分析得到的 **OTU/ASV 丰度表**、**分类注释表（taxonomy）** 与 **样本分组信息（metadata）** 整理为可直接用于绘制热图的 **CSV 矩阵**。

支持功能：

- 按指定分类水平汇总（`phylum / class / order / family / genus`）
- 可按样本输出或按分组输出（组内相对丰度均值）
- 基于**平均相对丰度**的低丰度过滤
- 热图常用标准化：`none / zscore / minmax`
- 两种输出方向：taxa×样本(组) 或 样本(组)×taxa
- 支持拼接上级分类注释（门/纲/目/科等链式信息）
- `taxonomy.txt` 同时支持 **7列标准格式** 与 **单列 lineage（分号分隔）格式** 自动解析

---

## 快速开始（推荐）

1）准备 3 个输入文件（与 `generate_heatmap_matrix.py` 同目录，或在 `CONFIG` 中填入路径）：

- `metadata.txt`
- `otutab.txt`
- `taxonomy.txt`

2）安装依赖：

```bash
pip install pandas numpy scipy
```

3）运行脚本：

```bash
python generate_heatmap_matrix.py
```

4）输出结果：

- 默认输出为 `heatmap_matrix.csv`（可通过 `CONFIG["output_path"]` 修改）

---

## 环境要求

- Python **3.8+**
- 依赖：`pandas`、`numpy`、`scipy`

也可用 conda：

```bash
conda install pandas numpy scipy
```

---

## 参数配置（CONFIG）

所有配置项都在 `generate_heatmap_matrix.py` 顶部的 `CONFIG` 字典中。

### 常用配置项

1）**选择分类水平（OTU/ASV 汇总到哪个层级）**

```python
"tax_level": "genus"
```

可选值：`phylum` | `class` | `order` | `family` | `genus`

2）**按样本输出 or 按分组输出**

```python
"aggregate_by": "group"
```

可选值：
- `sample`：按样本输出
- `group`：按分组输出（组内平均相对丰度）

3）**输出方向**

```python
"orientation": "vertical"
```

可选值：
- `vertical`（默认/推荐）：行=taxa；列=分类信息列 + 样本/分组列  
- `horizontal`：行=样本/分组；列=taxa

4）**标准化方式**

```python
"normalization": "zscore"
```

可选值：`none` | `zscore` | `minmax`

5）**低丰度过滤阈值**

```python
"abundance_threshold": 0.001
```

说明：对每个 taxa 计算其在所有列（样本/分组）中的**平均相对丰度**，低于阈值则过滤掉。  
- 设为 `0` 或 `None` 可关闭过滤。

6）**上级分类注释拼接（链式注释）**

```python
"higher_taxonomy_levels": ["phylum", "class", "order", "family"]
```

说明：
- 只会保留**高于** `tax_level` 的层级；无效/低于或等于当前层级的条目会被自动忽略
- 设为 `[]` 表示不拼接上级分类（只输出当前层级名称）

7）**输入/输出文件路径**

```python
"metadata_path": "metadata.txt",
"otutab_path": "otutab.txt",
"taxonomy_path": "taxonomy.txt",
"output_path": "heatmap_matrix.csv",
```

---

## 输入文件格式

所有输入文件均为 **Tab 分隔（TSV）**。

### 1）`metadata.txt`（样本信息表）

必须包含：
- `sampleID`
- `group`

列名大小写不敏感（脚本会自动识别并规范化）。允许存在额外表型列。

示例：

```
sampleID	group	age	sex
Sample1	Control	25	M
Sample2	Control	30	F
Sample3	Treatment	28	M
```

### 2）`otutab.txt`（OTU/ASV 丰度表）

- 第一列：OTU/ASV ID（列名可为空或 `#OTU ID`）
- 其余列：样本名（必须与 metadata 的 sampleID 对上）

示例：

```
#OTU ID	Sample1	Sample2	Sample3
OTU001	100	200	50
OTU002	0	300	120
```

### 3）`taxonomy.txt`（分类注释表）

支持两种格式：

**格式 A（推荐）：7 列标准格式**
- 第一列：OTU/ASV ID
- 分类列：`kingdom / phylum / class / order / family / genus / species`

示例：

```
#OTU ID	kingdom	phylum	class	order	family	genus	species
OTU001	k__Bacteria	p__Firmicutes	c__Bacilli	o__...	f__...	g__...	s__...
```

**格式 B：单列 lineage（分号分隔）**
如果 taxonomy 只有一列分类信息，脚本会自动解析类似：

```
k__Bacteria;p__Firmicutes;c__Bacilli;o__...;f__...;g__...;s__...
```

并拆成标准 7 列。

---

## 匹配规则（非常重要）

脚本会对 ID 做交集匹配，并把不匹配项剔除（会打印警告）：

- `metadata.sampleID` 与 `otutab` 的样本列名必须有���集，否则程序退出
- `otutab` 的 OTU/ASV ID 与 `taxonomy` 的 OTU/ASV ID 必须有交集，否则程序退出

因此：
- 只在其中一个文件出现的样本会被删除
- 没有分类注释的 OTU/ASV 会被删除

---

## 输出说明

输出为带 BOM 的 UTF-8 编码 CSV（`utf-8-sig`），更适合 Excel 直接打开。

数值含义：
- `none`：相对丰度（0~1）
- `zscore`：按 taxa（行）做 Z-score（可为负；方差=0 的行会置 0）
- `minmax`：按 taxa（行）做 Min-Max（0~1；range=0 的行会置 0）

### 输出方向 = `vertical`（默认）

- 行：taxa
- 列：先输出分类列（从拼接标签中拆分出来并去掉 `p__/g__` 等前缀），再输出样本/分组列

### 输出方向 = `horizontal`

- 行：样本/分组
- 列：taxa

---

## 兼容的热图绘图工具

- TBtools
- R `pheatmap`
- GraphPad Prism
- Origin
- Python `seaborn` / `matplotlib`

---

## 常见问题（FAQ）

**Q1：提示“文件不存在”**
- 检查 `CONFIG` 路径，建议用绝对路径
- Windows 下请用 `\` 或 `/`，不要直接写单个反斜杠

**Q2：提示“metadata 与 otutab 无任何共有样本”**
- 确保 `metadata.sampleID` 与 `otutab` 列名完全一致（区分大小写、空格、特殊字符）

**Q3：提示“过滤后无任何物种保留”**
- 降低 `abundance_threshold`（如 `0.0001`）或设为 `0` 关闭过滤

**Q4：Excel 打开中文乱码**
- 输出为 `utf-8-sig`，一般可直接识别；仍乱码则用“数据→自文本/CSV 导入”并选择 UTF-8

**Q5：热图里出现很多 NaN**
- Z-score 时若某 taxa 行方差为 0，脚本会自动将该行置 0，不影响绘图

---

## 方法说明 / 合规性

- 相对丰度：按样本总 reads 归一化（与常见流程一致）
- Z-score：按 taxa（行）维度标准化，使用样本标准差（`ddof=1`）
- 低丰度过滤：基于 taxa 平均相对丰度