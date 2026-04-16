# 16S Amplicon Taxa Abundance Heatmap Matrix Generation Tool

[English](README.md) | [中文](README.zh-CN.md)

A lightweight Python script that converts 16S amplicon OTU/ASV tables + taxonomy annotation + sample metadata into a **heatmap-ready abundance matrix (CSV)**.

It supports:

- Aggregation to a chosen taxonomic level (**phylum/class/order/family/genus**)
- Optional aggregation by group (mean relative abundance within groups)
- Low-abundance filtering based on **mean relative abundance**
- Normalization for heatmap visualization (**none / zscore / minmax**)
- Two output orientations (taxa × samples/groups, or samples/groups × taxa)
- Taxonomy annotation chains (prepend higher taxonomy levels)
- `taxonomy.txt` in either **7-column format** or **single-column semicolon-separated lineage** format

---

## Quick Start

1. Put these files in the same directory as `generate_heatmap_matrix.py` (or set their paths in `CONFIG`):

- `metadata.txt`
- `otutab.txt`
- `taxonomy.txt`

2. Install dependencies:

```bash
pip install pandas numpy scipy
```

3. Run:

```bash
python generate_heatmap_matrix.py
```

4. Result:

- Output file defaults to `heatmap_matrix.csv` (configurable via `CONFIG["output_path"]`).

---

## Environment

- Python **3.8+**
- Dependencies: `pandas`, `numpy`, `scipy`

Conda alternative:

```bash
conda install pandas numpy scipy
```

---

## Configuration (CONFIG)

All parameters are in the `CONFIG` dictionary at the top of `generate_heatmap_matrix.py`.

### Common options

1) **Switch taxonomic level**

```python
"tax_level": "phylum"
```

Available: `phylum` | `class` | `order` | `family` | `genus`

2) **Aggregate by sample or group**

```python
"aggregate_by": "group"
```

Available: `sample` | `group`

3) **Set output orientation**

```python
"orientation": "vertical"
```

Available:

- `vertical`: rows are taxa; columns are taxonomy columns + samples/groups (default; recommended)
- `horizontal`: rows are samples/groups; columns are taxa

4) **Normalization**

```python
"normalization": "zscore"
```

Available: `none` | `zscore` | `minmax`

5) **Low-abundance filtering**

```python
"abundance_threshold": 0.001
```

Rule: remove taxa whose **mean relative abundance across all columns** is below the threshold.

- Set to `0` or `None` to disable.

6) **Higher taxonomy annotation chain**

```python
"higher_taxonomy_levels": ["phylum", "class", "order", "family"]
```

- Only levels higher than `tax_level` are kept; invalid or lower/equal levels are ignored.
- Use `[]` to disable and output only the current level label.

7) **Input/output paths**

```python
"metadata_path": "metadata.txt",
"otutab_path": "otutab.txt",
"taxonomy_path": "taxonomy.txt",
"output_path": "heatmap_matrix.csv",
```

---

## Input file formats

All inputs are **tab-delimited**.

### 1) `metadata.txt`

Must include columns:

- `sampleID`
- `group`

Column names are matched case-insensitively. Extra columns are allowed.

Example:

```
sampleID	group	age	sex
Sample1	Control	25	M
Sample2	Control	30	F
Sample3	Treatment	28	M
```

### 2) `otutab.txt`

- First column: OTU/ASV ID (header can be empty or `#OTU ID`)
- Remaining columns: sample names

Example:

```
#OTU ID	Sample1	Sample2	Sample3
OTU001	100	200	50
OTU002	0	300	120
```

### 3) `taxonomy.txt`

Supported formats:

**Format A (recommended): 7 columns**

- First column: OTU/ASV ID
- Columns: `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`

Example:

```
#OTU ID	kingdom	phylum	class	order	family	genus	species
OTU001	k__Bacteria	p__Firmicutes	c__Bacilli	o__...	f__...	g__...	s__...
```

**Format B: single-column lineage string**

If `taxonomy.txt` has only one taxonomy column, the script will parse strings like:

```
k__Bacteria;p__Firmicutes;c__Bacilli;o__...;f__...;g__...;s__...
```

into the standard 7 columns automatically.

---

## Matching rules (important)

The script will keep the intersection of identifiers and drop unmatched entries with warnings:

- `metadata.sampleID` ∩ `otutab` sample columns must have overlap; otherwise it exits.
- `otutab` OTU/ASV IDs ∩ `taxonomy` OTU/ASV IDs must have overlap; otherwise it exits.

This means:

- Samples present in only one file are removed.
- OTUs/ASVs without taxonomy annotation are removed.

---

## Output

The output is a **UTF-8 with BOM** CSV (`utf-8-sig`) for better Excel compatibility.

- Default path: `heatmap_matrix.csv`
- Values depend on normalization:
  - `none`: relative abundance in [0, 1]
  - `zscore`: Z-scores (per-taxon row-wise; can be negative)
  - `minmax`: scaled to [0, 1] (per-taxon row-wise)

### Orientation = `vertical` (default)

- Rows: taxa
- Columns:
  - taxonomy columns (split from the concatenated label, prefixes removed)
  - then sample/group columns

### Orientation = `horizontal`

- Rows: samples/groups
- Columns: taxa

### Output example (vertical)

A simplified example (your taxonomy columns depend on `higher_taxonomy_levels` and `tax_level`):

```
phylum,family,genus,Control,Treatment
Firmicutes,Lachnospiraceae,Blautia,0.12,0.34
Bacteroidota,Bacteroidaceae,Bacteroides,0.08,0.02
```

---

## Compatible heatmap tools

- TBtools
- R `pheatmap`
- GraphPad Prism
- Origin
- Python `seaborn` / `matplotlib`

---

## Troubleshooting (FAQ)

**Q1: "File does not exist"**

- Check paths in `CONFIG`. Absolute paths are recommended.
- On Windows, use `\` or `/` instead of `\`.

**Q2: "No shared samples between metadata and otutab"**

- Ensure `metadata.sampleID` exactly matches `otutab` column names (case/whitespace sensitive).

**Q3: "No taxa retained after filtering"**

- Lower `abundance_threshold` (e.g. `0.0001`) or set it to `0` to disable filtering.

**Q4: Chinese garbled characters in Excel**

- The output uses `utf-8-sig`. If Excel still garbles, import via Excel "Data → From Text/CSV" and select UTF-8.

**Q5: Many NaN values in heatmap**

- For Z-score normalization, taxa with zero variance are set to 0 automatically.

---

## Methods / Compliance notes

The calculations follow common practices used in QIIME2 / USEARCH workflows:

- Relative abundance: per-sample normalization, similar to QIIME2 `feature-table relative-frequency`
- Z-score: computed per taxon (row-wise) using sample standard deviation (`ddof=1`)
- Low-abundance filtering: based on mean relative abundance