# 16S Amplicon Species Abundance Heatmap Matrix Generation Tool - User Manual

***

## 1. Environment Configuration

Python 3.8+ is required. Install the following dependency packages:

```bash
pip install pandas numpy scipy
```

Or use conda:

```bash
conda install pandas numpy scipy
```

***

## 2. Parameter Configuration Guide

All configuration items can be modified in the `CONFIG` dictionary at the top of the code. Common modification examples are as follows:

1.  **Switch taxonomic level (e.g. change to phylum level)**
    ```python
    "tax_level": "phylum"
    ```
    Available options: `phylum` | `class` | `order` | `family` | `genus`

2.  **Output by group (calculate mean relative abundance within groups)**
    ```python
    "aggregate_by": "group"
    ```
    Available options: `sample` | `group`

3.  **Turn off normalization (output raw relative abundance)**
    ```python
    "normalization": "none"
    ```
    Available options: `none` | `zscore` | `minmax`

4.  **Turn off low-abundance filtering**
    ```python
    "abundance_threshold": 0
    ```

5.  **Change input file path**
    ```python
    "metadata_path": "path/to/your/metadata.txt"
    ```

6.  **Disable higher taxonomy annotation**
    ```python
    "add_higher_taxonomy": False
    ```

***

## 3. Input File Format Requirements

### 3.1 metadata.txt (tab-delimited)

Must contain the `sampleID` column and `group` column. Additional columns can be added freely:

    sampleID    group       age     sex
    Sample1     Control     25      M
    Sample2     Control     30      F
    Sample3     Treatment   28      M

### 3.2 otutab.txt (tab-delimited)

The first column is OTU ID (with no column name or column name `#OTU ID` are both acceptable), and the remaining columns are sample names:

    #OTU ID    Sample1    Sample2    Sample3
    OTU001     100        200        50
    OTU002     0          300        120

### 3.3 taxonomy.txt (tab-delimited)

The first column is OTU ID, and the column names correspond to the 7 taxonomic levels (kingdom/phylum/class/order/family/genus/species):

    #OTU ID    kingdom      phylum          class           ...
    OTU001     k__Bacteria  p__Firmicutes   c__Bacilli      ...

Note: Column names are case-insensitive and will be automatically recognized by the program.

***

## 4. Output File Description

The output is a UTF-8 encoded CSV file, with the first column being species/sample names, and the remaining columns being sample/group names.

*   When normalization method is `none`: Values are relative abundance (decimal between 0 and 1)
*   When normalization method is `zscore`: Values are Z-scores (can be negative, suitable for heatmap color mapping)
*   When normalization method is `minmax`: Values are scaled between 0 and 1 (an alternative normalization method)

Compatible with the following heatmap plotting tools:

*   ✅ TBtools (directly import CSV for heatmap plotting)
*   ✅ R pheatmap (directly read via `read.csv()`, convert to matrix for plotting)
*   ✅ GraphPad Prism (select heatmap after importing CSV)
*   ✅ Origin (import CSV and insert heatmap)
*   ✅ Python seaborn/matplotlib (directly read via `pd.read_csv()`)

***

## 5. Troubleshooting (FAQs)

**Q1: Error "File does not exist"**
A: Check if the file path is correct. Absolute paths are recommended. Note that backslashes in Windows paths need to be written as double backslashes (`\\`) or replaced with forward slashes (`/`).

**Q2: Error "No shared samples between metadata and otutab"**
A: Check if the column names of otutab and the sampleID in metadata are completely consistent (case and space sensitive).

**Q3: Error "No taxa retained after filtering"**
A: Lower the `abundance_threshold` value, e.g. change to 0.0001 (0.01%), or set to 0 to turn off filtering.

**Q4: Chinese garbled characters after importing into Excel**
A: The output file uses utf-8-sig encoding, which can be normally recognized by Excel. If garbled characters still occur, please use Excel's "Data → Import Text" function and select UTF-8 encoding.

**Q5: Non-standard taxonomy column names (e.g. Phylum instead of phylum)**
A: The program automatically recognizes case, but if the column names are completely different (e.g. using Chinese), you need to manually unify them.

**Q6: A large number of NaN values appear in the heatmap**
A: Usually caused by the same abundance of a taxon in all samples (variance = 0) during Z-score normalization. The program has automatically set these rows to 0, which does not affect plotting.

***

## 6. Citation & Compliance

The calculation logic of this tool complies with the industry standards of QIIME2 and USEARCH:

*   **Relative abundance**: Normalized by sample, consistent with QIIME2 `feature-table relative-frequency`
*   **Z-score**: Calculated by taxon (row) dimension, using sample standard deviation (ddof=1)
*   **Low-abundance filtering**: Based on the mean relative abundance of each sample, consistent with QIIME2 `filter-features`

***

