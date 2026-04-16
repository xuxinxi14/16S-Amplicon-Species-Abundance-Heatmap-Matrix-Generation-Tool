#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
16S扩增子分析 - 物种丰度热图矩阵生成工具
===========================================
功能：读取OTU丰度表、物种注释表和元数据，生成适用于各种热图绘制工具的CSV矩阵文件
作者：徐新玺
版本：1.1.0
日期：2026-04-04

符合QIIME2/USEARCH行业标准，适用于TBtools、R pheatmap、GraphPad Prism、Origin等工具
"""

import os
import sys
import pandas as pd
import numpy as np
from scipy import stats

# ============================================================
# CONFIG 配置块 - 所有可配置项均在此处修改，无需改动核心逻辑
# ============================================================
CONFIG = {
    # ---- 输入文件路径配置 ----
    # metadata文件路径（包含sampleID和group列的制表符分隔txt文件）
    "metadata_path": "metadata.txt",

    # OTU丰度表路径（行为OTU ID，列为样本名的制表符分隔txt文件）
    "otutab_path": "otutab.txt",

    # 物种注释表路径（行为OTU ID，列为7个分类水平的制表符分隔txt文件）
    "taxonomy_path": "taxonomy.txt",

    # ---- 输出文件路径配置 ----
    # 输出CSV文件路径（含文件名）
    "output_path": "heatmap_matrix_new.csv",

    # ---- 分类水平配置 ----
    # 可选值："phylum"(门) | "class"(纲) | "order"(目) | "family"(科) | "genus"(属)
    # 默认：属水平（genus）
    "tax_level": "genus",

    # ---- 聚合维度配置 ----
    # 可选值："sample"(按样本) | "group"(按分组，计算组内相对丰度总和)
    # 默认：按样本（sample）
    "aggregate_by": "group",

    # ---- 数据方向配置 ----
    # 可选值："vertical"(纵向：行=物种，列=样本/分组，热图绘制最常用)
    #        "horizontal"(横向：行=样本/分组，列=物种)
    # 默认：纵向（vertical）
    "orientation": "vertical",

    # ---- 上级分类注释配置 ----
    # 指定需要附加在分类名称前的上级水平列表，顺序即输出顺序
    # 可选水平："phylum"(门) | "class"(纲) | "order"(目) | "family"(科)
    # 规则：只能选比当前 tax_level 更高的水平，低于或等于当前水平的条目会被自动忽略
    # 示例（当前水平为 genus）：
    #   ["phylum", "class", "order", "family"]  → 输出 p__xx;c__xx;o__xx;f__xx;g__xx（完整上级链）
    #   ["phylum", "family"]                    → 只标注门和科，输出 p__xx;f__xx;g__xx
    #   ["phylum"]                              → 只标注门，输出 p__xx;g__xx
    #   []                                      → 不附加任何上级，仅输出当前水平名称
    "higher_taxonomy_levels": ["phylum", "class", "order", "family"],

    # ---- 低丰度过滤配置 ----
    # 过滤所有样本中平均相对丰度低于该阈值的物种
    # 默认：0.001（即0.1%）；设为0或None则关闭过滤
    "abundance_threshold": 0.001,

    # ---- 标准化方法配置 ----
    # 可选值："none"(无标准化) | "zscore"(Z-score标准化，按物种维度) | "minmax"(Min-Max标准化)
    # 默认：Z-score标准化（zscore），热图绘制推荐
    "normalization": "zscore",
}

# ============================================================
# 分类水平映射表（内部使用，勿修改）
# ============================================================
TAX_LEVEL_MAP = {
    "phylum": {"col": "phylum", "prefix": "p__", "idx": 1, "higher": []},
    "class":  {"col": "class",  "prefix": "c__", "idx": 2, "higher": ["phylum"]},
    "order":  {"col": "order",  "prefix": "o__", "idx": 3, "higher": ["phylum", "class"]},
    "family": {"col": "family", "prefix": "f__", "idx": 4, "higher": ["phylum", "class", "order"]},
    "genus":  {"col": "genus",  "prefix": "g__", "idx": 5, "higher": ["phylum", "class", "order", "family"]},
}

TAX_COLUMNS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

TAX_PREFIX_MAP = {
    "kingdom": "k__",
    "phylum":  "p__",
    "class":   "c__",
    "order":   "o__",
    "family":  "f__",
    "genus":   "g__",
    "species": "s__",
}


# ============================================================
# 模块一：日志工具
# ============================================================
def log(msg: str, level: str = "INFO"):
    """统一日志输出，带级别前缀"""
    prefix = {"INFO": "[INFO]", "WARN": "[WARN]", "ERROR": "[ERROR]", "STEP": "[STEP]"}.get(level, "[INFO]")
    print(f"{prefix} {msg}", flush=True)


def log_separator(title: str = ""):
    """打印分隔线"""
    if title:
        print(f"\n{'='*60}", flush=True)
        print(f"  {title}", flush=True)
        print(f"{'='*60}", flush=True)
    else:
        print(f"{'─'*60}", flush=True)


# ============================================================
# 模块二：文件读取
# ============================================================
def read_input_files(cfg: dict) -> tuple:
    """
    读取3个输入文件，自动处理编码问题
    返回：(metadata_df, otutab_df, taxonomy_df)
    """
    log_separator("步骤1：读取输入文件")

    def safe_read(path: str, name: str, index_col=0) -> pd.DataFrame:
        """安全读取制表符分隔txt文件，自动尝试多种编码"""
        abs_path = os.path.abspath(path)
        if not os.path.exists(abs_path):
            log(f"文件不存在：{abs_path}", "ERROR")
            log(f"请检查 CONFIG 中 '{name}' 的路径配置是否正确", "ERROR")
            sys.exit(1)

        for encoding in ["utf-8", "utf-8-sig", "gbk", "gb2312", "latin1"]:
            try:
                df = pd.read_csv(abs_path, sep="\t", index_col=index_col, encoding=encoding)
                log(f"成功读取 {name}：{abs_path}（编码：{encoding}，{df.shape[0]}行 x {df.shape[1]}列）")
                return df
            except UnicodeDecodeError:
                continue
            except Exception as e:
                log(f"读取 {name} 时发生错误：{e}", "ERROR")
                sys.exit(1)

        log(f"无法识别 {name} 的文件编码，请将文件转存为UTF-8格式后重试", "ERROR")
        sys.exit(1)

    metadata = safe_read(cfg["metadata_path"], "metadata.txt", index_col=None)
    otutab   = safe_read(cfg["otutab_path"],   "otutab.txt",   index_col=0)
    taxonomy = safe_read(cfg["taxonomy_path"], "taxonomy.txt", index_col=0)

    # 清理列名前后空白
    metadata.columns = metadata.columns.str.strip()
    otutab.columns   = otutab.columns.str.strip()
    taxonomy.columns = taxonomy.columns.str.strip()

    # ── 自动解析 taxonomy 单列合并格式 ──
    # 支持两种格式：
    #   格式A（7列分开）：kingdom phylum class order family genus species
    #   格式B（单列合并）：k__Bacteria;p__Firmicutes;c__Bacilli;...（分号分隔的完整分类串）
    if taxonomy.shape[1] == 1:
        log("检测到 taxonomy.txt 为单列合并格式（分号分隔），自动解析拆分为7列...")
        taxonomy = _parse_taxonomy_single_column(taxonomy)

    return metadata, otutab, taxonomy


def _parse_taxonomy_single_column(taxonomy: pd.DataFrame) -> pd.DataFrame:
    """
    将单列合并格式的 taxonomy（如 k__Bacteria;p__Firmicutes;...）
    解析拆分为标准7列格式（kingdom/phylum/class/order/family/genus/species）
    """
    PREFIX_TO_COL = {
        "k__": "kingdom", "p__": "phylum", "c__": "class",
        "o__": "order",   "f__": "family", "g__": "genus", "s__": "species",
    }

    def parse_row(tax_str: str) -> dict:
        result = {col: "" for col in TAX_COLUMNS}
        if not isinstance(tax_str, str):
            return result
        for part in tax_str.strip().split(";"):
            part = part.strip()
            for prefix, col in PREFIX_TO_COL.items():
                if part.startswith(prefix):
                    result[col] = part
                    break
        return result

    raw_col = taxonomy.iloc[:, 0]
    parsed = raw_col.apply(parse_row)
    parsed_df = pd.DataFrame(list(parsed), index=taxonomy.index, columns=TAX_COLUMNS)
    log(f"taxonomy 解析完成：{len(parsed_df)} 个OTU，7个分类水平列")
    return parsed_df


# ============================================================
# 模块三：数据校验
# ============================================================
def validate_data(metadata: pd.DataFrame, otutab: pd.DataFrame,
                  taxonomy: pd.DataFrame, cfg: dict) -> tuple:
    """
    完整数据校验流程，返回校验后的有效数据
    返回：(metadata_valid, otutab_valid, taxonomy_valid)
    """
    log_separator("步骤2：数据校验")

    # ── 2.1 校验metadata必要列（大小写不敏感匹配） ──
    log("校验 metadata 核心列...", "STEP")
    col_lower_map = {c.lower(): c for c in metadata.columns}
    for required in ["sampleid", "group"]:
        if required not in col_lower_map:
            log(f"metadata.txt 缺少必要列：'{required}'（大小写不敏感），请检查列名", "ERROR")
            log(f"当前列名为：{list(metadata.columns)}", "ERROR")
            sys.exit(1)
    # 统一重命名为标准小写列名，保证后续逻辑一致
    rename_meta = {}
    if col_lower_map["sampleid"] != "sampleID":
        rename_meta[col_lower_map["sampleid"]] = "sampleID"
    if col_lower_map["group"] != "group":
        rename_meta[col_lower_map["group"]] = "group"
    if rename_meta:
        metadata = metadata.rename(columns=rename_meta)
        log(f"metadata 列名自动规范化：{rename_meta}")

    # 去重sampleID，保留唯一值
    dup_samples = metadata[metadata.duplicated("sampleID")]["sampleID"].tolist()
    if dup_samples:
        log(f"metadata 中 sampleID 存在重复值，将保留第一条记录：{dup_samples}", "WARN")
        metadata = metadata.drop_duplicates(subset="sampleID", keep="first")

    metadata = metadata.set_index("sampleID")
    log(f"metadata 有效样本数：{len(metadata)}，分组信息：{metadata['group'].value_counts().to_dict()}")

    # ── 2.2 校验样本匹配 ──
    log("校验 otutab 样本名 与 metadata sampleID 的匹配情况...", "STEP")
    meta_samples = set(metadata.index)
    otu_samples  = set(otutab.columns)

    only_in_meta = meta_samples - otu_samples
    only_in_otu  = otu_samples  - meta_samples
    common_samples = meta_samples & otu_samples

    if only_in_meta:
        log(f"以下样本在 metadata 中有记录但在 otutab 中无丰度数据，将被剔除：{sorted(only_in_meta)}", "WARN")
    if only_in_otu:
        log(f"以下样本在 otutab 中有丰度数据但在 metadata 中无记录，将被剔除：{sorted(only_in_otu)}", "WARN")
    if not common_samples:
        log("metadata 与 otutab 无任何共有样本，请检查样本名是否一致（区分大小写）", "ERROR")
        sys.exit(1)

    log(f"有效共有样本数：{len(common_samples)}")
    # 按metadata顺序保留有效样本
    valid_samples = [s for s in metadata.index if s in common_samples]
    metadata_valid = metadata.loc[valid_samples]
    otutab_valid   = otutab[valid_samples]

    # ── 2.3 校验taxonomy OTU匹配 ──
    log("校验 taxonomy OTU ID 与 otutab OTU ID 的匹配情况...", "STEP")
    otu_ids_otu  = set(otutab_valid.index)
    otu_ids_tax  = set(taxonomy.index)
    no_tax_otus  = otu_ids_otu - otu_ids_tax

    if no_tax_otus:
        log(f"以下OTU在 otutab 中存在但在 taxonomy 中无注释（共{len(no_tax_otus)}个），将被剔除：", "WARN")
        if len(no_tax_otus) <= 10:
            log(f"  {sorted(no_tax_otus)}", "WARN")
        else:
            log(f"  前10个：{sorted(no_tax_otus)[:10]} ...（共{len(no_tax_otus)}个）", "WARN")

    valid_otu_ids = sorted(otu_ids_otu & otu_ids_tax)
    if not valid_otu_ids:
        log("otutab 与 taxonomy 无任何共有OTU ID，请检查两个文件是否来自同一分析流程", "ERROR")
        sys.exit(1)

    otutab_valid  = otutab_valid.loc[valid_otu_ids]
    taxonomy_valid = taxonomy.loc[valid_otu_ids]
    log(f"有效OTU数：{len(valid_otu_ids)}")

    # ── 2.4 校验分类水平配置 ──
    log(f"校验分类水平配置：'{cfg['tax_level']}'...", "STEP")
    if cfg["tax_level"] not in TAX_LEVEL_MAP:
        log(f"不支持的分类水平：'{cfg['tax_level']}'，可选值为：{list(TAX_LEVEL_MAP.keys())}", "ERROR")
        sys.exit(1)

    # ── 2.5 校验taxonomy列名 ──
    # 标准化taxonomy列名（允许首字母大小写不同）
    tax_col_lower = {c.lower(): c for c in taxonomy_valid.columns}
    rename_map = {}
    for std_col in TAX_COLUMNS:
        if std_col in taxonomy_valid.columns:
            continue
        elif std_col in tax_col_lower:
            rename_map[tax_col_lower[std_col]] = std_col
        else:
            log(f"taxonomy.txt 缺少分类水平列：'{std_col}'，请检查列名格式", "WARN")

    if rename_map:
        taxonomy_valid = taxonomy_valid.rename(columns=rename_map)
        log(f"列名自动映射：{rename_map}", "WARN")

    # ── 2.6 校验选定分类水平的注释缺失 ──
    target_col = TAX_LEVEL_MAP[cfg["tax_level"]]["col"]
    if target_col not in taxonomy_valid.columns:
        log(f"taxonomy.txt 中缺少选定分类水平的列：'{target_col}'", "ERROR")
        sys.exit(1)

    # 过滤该分类水平注释缺失的OTU（空值、纯前缀、"unclassified"等）
    prefix = TAX_LEVEL_MAP[cfg["tax_level"]]["prefix"]
    tax_col_data = taxonomy_valid[target_col].astype(str).str.strip()
    invalid_mask = (
        tax_col_data.isin(["", "nan", prefix, prefix.rstrip("_"), "unclassified", "Unclassified"])
        | tax_col_data.str.match(r"^[a-z]__\s*$")
    )
    n_invalid = invalid_mask.sum()
    if n_invalid > 0:
        log(f"分类水平 '{cfg['tax_level']}' 中有 {n_invalid} 个OTU注释缺失，将被剔除", "WARN")
        taxonomy_valid = taxonomy_valid[~invalid_mask]
        otutab_valid   = otutab_valid.loc[taxonomy_valid.index]

    log(f"校验完成：有效样本 {len(valid_samples)} 个，有效OTU {len(taxonomy_valid)} 个")
    log_separator()

    return metadata_valid, otutab_valid, taxonomy_valid


# ============================================================
# 模块四：构建分类名称（含上级分类注释）
# ============================================================
def build_taxon_name(row: pd.Series, tax_level: str, higher_levels: list) -> str:
    """
    构建单个OTU的分类名称，按 higher_levels 指定的上级水平拼接前缀
    - higher_levels=[]：仅输出当前水平名称，如 g__Blautia
    - higher_levels=["phylum"]：输出 p__Firmicutes;g__Blautia
    - higher_levels=["phylum","family"]：输出 p__Firmicutes;f__Lachnospiraceae;g__Blautia
    只有比 tax_level 更高的水平才会被纳入（低于或等于当前水平的自动跳过）
    """
    level_info  = TAX_LEVEL_MAP[tax_level]
    target_col  = level_info["col"]
    target_val  = str(row.get(target_col, "")).strip()
    current_idx = level_info["idx"]

    # 过滤掉不合法的上级水平（等于或低于当前水平）
    valid_higher = [
        lv for lv in higher_levels
        if lv in TAX_LEVEL_MAP and TAX_LEVEL_MAP[lv]["idx"] < current_idx
    ]

    if not valid_higher:
        return target_val

    parts = []
    for higher_col in valid_higher:
        val = str(row.get(higher_col, "")).strip()
        if val and val not in ("nan", "", TAX_PREFIX_MAP.get(higher_col, ""), "unclassified"):
            parts.append(val)
        else:
            parts.append(f"{TAX_PREFIX_MAP.get(higher_col, '')}unclassified")

    parts.append(target_val)
    return ";".join(parts)


# ============================================================
# 模块五：生成原始丰度矩阵（OTU → 分类水平汇总）
# ============================================================
def build_abundance_matrix(otutab: pd.DataFrame, taxonomy: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """
    将OTU丰度表按选定分类水平汇总，生成原始丰度矩阵
    行=物种分类，列=样本
    """
    log_separator("步骤3：生成原始丰度矩阵")
    tax_level     = cfg["tax_level"]
    higher_levels = cfg.get("higher_taxonomy_levels", [])

    # 校验 higher_taxonomy_levels 中的无效条目并提示
    invalid_lvs = [lv for lv in higher_levels if lv not in TAX_LEVEL_MAP]
    if invalid_lvs:
        log(f"higher_taxonomy_levels 中包含无效水平：{invalid_lvs}，将被忽略", "WARN")
    skipped_lvs = [
        lv for lv in higher_levels
        if lv in TAX_LEVEL_MAP and TAX_LEVEL_MAP[lv]["idx"] >= TAX_LEVEL_MAP[tax_level]["idx"]
    ]
    if skipped_lvs:
        log(f"higher_taxonomy_levels 中以下水平不高于当前水平 '{tax_level}'，将被忽略：{skipped_lvs}", "WARN")

    valid_higher = [
        lv for lv in higher_levels
        if lv in TAX_LEVEL_MAP and TAX_LEVEL_MAP[lv]["idx"] < TAX_LEVEL_MAP[tax_level]["idx"]
    ]
    if valid_higher:
        log(f"上级分类注释：将附加 {valid_higher} 水平信息")
    else:
        log("上级分类注释：未指定，仅输出当前水平名称")

    # 构建每个OTU的分类标签
    tax_labels = taxonomy.apply(
        lambda row: build_taxon_name(row, tax_level, valid_higher), axis=1
    )

    # 将OTU丰度表与分类标签合并，按分类标签汇总（sum）
    otutab_labeled = otutab.copy()
    otutab_labeled.insert(0, "_taxon_", tax_labels)
    abundance_matrix = otutab_labeled.groupby("_taxon_").sum()

    log(f"按 '{tax_level}' 水平汇总完成：{len(abundance_matrix)} 个分类单元，{len(abundance_matrix.columns)} 个样本")
    return abundance_matrix


# ============================================================
# 模块六：转换为相对丰度
# ============================================================
def to_relative_abundance(abundance_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    将原始丰度矩阵转换为相对丰度（每个样本内归一化到1）
    遵循：每个物种的相对丰度 = 该物种reads / 该样本总reads
    """
    log_separator("步骤4：转换为相对丰度")
    sample_totals = abundance_matrix.sum(axis=0)

    zero_total_samples = sample_totals[sample_totals == 0].index.tolist()
    if zero_total_samples:
        log(f"以下样本总reads为0，将在相对丰度中保持为0（请检查数据完整性）：{zero_total_samples}", "WARN")

    rel_abundance = abundance_matrix.div(sample_totals.replace(0, np.nan), axis=1).fillna(0)
    log(f"相对丰度转换完成，各样本总丰度均为1.0（reads为0的样本除外）")
    return rel_abundance


# ============================================================
# 模块七：按分组聚合（可选）
# ============================================================
def aggregate_by_group(rel_abundance: pd.DataFrame, metadata: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """
    若 aggregate_by="group"，则按分组计算组内相对丰度总和
    """
    if cfg["aggregate_by"] == "sample":
        log_separator("步骤5：聚合维度 = 按样本，跳过分组聚合")
        return rel_abundance

    log_separator("步骤5：按分组计算相对丰度总和")
    # 转置：行=样本，列=物种
    rel_t = rel_abundance.T
    rel_t.index.name = "sampleID"
    rel_t = rel_t.join(metadata[["group"]], how="left")

    # 按group分组求总和
    group_sum = rel_t.groupby("group").sum()

    log(f"按分组聚合完成：{len(group_sum)} 个分组")
    for grp, cnt in metadata["group"].value_counts().items():
        log(f"  分组 '{grp}'：{cnt} 个样本")

    # 转回：行=物种，列=分组
    return group_sum.T


# ============================================================
# 模块八：低丰度过滤
# ============================================================
def filter_low_abundance(matrix: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """
    过滤所有样本/分组中平均相对丰度低于阈值的物种
    注意：过滤基于相对丰度（标准化前），行=物种，列=样本/分组
    """
    log_separator("步骤6（前置）：低丰度物种过滤")
    threshold = cfg.get("abundance_threshold") or 0

    if not threshold:
        log("低丰度过滤已关闭，保留全部物种")
        return matrix

    n_before = len(matrix)
    row_means = matrix.mean(axis=1)
    keep_mask = row_means >= threshold
    matrix_filtered = matrix[keep_mask]
    n_after = len(matrix_filtered)

    log(f"过滤阈值：平均相对丰度 >= {threshold:.4f}（{threshold*100:.2f}%）")
    log(f"过滤前物种数：{n_before}，过滤后物种数：{n_after}，共过滤 {n_before - n_after} 个低丰度物种")

    if n_after == 0:
        log("警告：过滤后无任何物种保留！请降低 abundance_threshold 阈值后重试", "ERROR")
        sys.exit(1)

    return matrix_filtered


# ============================================================
# 模块九：标准化处理
# ============================================================
def normalize_matrix(matrix: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """
    对丰度矩阵进行标准化处理
    - zscore：按物种维度（行）Z-score标准化，符合16S行业规范
    - minmax：按物种维度（行）Min-Max归一化
    - none：不做任何标准化
    """
    log_separator("步骤7：标准化处理")
    method = cfg["normalization"]

    if method == "none":
        log("标准化方法：无标准化，直接使用相对丰度矩阵")
        return matrix

    elif method == "zscore":
        log("标准化方法：Z-score标准化（按物种维度，行标准化）")
        # scipy.stats.zscore 按行计算，ddof=1使用样本标准差
        normalized = matrix.apply(
            lambda row: stats.zscore(row, ddof=1) if row.std() > 0 else row * 0,
            axis=1,
            result_type="expand"
        )
        normalized.columns = matrix.columns
        normalized.index   = matrix.index
        log("Z-score标准化完成（方差为0的物种对应行全部置0）")
        return normalized

    elif method == "minmax":
        log("标准化方法：Min-Max标准化（按物种维度，行标准化）")
        row_min = matrix.min(axis=1)
        row_max = matrix.max(axis=1)
        row_range = row_max - row_min
        # 避免除以0（范围为0时，结果为0）
        normalized = matrix.sub(row_min, axis=0).div(row_range.replace(0, np.nan), axis=0).fillna(0)
        log("Min-Max标准化完成")
        return normalized

    else:
        log(f"不支持的标准化方法：'{method}'，可选值为：none | zscore | minmax", "ERROR")
        sys.exit(1)


# ============================================================
# 模块十：调整输出方向并导出CSV
# ============================================================
def _split_taxon_index(matrix: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """
    将拼接在行索引中的分类信息拆分为独立的分类列，插入数值列之前。
    每列只保留 __ 后面的内容（去掉 p__/c__/g__ 等前缀）。
    仅在 orientation="vertical"（行=物种）时调用。

    示例：
      原索引：p__Firmicutes;f__Lachnospiraceae;g__Blautia
      → 新增列 phylum=Firmicutes, family=Lachnospiraceae, genus=Blautia
    """
    tax_level     = cfg["tax_level"]
    higher_levels = cfg.get("higher_taxonomy_levels", [])

    # 按实际写入顺序确定所有分类列（顺序与构建时一致）
    current_idx  = TAX_LEVEL_MAP[tax_level]["idx"]
    valid_higher = [
        lv for lv in higher_levels
        if lv in TAX_LEVEL_MAP and TAX_LEVEL_MAP[lv]["idx"] < current_idx
    ]
    all_levels = valid_higher + [tax_level]  # 上级在前，当前水平在后

    # 前缀 → 水平名映射，用于从拼接串中识别各段
    prefix_to_level = {TAX_PREFIX_MAP[lv]: lv for lv in all_levels}

    def parse_index(taxon_str: str) -> dict:
        """将单条拼接索引解析为 {水平名: 去前缀值} 字典"""
        result = {lv: "" for lv in all_levels}
        if not isinstance(taxon_str, str):
            return result
        for part in taxon_str.split(";"):
            part = part.strip()
            for prefix, lv in prefix_to_level.items():
                if part.startswith(prefix):
                    result[lv] = part[len(prefix):]  # 只保留 __ 后内容
                    break
            else:
                # 无前缀时按当前水平处理（兼容仅有当前水平名的情况）
                if not any(part.startswith(p) for p in prefix_to_level):
                    result[tax_level] = part
        return result

    parsed  = matrix.index.map(parse_index)
    tax_df  = pd.DataFrame(list(parsed), columns=all_levels)  # 用整数索引，避免concat错位

    # 拼接：分类列在左，数值列在右，均使用整数索引对齐
    data_df = matrix.reset_index(drop=True)
    output_df = pd.concat([tax_df, data_df], axis=1)
    output_df.index = range(len(output_df))

    log(f"分类列已拆分为独立列：{all_levels}（内容已去除分类前缀）")
    return output_df


def export_matrix(matrix: pd.DataFrame, cfg: dict, stats_info: dict):
    """
    调整矩阵方向、拆分分类列并导出为CSV文件
    同时打印完整的运行统计信息
    """
    log_separator("步骤8：调整输出方向并导出CSV")

    orientation = cfg["orientation"]
    if orientation == "vertical":
        log("输出方向：纵向（行=物种，列=分类信息列+样本/分组列）")
        output_df = _split_taxon_index(matrix, cfg)
    elif orientation == "horizontal":
        output_df = matrix.T
        output_df.index.name = "sample_group"
        log("输出方向：横向（行=样本/分组，列=物种）")
    else:
        log(f"不支持的输出方向：'{orientation}'，可选值为：vertical | horizontal", "ERROR")
        sys.exit(1)

    # 导出CSV
    output_path = os.path.abspath(cfg["output_path"])
    output_df.to_csv(output_path, encoding="utf-8-sig", index=(orientation == "horizontal"))

    log(f"CSV文件已成功导出至：{output_path}")
    log(f"矩阵维度：{output_df.shape[0]} 行 x {output_df.shape[1]} 列")

    # 打印完整运行统计
    log_separator("运行统计摘要")
    for k, v in stats_info.items():
        log(f"  {k}：{v}")
    log_separator()


# ============================================================
# 主流程
# ============================================================
def main():
    log_separator("16S扩增子 - 物种丰度热图矩阵生成工具 v1.0.0")
    log(f"当前工作目录：{os.getcwd()}")
    log(f"分类水平：{CONFIG['tax_level']} | 聚合方式：{CONFIG['aggregate_by']} | "
        f"标准化：{CONFIG['normalization']} | 方向：{CONFIG['orientation']}")
    log_separator()

    # Step 1: 读取文件
    metadata, otutab, taxonomy = read_input_files(CONFIG)

    # Step 2: 数据校验（返回清洗后的有效数据）
    metadata, otutab, taxonomy = validate_data(metadata, otutab, taxonomy, CONFIG)

    # Step 3: 生成原始丰度矩阵（OTU → 分类水平汇总）
    abundance_matrix = build_abundance_matrix(otutab, taxonomy, CONFIG)

    # Step 4: 转换为相对丰度（每个样本归一化到1）
    rel_abundance = to_relative_abundance(abundance_matrix)

    # Step 5: 按分组聚合（可选）
    aggregated = aggregate_by_group(rel_abundance, metadata, CONFIG)

    # Step 6: 低丰度过滤（基于相对丰度，在标准化前过滤）
    filtered = filter_low_abundance(aggregated, CONFIG)

    # Step 7: 标准化处理
    normalized = normalize_matrix(filtered, CONFIG)

    # Step 8: 导出CSV，打印统计信息
    stats_info = {
        "有效样本数":   len(metadata),
        "有效OTU数":    len(otutab),
        f"{CONFIG['tax_level']}水平分类单元数（过滤前）": len(aggregated),
        f"{CONFIG['tax_level']}水平分类单元数（过滤后）": len(filtered),
        "低丰度过滤阈值": f"{CONFIG['abundance_threshold']*100:.2f}%" if CONFIG['abundance_threshold'] else "未启用",
        "标准化方法":   CONFIG["normalization"],
        "聚合维度":     CONFIG["aggregate_by"],
        "输出方向":     CONFIG["orientation"],
        "输出文件":     os.path.abspath(CONFIG["output_path"]),
    }
    export_matrix(normalized, CONFIG, stats_info)

    log("全部流程执行完成！")


# ============================================================
# 程序入口
# ============================================================
if __name__ == "__main__":
    main()


# ============================================================
# 使用说明（Usage Documentation）
# ============================================================
"""
==============================================================
    16S扩增子 物种丰度热图矩阵生成工具 - 使用说明
==============================================================

【一、环境配置】
    需要Python 3.8+，安装以下依赖包：
        pip install pandas numpy scipy

    或使用conda：
        conda install pandas numpy scipy

【二、参数修改方法】
    所有配置项均在代码顶部的 CONFIG 字典中修改，常用修改示例：

    1. 切换分类水平（如改为门水平）：
       "tax_level": "phylum"
       可选：phylum | class | order | family | genus

    2. 按分组输出（计算组内相对丰度总和）：
       "aggregate_by": "group"
       可选：sample | group

    3. 关闭标准化（输出原始相对丰度）：
       "normalization": "none"
       可选：none | zscore | minmax

    4. 关闭低丰度过滤：
       "abundance_threshold": 0

    5. 更改输入文件路径：
       "metadata_path": "path/to/your/metadata.txt"

    6. 关闭上级分类注释：
       "add_higher_taxonomy": False

【三、输入文件格式要求】

    1. metadata.txt（制表符分隔）
       必须包含 sampleID 列和 group 列，其余列可自由扩展：
       sampleID    group    age    sex
       Sample1     Control  25     M
       Sample2     Control  30     F
       Sample3     Treatment 28    M

    2. otutab.txt（制表符分隔）
       第一列为OTU ID（无列名或列名为#OTU ID均可），其余列为样本名：
       #OTU ID    Sample1    Sample2    Sample3
       OTU001     100        200        50
       OTU002     0          300        120

    3. taxonomy.txt（制表符分隔）
       第一列为OTU ID，列名为7个分类水平（kingdom/phylum/class/order/family/genus/species）：
       #OTU ID    kingdom      phylum          class           ...
       OTU001     k__Bacteria  p__Firmicutes   c__Bacilli      ...
       注意：列名大小写不敏感，程序会自动识别

【四、输出文件说明】
    输出为UTF-8编码的CSV文件，第一列为物种/样本名，其余列为样本/分组名。
    - 标准化方法为none时：值为相对丰度（0~1之间的小数）
    - 标准化方法为zscore时：值为Z-score（可为负数，适合热图配色）
    - 标准化方法为minmax时：值在0~1之间（另一种标准化方式）

    兼容以下热图绘制工具：
    ✓ TBtools（直接导入CSV绘制热图）
    ✓ R pheatmap（read.csv() 直接读取，转为矩阵后绘图）
    ✓ GraphPad Prism（导入CSV后选择热图）
    ✓ Origin（导入CSV，插入热图）
    ✓ Python seaborn/matplotlib（pd.read_csv() 直接读取）

【五、常见问题排查】

    Q1：报错 "文件不存在"
    A：检查文件路径是否正确，可使用绝对路径，注意Windows路径反斜杠需写为
       双反斜杠（\\）或改用正斜杠（/）

    Q2：报错 "metadata 与 otutab 无任何共有样本"
    A：检查otutab的列名和metadata的sampleID是否完全一致（区分大小写和空格）

    Q3：报错 "过滤后无任何物种保留"
    A：降低 abundance_threshold 阈值，如改为0.0001（0.01%），或设为0关闭过滤

    Q4：导入Excel后中文乱码
    A：输出文件使用utf-8-sig编码，Excel可正常识别；若仍乱码请用Excel"数据→
       导入文本"功能，选择UTF-8编码

    Q5：taxonomy列名不标准（如Phylum而非phylum）
    A：程序会自动识别大小写，但如果列名完全不同（如使用中文）则需手动统一

    Q6：热图中出现大量NaN
    A：通常因Z-score标准化时某物种在所有样本中丰度相同（方差为0），程序已自动
       将这些行置0，不影响绘图

【六、引用与合规】
    本工具计算逻辑符合QIIME2和USEARCH的行业标准：
    - 相对丰度：按样本归一化，与QIIME2 feature-table relative-frequency一致
    - Z-score：按物种（行）维度计算，使用样本标准差（ddof=1）
    - 低丰度过滤：基于各样本平均相对丰度，与QIIME2 filter-features一致

==============================================================
"""
