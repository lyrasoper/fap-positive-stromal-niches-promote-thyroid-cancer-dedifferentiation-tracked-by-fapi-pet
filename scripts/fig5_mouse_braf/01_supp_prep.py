# Cancer Research submission - figure code release
# Builds: Supplementary Fig. S21 (mouse single-cell source characteristics
#          and bulk-cohort validation). Internal names use "supp_fig6".

from __future__ import annotations
import os

from io import StringIO
from pathlib import Path
import shutil
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from PIL import Image


PROJECT_DIR = Path(os.environ.get("PROJECT_ROOT", "<PROJECT_ROOT>"))
FIG3_DATA_DIR = PROJECT_DIR / "data" / "fig3"
SUPP_FIG_DIR = PROJECT_DIR / "outputs" / "supplementary_figures"
SUPP_TABLE_DIR = PROJECT_DIR / "outputs" / "supplementary_tables"
SOURCE_DATA_DIR = SUPP_FIG_DIR / "source_data" / "supp_fig6"

BULK_PROJECT_DIR = Path((os.environ.get("EXTERNAL_DATA", "<EXTERNAL_DATA>") + "/mouse_bulk"))
BULK_DATA_DIR = BULK_PROJECT_DIR / "mTC"
BULK_SCRIPT_DIR = BULK_PROJECT_DIR / "scripts"

FULL_OBJ_PATH = FIG3_DATA_DIR / "epi_fib_imm_endo_seurat_object.rds"
FIB_OBJ_PATH = FIG3_DATA_DIR / "fibroblast_seurat_object.rds"

GSE118022_HEATMAP_PATH = BULK_PROJECT_DIR / "results" / "dataset_support_fap_dediff" / "gse118022_marker_heatmap.png"

SUPP_FIG6_PNG = SUPP_FIG_DIR / "Supplementary_Fig_6_Source_overview_and_extended_transcriptomic_validation_for_the_time_resolved_mouse_model.png"
SUPP_FIG6_PDF = SUPP_FIG_DIR / "Supplementary_Fig_6_Source_overview_and_extended_transcriptomic_validation_for_the_time_resolved_mouse_model.pdf"
TABLE7_PATH = SUPP_TABLE_DIR / "Supplementary_Table_7_mouse_scRNA_and_bulk_mTC_cohorts.xlsx"
TABLE1_PATH = SUPP_TABLE_DIR / "Supplementary_Table_1_Public_thyroid_cancer_scRNA_seq_cohort_metadata.xlsx"
TABLE2_PATH = SUPP_TABLE_DIR / "Supplementary_Table_2_Sample_level_major_lineage_counts_and_proportions.xlsx"
TABLE3_PATH = SUPP_TABLE_DIR / "Supplementary_Table_3_Fibroblast_CAF_reclustering_and_subtype_annotation_summary.xlsx"
TABLE4_PATH = SUPP_TABLE_DIR / "Supplementary_Table_4_Differentially_expressed_genes_in_FAP_positive_infCAF_vs_other_CAF_populations.xlsx"
TABLE5_PATH = SUPP_TABLE_DIR / "Supplementary_Table_5_Gene_sets_and_statistical_summary_for_Fig1_and_Supplementary_Figs_1_3.xlsx"
TABLE6_PATH = SUPP_TABLE_DIR / "Supplementary_Table_6_Spatial_transcriptomics_cohort_metadata_and_spot_cluster_summary.xlsx"
TABLE7_PATH = SUPP_TABLE_DIR / "Supplementary_Table_7_Mouse_scRNA_seq_and_bulk_mTC_cohort_summary.xlsx"
COMBINED_1_7_PATH = SUPP_TABLE_DIR / "Supplementary_Tables_1_7_Combined.xlsx"

TIMEPOINT_ORDER = {
    "mPTC_1month": 1,
    "mPTC_2month": 2,
    "mPTC_4month": 4,
}
TIMEPOINT_LABELS = {
    "mPTC_1month": "1 month",
    "mPTC_2month": "2 months",
    "mPTC_4month": "4 months",
}

RADAR_GROUPS = ["WT", "FTC", "ATC"]
RADAR_COLORS = {
    "WT": "#E07A6A",
    "FTC": "#7B97F3",
    "ATC": "#63A78A",
}
RADAR_GENES = [
    "Tpo",
    "Tg",
    "Serpine1",
    "Ctgf",
    "Tagln",
    "Acta2",
    "Lox",
    "Lum",
    "Dcn",
    "Postn",
    "Sparc",
    "Fn1",
    "Col3a1",
    "Col1a2",
    "Col1a1",
    "Fap",
    "Dio1",
    "Nkx2-1",
    "Tshr",
    "Pax8",
    "Slc5a5",
]

if not (BULK_SCRIPT_DIR / "integrate_thyroid_markers.py").exists():
    raise FileNotFoundError(
        f"Missing {BULK_SCRIPT_DIR / 'integrate_thyroid_markers.py'} — "
        "external bulk-mTC helper module (set EXTERNAL_DATA); see README."
    )
sys.path.insert(0, str(BULK_SCRIPT_DIR))
from integrate_thyroid_markers import (  # noqa: E402
    align_metadata,
    load_gse118022_metadata,
    load_gse30427_metadata,
    load_gse55933_metadata,
    load_matrix,
)


def run_r_csv(code: str) -> pd.DataFrame:
    if shutil.which("Rscript") is None:
        raise RuntimeError(
            "Rscript executable not found on PATH; install R to extract scRNA "
            "summaries from the Seurat objects."
        )
    completed = subprocess.run(
        ["Rscript", "-e", code],
        check=True,
        capture_output=True,
        text=True,
    )
    return pd.read_csv(StringIO(completed.stdout))


def fetch_scRNA_summaries() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for p in (FULL_OBJ_PATH, FIB_OBJ_PATH):
        if not p.exists():
            raise FileNotFoundError(
                f"Missing {p} — Seurat object produced by an upstream Fig. 3 "
                "mouse scRNA processing step; see README."
            )
    full_timepoint = run_r_csv(
        f"""
        obj <- readRDS('{FULL_OBJ_PATH.as_posix()}');
        md <- obj@meta.data;
        out <- as.data.frame(table(md$orig.ident));
        names(out) <- c('timepoint', 'n_cells');
        write.csv(out, stdout(), row.names = FALSE)
        """
    )
    full_celltype = run_r_csv(
        f"""
        obj <- readRDS('{FULL_OBJ_PATH.as_posix()}');
        md <- obj@meta.data;
        out <- as.data.frame(table(md$orig.ident, md$celltype));
        names(out) <- c('timepoint', 'celltype', 'n_cells');
        write.csv(out, stdout(), row.names = FALSE)
        """
    )
    fib_subtype = run_r_csv(
        f"""
        obj <- readRDS('{FIB_OBJ_PATH.as_posix()}');
        md <- obj@meta.data;
        out <- as.data.frame(table(md$orig.ident, md$Fib_clu1));
        names(out) <- c('timepoint', 'caf_state', 'n_cells');
        write.csv(out, stdout(), row.names = FALSE)
        """
    )

    for frame in (full_timepoint, full_celltype, fib_subtype):
        frame["timepoint_order"] = frame["timepoint"].map(TIMEPOINT_ORDER)
        frame["timepoint_label"] = frame["timepoint"].map(TIMEPOINT_LABELS)
        frame.sort_values(["timepoint_order"], inplace=True)

    return full_timepoint, full_celltype, fib_subtype


def format_counts(series: pd.Series) -> str:
    ordered = series.sort_index()
    return ", ".join(f"{group}={int(count)}" for group, count in ordered.items())


def build_bulk_summary() -> pd.DataFrame:
    bulk_inputs = [
        BULK_DATA_DIR / "GSE55933_expr_matrix.csv",
        BULK_DATA_DIR / "GSE118022_groupinf.csv",
        BULK_DATA_DIR / "GSE30427_expr.csv",
        BULK_DATA_DIR / "GSE55933_phenotype_info.csv",
        BULK_DATA_DIR / "GSE118022_expr.csv",
        BULK_DATA_DIR / "GSE30427_groupinf.csv",
    ]
    for p in bulk_inputs:
        if not p.exists():
            raise FileNotFoundError(
                f"Missing {p} — external bulk mTC cohort file (set EXTERNAL_DATA); "
                "see README."
            )
    gse55933_matrix = load_matrix(BULK_DATA_DIR / "GSE55933_expr_matrix.csv")
    gse118022_matrix = load_matrix(BULK_DATA_DIR / "GSE118022_groupinf.csv")
    gse30427_matrix = load_matrix(BULK_DATA_DIR / "GSE30427_expr.csv")

    gse55933_meta = align_metadata(
        load_gse55933_metadata(BULK_DATA_DIR / "GSE55933_phenotype_info.csv"),
        gse55933_matrix,
    )
    gse118022_meta = align_metadata(
        load_gse118022_metadata(BULK_DATA_DIR / "GSE118022_expr.csv"),
        gse118022_matrix,
    )
    gse30427_meta = align_metadata(
        load_gse30427_metadata(BULK_DATA_DIR / "GSE30427_groupinf.csv"),
        gse30427_matrix,
    )

    rows = [
        {
            "cohort": "GSE55933",
            "species": "Mouse",
            "total_samples": int(gse55933_meta.shape[0]),
            "all_groups": format_counts(gse55933_meta["group"].value_counts()),
            "displayed_groups": "mPTC, mATC",
            "analytical_role": "Independent bulk reference cohort for Braf-driven mPTC/mATC expression patterns.",
            "source_files": "GSE55933_expr_matrix.csv; GSE55933_phenotype_info.csv",
        },
        {
            "cohort": "GSE118022",
            "species": "Mouse",
            "total_samples": int(gse118022_meta.shape[0]),
            "all_groups": format_counts(gse118022_meta["group"].value_counts()),
            "displayed_groups": "mNT, mPTC, mATC, remATC",
            "analytical_role": "External module-score validation used in Fig. 3I and Supplementary Fig. 6C.",
            "source_files": "GSE118022_groupinf.csv; GSE118022_expr.csv",
        },
        {
            "cohort": "GSE30427",
            "species": "Mouse",
            "total_samples": int(gse30427_meta.shape[0]),
            "all_groups": format_counts(gse30427_meta["group"].value_counts()),
            "displayed_groups": "WT, FTC, ATC (display subset); full cohort also includes Pten_KO, p53_KO, and Pten_Kras",
            "analytical_role": "Independent validation of differentiation loss versus FAP/ECM activation in Supplementary Fig. 6D.",
            "source_files": "GSE30427_expr.csv; GSE30427_groupinf.csv",
        },
    ]
    return pd.DataFrame(rows)


def build_scRNA_overview(full_timepoint: pd.DataFrame, fib_subtype: pd.DataFrame) -> pd.DataFrame:
    fib_by_timepoint = (
        fib_subtype.groupby(["timepoint", "timepoint_label", "timepoint_order"], as_index=False)["n_cells"]
        .sum()
        .sort_values("timepoint_order")
    )
    merged = full_timepoint.merge(
        fib_by_timepoint[["timepoint", "n_cells"]].rename(columns={"n_cells": "fibroblast_caf_cells"}),
        on="timepoint",
        how="left",
    )
    merged["timepoint"] = merged["timepoint"].map(TIMEPOINT_LABELS)
    merged.rename(columns={"n_cells": "all_cells"}, inplace=True)
    return merged[["timepoint", "all_cells", "fibroblast_caf_cells"]]


def write_table7(
    full_timepoint: pd.DataFrame,
    full_celltype: pd.DataFrame,
    fib_subtype: pd.DataFrame,
    bulk_summary: pd.DataFrame,
) -> None:
    SUPP_TABLE_DIR.mkdir(parents=True, exist_ok=True)

    sc_overview = pd.DataFrame(
        [
            {
                "cohort": "Time-resolved mouse scRNA-seq cohort",
                "model": "Braf^LSL-V600E; Tpo-Cre",
                "timepoints": "mPTC_1month, mPTC_2month, mPTC_4month",
                "full_object_file": FULL_OBJ_PATH.name,
                "full_object_cells": int(full_timepoint["n_cells"].sum()),
                "fibroblast_object_file": FIB_OBJ_PATH.name,
                "fibroblast_object_cells": int(fib_subtype["n_cells"].sum()),
                "analytical_role": "Source of the mouse single-cell analyses shown in Fig. 3A-H.",
            }
        ]
    )

    full_timepoint_export = full_timepoint.copy()
    full_timepoint_export["timepoint"] = full_timepoint_export["timepoint"].map(TIMEPOINT_LABELS)
    full_timepoint_export = full_timepoint_export[["timepoint", "n_cells"]].rename(columns={"n_cells": "all_cells"})

    full_celltype_export = full_celltype.copy()
    full_celltype_export["timepoint"] = full_celltype_export["timepoint"].map(TIMEPOINT_LABELS)
    full_celltype_export = full_celltype_export[["timepoint", "celltype", "n_cells"]]

    fib_subtype_export = fib_subtype.copy()
    fib_subtype_export["timepoint"] = fib_subtype_export["timepoint"].map(TIMEPOINT_LABELS)
    fib_subtype_export = fib_subtype_export[["timepoint", "caf_state", "n_cells"]]

    with pd.ExcelWriter(TABLE7_PATH, engine="openpyxl") as writer:
        sc_overview.to_excel(writer, sheet_name="mouse_scRNA_overview", index=False)
        full_timepoint_export.to_excel(writer, sheet_name="mouse_scRNA_timepoints", index=False)
        full_celltype_export.to_excel(writer, sheet_name="mouse_scRNA_celltypes", index=False)
        fib_subtype_export.to_excel(writer, sheet_name="mouse_fibroblast_states", index=False)
        bulk_summary.to_excel(writer, sheet_name="bulk_mTC_cohorts", index=False)


def safe_sheet_name(prefix: str, name: str) -> str:
    raw = f"{prefix}_{name}".replace(" ", "_")
    return raw[:31]


def rebuild_combined_workbook() -> None:
    table_paths = [
        ("T1", TABLE1_PATH),
        ("T2", TABLE2_PATH),
        ("T3", TABLE3_PATH),
        ("T4", TABLE4_PATH),
        ("T5", TABLE5_PATH),
        ("T6", TABLE6_PATH),
        ("T7", TABLE7_PATH),
    ]

    with pd.ExcelWriter(COMBINED_1_7_PATH, engine="openpyxl") as writer:
        for prefix, table_path in table_paths:
            if not table_path.exists():
                continue
            excel = pd.ExcelFile(table_path)
            for sheet_name in excel.sheet_names:
                frame = pd.read_excel(table_path, sheet_name=sheet_name)
                frame.to_excel(writer, sheet_name=safe_sheet_name(prefix, sheet_name), index=False)


def draw_table_panel(ax: plt.Axes, title: str, dataframe: pd.DataFrame, footnote: str | None = None) -> None:
    ax.axis("off")
    ax.text(0, 1.08, title, transform=ax.transAxes, fontsize=12.5, fontweight="bold", va="bottom")

    cell_text = dataframe.astype(str).values.tolist()
    table = ax.table(
        cellText=cell_text,
        colLabels=list(dataframe.columns),
        loc="upper left",
        cellLoc="left",
        colLoc="left",
        bbox=[0, 0.10, 1, 0.85],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.3)

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#C8C8C8")
        if row == 0:
            cell.set_facecolor("#EDEDED")
            cell.set_text_props(weight="bold")
        else:
            cell.set_facecolor("white")

    if footnote:
        ax.text(0, 0.02, footnote, transform=ax.transAxes, fontsize=8.2, color="#4A4A4A", va="bottom")


def add_panel_letter(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.10, 1.12, letter, transform=ax.transAxes, fontsize=18, fontweight="bold", va="top")


def build_radar_panel() -> tuple[pd.DataFrame, pd.DataFrame]:
    for p in (BULK_DATA_DIR / "GSE30427_expr.csv", BULK_DATA_DIR / "GSE30427_groupinf.csv"):
        if not p.exists():
            raise FileNotFoundError(
                f"Missing {p} — external bulk mTC cohort file (set EXTERNAL_DATA); "
                "see README."
            )
    matrix = load_matrix(BULK_DATA_DIR / "GSE30427_expr.csv")
    metadata = align_metadata(
        load_gse30427_metadata(BULK_DATA_DIR / "GSE30427_groupinf.csv"),
        matrix,
    )
    metadata = metadata.loc[metadata["group"].isin(RADAR_GROUPS)].copy()
    metadata["group"] = pd.Categorical(metadata["group"], categories=RADAR_GROUPS, ordered=True)
    metadata.sort_values("group", inplace=True)

    group_means = pd.DataFrame(index=RADAR_GENES, columns=RADAR_GROUPS, dtype=float)
    for group in RADAR_GROUPS:
        sample_ids = metadata.loc[metadata["group"] == group, "sample_id"].tolist()
        group_means[group] = matrix.reindex(RADAR_GENES)[sample_ids].mean(axis=1)

    scaled = group_means.copy()
    for gene in scaled.index:
        values = group_means.loc[gene].astype(float)
        vmin = np.nanmin(values.to_numpy())
        vmax = np.nanmax(values.to_numpy())
        if np.isclose(vmax, vmin):
            scaled.loc[gene] = 0.0
        else:
            scaled.loc[gene] = (values - vmin) / (vmax - vmin) * 100.0

    return group_means.reset_index(names="gene"), scaled.reset_index(names="gene")


def draw_radar(ax: plt.Axes, scaled_radar: pd.DataFrame) -> None:
    genes = scaled_radar["gene"].tolist()
    angles = np.linspace(0, 2 * np.pi, len(genes), endpoint=False)
    angles = np.concatenate([angles, [angles[0]]])

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(genes, fontsize=8.3)
    ax.set_ylim(0, 100)
    ax.set_yticks([25, 50, 75, 100])
    ax.set_yticklabels(["25", "50", "75", "100 (%)"], fontsize=8)
    ax.grid(color="#CFCFCF", alpha=0.8)

    for group in RADAR_GROUPS:
        values = scaled_radar[group].to_numpy(dtype=float)
        values = np.concatenate([values, [values[0]]])
        ax.plot(angles, values, linewidth=2.0, color=RADAR_COLORS[group], label=group)
        ax.fill(angles, values, color=RADAR_COLORS[group], alpha=0.14)

    ax.legend(loc="upper right", bbox_to_anchor=(1.24, 1.12), frameon=False, fontsize=9)
    ax.set_title(
        "GSE30427: differentiation loss vs FAP/ECM activation",
        fontsize=12.2,
        fontweight="bold",
        pad=18,
    )


def build_supp_fig6(sc_overview: pd.DataFrame, bulk_summary: pd.DataFrame, scaled_radar: pd.DataFrame) -> None:
    SUPP_FIG_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DATA_DIR.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(14.5, 10.5))
    grid = GridSpec(2, 2, height_ratios=[1.05, 1.45], width_ratios=[1.05, 1], hspace=0.24, wspace=0.18, figure=fig)

    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1], polar=True)

    draw_table_panel(
        ax_a,
        "Mouse scRNA-seq source used for Fig. 3A-H",
        sc_overview.rename(
            columns={
                "timepoint": "Time point",
                "all_cells": "All cells",
                "fibroblast_caf_cells": "Fibroblast/CAF cells",
            }
        ),
        footnote="Braf^LSL-V600E; Tpo-Cre pooled scRNA-seq. Full object: 25,095 cells; fibroblast object: 6,880 cells.",
    )
    add_panel_letter(ax_a, "A")

    panel_b = pd.DataFrame(
        [
            {
                "Cohort": "GSE55933",
                "n": "10",
                "Displayed groups": "mPTC=5\nmATC=5",
                "Primary use": "reference cohort",
            },
            {
                "Cohort": "GSE118022",
                "n": "24",
                "Displayed groups": "mNT=6\nmPTC=4\nmATC=5\nremATC=9",
                "Primary use": "Fig. 3I / Supp Fig. 6C",
            },
            {
                "Cohort": "GSE30427",
                "n": "31",
                "Displayed groups": "WT=6\nFTC=5\nATC=8\n(+ Pten_KO / p53_KO / Pten_Kras)",
                "Primary use": "Supp Fig. 6D",
            },
        ]
    )
    draw_table_panel(
        ax_b,
        "Independent mouse bulk transcriptomic cohorts",
        panel_b,
        footnote="Detailed cohort composition, source files, and analytical roles are listed in Supplementary Table 7.",
    )
    add_panel_letter(ax_b, "B")

    if not GSE118022_HEATMAP_PATH.exists():
        raise FileNotFoundError(
            f"Missing {GSE118022_HEATMAP_PATH} — GSE118022 marker heatmap produced by "
            "an upstream bulk-integration step (set EXTERNAL_DATA); see README."
        )
    heatmap_img = np.asarray(Image.open(GSE118022_HEATMAP_PATH))
    ax_c.imshow(heatmap_img)
    ax_c.axis("off")
    ax_c.set_title(
        "GSE118022: TDS versus FAP/ECM genes",
        fontsize=12.2,
        fontweight="bold",
        pad=12,
    )
    add_panel_letter(ax_c, "C")

    draw_radar(ax_d, scaled_radar)
    add_panel_letter(ax_d, "D")

    fig.suptitle(
        "Supplementary Fig. 6 | Source overview and extended transcriptomic validation for the time-resolved mouse model",
        fontsize=15,
        fontweight="bold",
        y=0.985,
    )
    fig.subplots_adjust(top=0.93, left=0.05, right=0.97, bottom=0.05)
    fig.savefig(SUPP_FIG6_PNG, dpi=320, bbox_inches="tight")
    fig.savefig(SUPP_FIG6_PDF, bbox_inches="tight")
    plt.close(fig)

    sc_overview.to_csv(SOURCE_DATA_DIR / "supp_fig6_panelA_mouse_scRNA_overview.csv", index=False)
    bulk_summary.to_csv(SOURCE_DATA_DIR / "supp_fig6_panelB_bulk_mTC_overview.csv", index=False)
    scaled_radar.to_csv(SOURCE_DATA_DIR / "supp_fig6_panelD_gse30427_radar_scaled.csv", index=False)


def main() -> None:
    full_timepoint, full_celltype, fib_subtype = fetch_scRNA_summaries()
    bulk_summary = build_bulk_summary()
    sc_overview = build_scRNA_overview(full_timepoint, fib_subtype)
    radar_raw, radar_scaled = build_radar_panel()

    write_table7(full_timepoint, full_celltype, fib_subtype, bulk_summary)
    rebuild_combined_workbook()
    build_supp_fig6(sc_overview, bulk_summary, radar_scaled)

    SOURCE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    radar_raw.to_csv(SOURCE_DATA_DIR / "supp_fig6_panelD_gse30427_radar_raw_group_means.csv", index=False)

    print(f"Supplementary Table 7 written to: {TABLE7_PATH}")
    print(f"Supplementary Fig. 6 written to: {SUPP_FIG6_PNG}")
    print(f"Combined workbook written to: {COMBINED_1_7_PATH}")


if __name__ == "__main__":
    main()
