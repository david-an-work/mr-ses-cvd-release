"""make_figures.py

Draws the two manuscript figures from the saved results tables (no estimation).
  Figure 1: results/figures/mr_primary_grid_forest.{pdf,png}  <- results/tables/mr_primary_grid.tsv
  Figure 2: results/figures/mvmr_appraisal_cad.{pdf,png}      <- results/tables/mvmr_appraisal.tsv (CAD, Aragam EUR)
Usage (repository root): python3 manuscript/make_figures.py
Python 3.11, matplotlib 3.x, pandas.
"""
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import pandas as pd

OUT = "results/figures"
SIZES = (8, 7, 6)  # base, secondary, tick
COL = {"EA": "#1b6ca8", "Income": "#d17a22", "Income_Hill": "#d17a22", "TDI": "#6a6a6a",
       "TDI_EUR": "#6a6a6a", "Occupation": "#3a8f5c", "Occupation_ISEI": "#3a8f5c"}


def style():
    base, sec, tick = SIZES
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.size": base, "axes.labelsize": base, "axes.titlesize": base,
        "legend.fontsize": sec, "xtick.labelsize": tick, "ytick.labelsize": tick, "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "savefig.dpi": 300, "savefig.bbox": "tight", "axes.titlelocation": "left",
        "pdf.fonttype": 42, "ps.fonttype": 42})


def log_ticks(ax, ticks):
    ax.set_xscale("log")
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())


def figure1():
    P = pd.read_csv("results/tables/mr_primary_grid.tsv", sep="\t")
    sets = ["EA4", "Kweon", "TDI", "CAMSIS"]
    ecol = {"EA4": COL["EA"], "Kweon": COL["Income"], "TDI": COL["TDI"], "CAMSIS": COL["Occupation"]}
    labs = {"EA4": "Education (EA4)", "Kweon": "Income (Kweon)", "TDI": "Low deprivation (TDI)",
            "CAMSIS": "Occupation (CAMSIS)"}
    outl = {"cad_aragam_eur": "CAD\n(Aragam EUR)", "pad_mvp_eur": "PAD\n(MVP EUR)",
            "hf_hermes2": "Heart failure\n(HERMES 2.0)", "stroke_gigastroke_multi": "Any stroke\n(GIGASTROKE multi-ancestry)"}
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.9), sharey=True, gridspec_kw={"wspace": 0.08})
    for ax, o in zip(axes, outl):
        for i, s in enumerate(sets):
            r = P[(P.set == s) & (P.outcome == o)].iloc[0]
            y = len(sets) - 1 - i
            ax.plot([r.or_lo_bonf, r.or_hi_bonf], [y, y], color="0.45", lw=0.9, zorder=1)
            ax.plot([r.or_lo, r.or_hi], [y, y], color=ecol[s], lw=2.4, zorder=2)
            ax.plot(r["or"], y, "o", color=ecol[s], ms=4.5, mew=1.2, zorder=3,
                    mfc=ecol[s] if r.bonferroni_sig else "white")
        ax.axvline(1, color="0.5", lw=0.7, ls=":")
        log_ticks(ax, [0.2, 0.5, 0.7, 1.0])
        ax.set_xlim(0.14, 1.45)
        ax.set_title(outl[o], loc="left")
        ax.set_ylim(-0.6, 3.6)
    axes[0].set_yticks(range(4))
    axes[0].set_yticklabels([labs[s] for s in sets[::-1]])
    fig.supxlabel("Odds ratio per SD higher genetically predicted SES (IVW, multiplicative random effects; log scale)",
                  y=0.0, fontsize=SIZES[0])
    key = [Line2D([], [], marker="o", color="0.3", mfc="0.3", ls="", ms=4.5, label="p < 0.003125 (Bonferroni)"),
           Line2D([], [], marker="o", color="0.3", mfc="white", ls="", ms=4.5, label="not significant at 0.003125"),
           Line2D([], [], color="0.3", lw=2.4, label="95% CI"),
           Line2D([], [], color="0.45", lw=0.9, label="99.69% CI (Bonferroni level)")]
    fig.legend(handles=key, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.13), handlelength=1.6,
               columnspacing=1.2)
    fig.savefig(f"{OUT}/mr_primary_grid_forest.png")
    fig.savefig(f"{OUT}/mr_primary_grid_forest.pdf")
    return fig


def figure2():
    R = pd.read_csv("results/tables/mvmr_appraisal.tsv", sep="\t")
    ms = ["EA+Income+TDI+Occupation", "EA+Income_Hill+TDI+Occupation", "EA+Income", "EA+TDI", "EA+Occupation",
          "LiD_S1", "EA+Occupation_ISEI"]
    mlab = {"EA+Income+TDI+Occupation": "A1 four-measure (Kweon income)",
            "EA+Income_Hill+TDI+Occupation": "A5 four-measure (Hill income)",
            "EA+Income": "A2 EA4 + Kweon income", "EA+TDI": "A3 EA4 + deprivation (TDI)",
            "EA+Occupation": "A4 EA4 + occupation (CAMSIS)",
            "LiD_S1": "B1/B2 EA4 + Hill income", "EA+Occupation_ISEI": "B3 EA4 + occupation (ISEI)"}
    elab = {"EA": "Education", "Income": "Income", "Income_Hill": "Income (Hill)\u2020", "TDI": "Low deprivation",
            "Occupation": "Occupation", "Occupation_ISEI": "Occupation (ISEI)"}
    c = R[(R.outcome == "cad_aragam_eur") & (R.model.isin(ms))]
    rows, yt, yl, y = [], [], [], 0
    for m in ms:
        for _, r in c[c.model == m].iterrows():
            rows.append((y, r)); yt.append(y); yl.append(f"{mlab[m]} \u00b7 {elab[r.exposure]}"); y -= 1
        y -= 0.6
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 5.0), sharey=True,
                                 gridspec_kw={"wspace": 0.06, "width_ratios": [1, 1.15]})
    for y, r in rows:
        col = COL[r.exposure]
        a1.plot([r.condF_min, r.condF_max], [y, y], color=col, lw=3, solid_capstyle="butt")
        a1.plot([r.condF_x1], [y], "|", color="k", ms=7)
        a2.plot([r.or_lo, r.or_hi], [y, y], color=col, lw=2)
        a2.plot(r["or"], y, "o", color=col, mfc="white", ms=4.5, mew=1.2)
    a1.axvline(10, color="#b2182b", lw=0.9, ls="--")
    a1.text(9.7, rows[0][0] + 1.0, "F = 10 (study reporting requirement)", ha="right", va="bottom",
            fontsize=SIZES[2], color="#b2182b")
    a1.set_xlim(0, 11)
    a1.set_xlabel("Conditional F\n(bar: 0\u20131.5 \u00d7 LDSC correlation; tick: 1\u00d7)")
    a1.set_title("Every coefficient is below F = 10", loc="left")
    a2.axvline(1, color="0.5", lw=0.7, ls=":")
    log_ticks(a2, [0.6, 0.8, 1.0, 1.2, 1.6])
    a2.set_xlabel("IVW-MVMR odds ratio (95% CI; log scale)\nper SD; \u2020Hill: per released unit")
    a2.set_title("Weak-instrument-limited: not effect sizes", loc="left")
    a1.set_yticks(yt); a1.set_yticklabels(yl)
    a1.set_ylim(rows[-1][0] - 0.8, rows[0][0] + 2.2)
    for ax in (a1, a2):
        ax.margins(x=0.04)
    fig.savefig(f"{OUT}/mvmr_appraisal_cad.png")
    fig.savefig(f"{OUT}/mvmr_appraisal_cad.pdf")
    return fig


def overlaps(fig):
    r = fig.canvas.get_renderer()
    t = [(x, x.get_window_extent(r)) for x in fig.findobj(mpl.text.Text) if x.get_text().strip() and x.get_visible()]
    return [(a.get_text()[:30], b.get_text()[:30]) for i, (a, ba) in enumerate(t) for b, bb in t[i + 1:] if ba.overlaps(bb)]


if __name__ == "__main__":
    style()
    for f in (figure1(), figure2()):
        print("text overlaps:", overlaps(f))
