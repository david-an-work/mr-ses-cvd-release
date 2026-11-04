"""verify_supplementary_tables.py

Reconciles every sheet of manuscript/supplementary_tables.xlsx with its source file (named on the Contents sheet):
same columns, same row count, and cell-by-cell equality (numbers to a relative tolerance of 1e-12, strings exactly,
missing values matched, booleans compared case-insensitively). Writes manuscript/verification/supplementary_tables_reconciliation.tsv.
Usage (repository root): python3 manuscript/verify_supplementary_tables.py
"""
import numpy as np
import pandas as pd

X = "manuscript/supplementary_tables.xlsx"
contents = pd.read_excel(X, sheet_name="Contents")
out = []
for _, r in contents.iterrows():
    sheet, src = r.iloc[0], r.iloc[2]
    xl = pd.read_excel(X, sheet_name=sheet[:31], dtype=object)
    sv = pd.read_csv(src, sep="\t", dtype=object)
    res = {"sheet": sheet, "source": src, "rows_sheet": len(xl), "rows_source": len(sv),
           "cols_equal": list(xl.columns) == list(sv.columns), "n_cell_mismatch": None, "first_mismatch": ""}
    if res["cols_equal"] and len(xl) == len(sv):
        bad = []
        for c in sv.columns:
            a, b = xl[c], sv[c]
            an, bn = pd.to_numeric(a, errors="coerce"), pd.to_numeric(b, errors="coerce")
            both_num = an.notna() & bn.notna()
            num_ok = np.isclose(an[both_num].astype(float), bn[both_num].astype(float), rtol=1e-12, atol=0)
            miss = a.isna() | (a.astype(str) == "")
            missb = b.isna() | (b.astype(str) == "")
            str_mask = ~both_num & ~(miss & missb)
            norm = lambda v: v.astype(str).str.strip().str.lower()  # Excel stores TRUE/FALSE as booleans
            str_ok = norm(a[str_mask]).values == norm(b[str_mask]).values
            nb = int((~num_ok).sum() + (~str_ok).sum())
            if nb:
                idx = list(a[both_num].index[~num_ok]) + list(a[str_mask].index[~str_ok])
                bad.append((c, nb, idx[0], a[idx[0]], b[idx[0]]))
        res["n_cell_mismatch"] = sum(x[1] for x in bad)
        res["first_mismatch"] = "; ".join(f"{c} row {i}: sheet={v1!r} source={v2!r} (n={n})" for c, n, i, v1, v2 in bad[:3])
    out.append(res)
R = pd.DataFrame(out)
import os
os.makedirs("manuscript/verification", exist_ok=True)
R.to_csv("manuscript/verification/supplementary_tables_reconciliation.tsv", sep="\t", index=False)
print(R[["sheet", "rows_sheet", "rows_source", "cols_equal", "n_cell_mismatch"]].to_string(index=False))
print(R.loc[R.n_cell_mismatch.fillna(1) > 0, "first_mismatch"].to_string())
