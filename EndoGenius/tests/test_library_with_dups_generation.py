# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 10:52:26 2025

@author: lafields2
"""

import pandas as pd
from library_with_dups_generation import lib_w_dubs_generate

def test_lib_with_dups_generation(tmp_path):
    src = tmp_path / "src"
    dst = tmp_path / "dst"
    src.mkdir()
    dst.mkdir()

    # Create two small CSVs
    df1 = pd.DataFrame({"Sequence": ["AAA", "BBB"], "hyperscore": [10, 20]})
    df2 = pd.DataFrame({"Sequence": ["AAA", "CCC"], "hyperscore": [30, 40]})
    (src / "lib1.csv").write_text(df1.to_csv(index=False))
    (src / "lib2.csv").write_text(df2.to_csv(index=False))

    out_path = lib_w_dubs_generate(str(src), str(dst))
    assert (tmp_path / "dst" / "library_with_duplicates.csv").exists()
    assert out_path.endswith("library_with_duplicates.csv")

    merged = pd.read_csv(out_path)
    # Expect 4 rows = 2 + 2, with origin labels present
    assert len(merged) == 4
    assert "Origin library" in merged.columns
    assert set(merged["Origin library"]) == {"lib1.csv", "lib2.csv"}