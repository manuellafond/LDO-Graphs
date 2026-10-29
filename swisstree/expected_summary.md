# Results

Ortholog pairs, precision / recall / F1, gene pairs pooled over the families (micro average).
Truth: the curated SwissTree trees, 8 families (ST006 is scored separately below: its tree file has no event tags).
Clusters: cluster precision / recall of the reconstructed trees against the curated trees, all 9 families.

| Run | Curated, 8 families | Clusters P / R |
|---|---|---|
| LDO, no repair | 0.843 / 0.999 / 0.914 | 0.726 / 0.684 |
| LDO, count | 0.898 / 0.997 / 0.945 | 0.805 / 0.780 |
| LDO, count+sparse | 0.947 / 0.987 / 0.967 | 0.831 / 0.818 |
| LDO, count+sparse + age test (--count-age 1.0) | 0.970 / 0.982 / 0.976 | 0.820 / 0.799 |
| OrthoFinder, pooled, true species tree | 0.991 / 0.819 / 0.897 | - |
| OrthoFinder, pooled, patched normalisation | 0.996 / 0.962 / 0.979 | - |
| OrthoFinder, per family, true species tree | 0.990 / 0.753 / 0.855 | - |
| OrthoFinder, per family, patched normalisation | 0.996 / 0.962 / 0.979 | - |

## ST006, curators' pair list

ST006 against the ortholog pairs of the curators' gene_relationships.txt. st: plain conversion; stf: leaves also matched by gene name and fragments removed (OrthoFinder was run on st only).

| Run | st | stf |
|---|---|---|
| LDO, no repair | 0.225 / 1.000 / 0.367 | 0.382 / 1.000 / 0.552 |
| LDO, count | 0.501 / 0.962 / 0.659 | 0.594 / 1.000 / 0.745 |
| LDO, count+sparse | 0.623 / 0.962 / 0.756 | 0.684 / 0.960 / 0.799 |
| LDO, count+sparse + age test (--count-age 1.0) | 0.655 / 0.928 / 0.768 | 0.719 / 0.951 / 0.819 |
| OrthoFinder, pooled, true species tree | 0.841 / 0.959 / 0.896 | - |
| OrthoFinder, pooled, patched normalisation | 0.864 / 0.935 / 0.898 | - |
| OrthoFinder, per family, true species tree | 0.851 / 0.877 / 0.864 | - |
| OrthoFinder, per family, patched normalisation | 0.841 / 0.959 / 0.896 | - |

## Per tree, curated truth

Precision / recall.

| Tree | LDO none | LDO count | LDO count+sparse | LDO count+sparse+age | OF unpatched | OF patched |
|---|---|---|---|---|---|---|
| ST001 | 0.626 / 1.000 | 0.798 / 1.000 | 1.000 / 1.000 | 0.968 / 0.974 | 0.927 / 1.000 | 1.000 / 1.000 |
| ST002 | 0.452 / 1.000 | 0.546 / 1.000 | 0.550 / 1.000 | 0.728 / 0.959 | 1.000 / 1.000 | 1.000 / 1.000 |
| ST003 | 0.976 / 0.993 | 1.000 / 0.993 | 1.000 / 0.993 | 1.000 / 0.993 | 1.000 / 0.906 | 1.000 / 0.997 |
| ST004 | 0.927 / 1.000 | 0.927 / 1.000 | 0.992 / 0.983 | 0.986 / 0.980 | 0.993 / 0.793 | 0.995 / 1.000 |
| ST005 | 0.555 / 1.000 | 0.613 / 1.000 | 0.701 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 |
| ST008 | 0.994 / 1.000 | 0.994 / 1.000 | 0.994 / 1.000 | 0.994 / 1.000 | 0.994 / 1.000 | 0.994 / 1.000 |
| ST009 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 0.413 | 1.000 / 0.607 |
| ST010 | 0.469 / 1.000 | 0.918 / 0.918 | 0.918 / 0.918 | 0.918 / 0.918 | 1.000 / 1.000 | 1.000 / 1.000 |
| **All** | 0.843 / 0.999 | 0.898 / 0.997 | 0.947 / 0.987 | 0.970 / 0.982 | 0.991 / 0.819 | 0.996 / 0.962 |
