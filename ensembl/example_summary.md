# Ensembl Compara: LDO reconstruction

158 families, 5353 genes (0 leaves removed: no aligned sequence or species not in the species tree). Internal nodes of the Ensembl trees: 4391 speciation (665 that the species tree cannot explain), 531 duplication, 267 dubious (54 of them would be valid speciations), 6 other. 83 families have no dubious node.

Truth: the Ensembl trees. Dubious nodes read as duplications or as speciations; on the families without dubious nodes both readings agree.

## Ortholog pairs, pooled (precision / recall / F1)

| Method | dubious = Dup, 158 families | dubious = Spec, 158 families | no dubious node, 83 families |
|---|---|---|---|
| LDO, no repair | 0.240 / 0.998 / 0.387 | 0.277 / 0.999 / 0.434 | 0.893 / 0.996 / 0.942 |
| LDO, count | 0.512 / 0.981 / 0.673 | 0.588 / 0.974 / 0.733 | 0.942 / 0.987 / 0.964 |
| LDO, count+sparse | 0.522 / 0.980 / 0.681 | 0.598 / 0.974 / 0.741 | 0.943 / 0.987 / 0.964 |
| LDO, count+sparse, age test (--count-age 1.0) | 0.598 / 0.972 / 0.740 | 0.684 / 0.963 / 0.800 | 0.944 / 0.986 / 0.965 |
| LDO, count+sparse+mutual | 0.603 / 0.977 / 0.746 | 0.691 / 0.969 / 0.807 | 0.943 / 0.987 / 0.964 |
| LDO, count+sparse+mutual, age test | 0.634 / 0.970 / 0.767 | 0.725 / 0.962 / 0.827 | 0.944 / 0.986 / 0.965 |
| OrthoFinder, pooled, true species tree | 0.732 / 0.923 / 0.816 | 0.837 / 0.915 / 0.874 | 0.975 / 0.977 / 0.976 |
| OrthoFinder, pooled, patched normalisation | 0.732 / 0.923 / 0.816 | 0.837 / 0.915 / 0.874 | 0.975 / 0.977 / 0.976 |

## Ortholog pairs, average per family (precision / recall)

Families without any true ortholog pair are left out; a family without predicted pair counts as precision 1.

| Method | dubious = Dup, 158 families | dubious = Spec, 158 families | no dubious node, 83 families |
|---|---|---|---|
| LDO, no repair | 0.785 / 0.974 | 0.868 / 0.974 | 0.950 / 0.951 |
| LDO, count | 0.829 / 0.965 | 0.918 / 0.964 | 0.964 / 0.948 |
| LDO, count+sparse | 0.830 / 0.965 | 0.919 / 0.964 | 0.964 / 0.948 |
| LDO, count+sparse, age test (--count-age 1.0) | 0.834 / 0.959 | 0.923 / 0.958 | 0.965 / 0.947 |
| LDO, count+sparse+mutual | 0.831 / 0.961 | 0.920 / 0.960 | 0.964 / 0.948 |
| LDO, count+sparse+mutual, age test | 0.835 / 0.958 | 0.924 / 0.957 | 0.965 / 0.947 |
| OrthoFinder, pooled, true species tree | 0.853 / 0.913 | 0.942 / 0.910 | 0.984 / 0.919 |
| OrthoFinder, pooled, patched normalisation | 0.853 / 0.913 | 0.942 / 0.910 | 0.984 / 0.919 |

## Ortholog pairs by family size, dubious = Dup (precision / recall / F1)

Gene pairs pooled within each size bin (genes per family), so that no single family dominates.

| Genes per family | families | true pairs | LDO, count+sparse | LDO, count+sparse, age test (--count-age 1.0) | LDO, count+sparse+mutual, age test | OrthoFinder, pooled, true species tree |
|---|---|---|---|---|---|---|
| 2-10 | 22 | 332 | 0.883 / 0.822 / 0.852 | 0.883 / 0.822 / 0.852 | 0.883 / 0.822 / 0.852 | 0.902 / 0.855 / 0.878 |
| 11-25 | 88 | 17804 | 0.853 / 0.998 / 0.920 | 0.853 / 0.998 / 0.920 | 0.853 / 0.998 / 0.920 | 0.861 / 0.972 / 0.913 |
| 26-50 | 31 | 11355 | 0.793 / 0.974 / 0.874 | 0.802 / 0.972 / 0.879 | 0.802 / 0.972 / 0.879 | 0.851 / 0.926 / 0.887 |
| 51-100 | 10 | 7920 | 0.730 / 0.998 / 0.844 | 0.821 / 0.975 / 0.891 | 0.821 / 0.975 / 0.891 | 0.831 / 0.859 / 0.845 |
| 101-200 | 5 | 3847 | 0.207 / 0.977 / 0.341 | 0.378 / 0.940 / 0.539 | 0.455 / 0.936 / 0.613 | 0.637 / 0.931 / 0.757 |
| > 200 | 2 | 1719 | 0.083 / 0.804 / 0.150 | 0.084 / 0.785 / 0.152 | 0.097 / 0.765 / 0.172 | 0.156 / 0.693 / 0.255 |

## Ortholog pairs by family size, dubious = Spec (precision / recall / F1)

Gene pairs pooled within each size bin (genes per family), so that no single family dominates.

| Genes per family | families | true pairs | LDO, count+sparse | LDO, count+sparse, age test (--count-age 1.0) | LDO, count+sparse+mutual, age test | OrthoFinder, pooled, true species tree |
|---|---|---|---|---|---|---|
| 2-10 | 22 | 363 | 0.984 / 0.837 / 0.905 | 0.984 / 0.837 / 0.905 | 0.984 / 0.837 / 0.905 | 1.000 / 0.868 / 0.929 |
| 11-25 | 88 | 20130 | 0.965 / 0.998 / 0.981 | 0.965 / 0.998 / 0.981 | 0.965 / 0.998 / 0.981 | 0.975 / 0.973 / 0.974 |
| 26-50 | 31 | 12441 | 0.869 / 0.974 / 0.919 | 0.879 / 0.973 / 0.924 | 0.879 / 0.973 / 0.924 | 0.928 / 0.921 / 0.925 |
| 51-100 | 10 | 8909 | 0.820 / 0.997 / 0.900 | 0.923 / 0.975 / 0.948 | 0.923 / 0.975 / 0.948 | 0.946 / 0.869 / 0.906 |
| 101-200 | 5 | 4734 | 0.255 / 0.980 / 0.405 | 0.456 / 0.922 / 0.610 | 0.548 / 0.916 / 0.686 | 0.755 / 0.896 / 0.819 |
| > 200 | 2 | 3044 | 0.138 / 0.760 / 0.234 | 0.141 / 0.742 / 0.238 | 0.163 / 0.726 / 0.266 | 0.266 / 0.668 / 0.380 |

## Clusters (precision / recall)

As cluster_analysis.py, pooled; nested duplications depend on the reading of the dubious nodes.

| Repair | dubious = Dup, 158 families | dubious = Spec, 158 families | no dubious node, 83 families |
|---|---|---|---|
| none | 0.505 / 0.495 | 0.505 / 0.474 | 0.684 / 0.674 |
| count | 0.555 / 0.569 | 0.555 / 0.547 | 0.715 / 0.708 |
| count+sparse | 0.554 / 0.570 | 0.554 / 0.548 | 0.716 / 0.709 |
| count+sparse+age | 0.552 / 0.575 | 0.552 / 0.552 | 0.716 / 0.709 |
| count+sparse+mutual | 0.551 / 0.573 | 0.551 / 0.550 | 0.716 / 0.709 |
| count+sparse+mutual+age | 0.551 / 0.576 | 0.551 / 0.553 | 0.716 / 0.709 |
