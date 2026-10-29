# simphy: LDO reconstruction vs OrthoFinder

Truth: the simphy gene trees labelled with simphy's own events (simphy_prepare.py --labels simphy). Ortholog pairs, precision / recall / F1, gene pairs pooled over the families of each setting.

25 replicates, 2500 families.

| Setting | replicates | families | genes / family | Dup nodes / family |
|---|---|---|---|---|
| dl0 | 5 | 500 | 25.0 | 0.00 |
| dl1 | 5 | 500 | 28.2 | 4.24 |
| dl2 | 5 | 500 | 44.9 | 14.80 |
| dl3 | 5 | 500 | 49.4 | 20.31 |
| dl5 | 5 | 500 | 71.2 | 37.88 |

## Ortholog pairs, all families

Each method on all the families it was run on (per-family OrthoFinder may cover a subset only).

| Method | dl0 | dl1 | dl2 | dl3 | dl5 | All |
|---|---|---|---|---|---|---|
| LDO, no repair | 1.000 / 1.000 / 1.000 | 0.545 / 1.000 / 0.705 | 0.253 / 0.999 / 0.403 | 0.166 / 0.999 / 0.284 | 0.129 / 0.999 / 0.229 | 0.199 / 0.999 / 0.331 |
| LDO, count | 1.000 / 1.000 / 1.000 | 0.581 / 0.944 / 0.719 | 0.272 / 0.959 / 0.424 | 0.169 / 0.965 / 0.288 | 0.134 / 0.967 / 0.235 | 0.207 / 0.966 / 0.340 |
| LDO, count+sparse | 1.000 / 1.000 / 1.000 | 0.692 / 0.938 / 0.796 | 0.346 / 0.945 / 0.507 | 0.221 / 0.954 / 0.358 | 0.195 / 0.952 / 0.323 | 0.280 / 0.955 / 0.434 |
| LDO, count+sparse, age test (--count-age 1.0) | 1.000 / 1.000 / 1.000 | 0.799 / 0.861 / 0.829 | 0.475 / 0.799 / 0.595 | 0.440 / 0.818 / 0.572 | 0.298 / 0.751 / 0.427 | 0.447 / 0.820 / 0.578 |
| OrthoFinder, pooled per replicate | 1.000 / 0.651 / 0.789 | 0.873 / 0.671 / 0.759 | 0.680 / 0.646 / 0.663 | 0.503 / 0.617 / 0.554 | 0.388 / 0.562 / 0.459 | 0.552 / 0.616 / 0.582 |
| OrthoFinder, per family | 1.000 / 0.000 / 0.000 | 0.735 / 0.423 / 0.537 | 0.556 / 0.510 / 0.532 | 0.322 / 0.704 / 0.442 | 0.254 / 0.566 / 0.351 | 0.346 / 0.499 / 0.409 |
| OrthoFinder, per family, patched normalisation | 1.000 / 0.643 / 0.782 | 0.830 / 0.743 / 0.784 | 0.770 / 0.483 / 0.594 | 0.349 / 0.721 / 0.470 | 0.246 / 0.525 / 0.335 | 0.421 / 0.612 / 0.499 |

## Ortholog pairs, families with per-family OrthoFinder runs

All methods restricted to the 248 families that have per-family OrthoFinder results.

| Method | dl0 | dl1 | dl2 | dl3 | dl5 | All |
|---|---|---|---|---|---|---|
| LDO, no repair | 1.000 / 1.000 / 1.000 | 0.361 / 0.999 / 0.531 | 0.272 / 1.000 / 0.428 | 0.053 / 0.999 / 0.101 | 0.083 / 0.999 / 0.153 | 0.106 / 0.999 / 0.192 |
| LDO, count | 1.000 / 1.000 / 1.000 | 0.362 / 0.964 / 0.526 | 0.280 / 0.980 / 0.436 | 0.052 / 0.966 / 0.099 | 0.083 / 0.873 / 0.152 | 0.107 / 0.947 / 0.192 |
| LDO, count+sparse | 1.000 / 1.000 / 1.000 | 0.608 / 0.954 / 0.743 | 0.361 / 0.973 / 0.527 | 0.061 / 0.953 / 0.114 | 0.122 / 0.839 / 0.213 | 0.136 / 0.931 / 0.238 |
| LDO, count+sparse, age test (--count-age 1.0) | 1.000 / 1.000 / 1.000 | 0.712 / 0.910 / 0.799 | 0.392 / 0.864 / 0.539 | 0.411 / 0.626 / 0.496 | 0.156 / 0.600 / 0.247 | 0.341 / 0.753 / 0.470 |
| OrthoFinder, pooled per replicate | 1.000 / 0.643 / 0.782 | 0.813 / 0.725 / 0.767 | 0.660 / 0.606 / 0.632 | 0.420 / 0.802 / 0.552 | 0.280 / 0.635 / 0.388 | 0.454 / 0.686 / 0.546 |
| OrthoFinder, per family | 1.000 / 0.000 / 0.000 | 0.735 / 0.423 / 0.537 | 0.556 / 0.510 / 0.532 | 0.322 / 0.704 / 0.442 | 0.254 / 0.566 / 0.351 | 0.346 / 0.499 / 0.409 |
| OrthoFinder, per family, patched normalisation | 1.000 / 0.643 / 0.782 | 0.830 / 0.743 / 0.784 | 0.770 / 0.483 / 0.594 | 0.349 / 0.721 / 0.470 | 0.246 / 0.525 / 0.335 | 0.421 / 0.612 / 0.499 |

## Average per family (precision / recall)

Each family scored on its own, then averaged; families without any ortholog pair (all genes in one species) are left out for every method. A family with no predicted pair counts as precision 1.

| Method | families | dl0 | dl1 | dl2 | dl3 | dl5 | All |
|---|---|---|---|---|---|---|---|
| LDO, no repair | 2479 | 1.000 / 1.000 | 0.825 / 1.000 | 0.684 / 0.998 | 0.614 / 0.998 | 0.498 / 0.997 | 0.725 / 0.999 |
| LDO, count | 2479 | 1.000 / 1.000 | 0.852 / 0.982 | 0.713 / 0.977 | 0.650 / 0.969 | 0.528 / 0.968 | 0.750 / 0.979 |
| LDO, count+sparse | 2479 | 1.000 / 1.000 | 0.875 / 0.978 | 0.753 / 0.967 | 0.681 / 0.960 | 0.566 / 0.955 | 0.776 / 0.972 |
| LDO, count+sparse, age test (--count-age 1.0) | 2479 | 1.000 / 1.000 | 0.892 / 0.942 | 0.789 / 0.903 | 0.713 / 0.893 | 0.599 / 0.851 | 0.800 / 0.918 |
| OrthoFinder, pooled per replicate | 2479 | 1.000 / 0.651 | 0.934 / 0.734 | 0.858 / 0.719 | 0.795 / 0.702 | 0.709 / 0.703 | 0.860 / 0.702 |
| OrthoFinder, per family | 248 | 1.000 / 0.000 | 0.949 / 0.350 | 0.858 / 0.388 | 0.803 / 0.487 | 0.633 / 0.573 | 0.849 / 0.359 |
| OrthoFinder, per family, patched normalisation | 248 | 1.000 / 0.643 | 0.950 / 0.806 | 0.883 / 0.707 | 0.829 / 0.632 | 0.662 / 0.609 | 0.865 / 0.679 |

## Average per family, families with per-family OrthoFinder runs

As above, restricted to the families with per-family OrthoFinder results.

| Method | families | dl0 | dl1 | dl2 | dl3 | dl5 | All |
|---|---|---|---|---|---|---|---|
| LDO, no repair | 248 | 1.000 / 1.000 | 0.863 / 1.000 | 0.723 / 1.000 | 0.607 / 0.997 | 0.454 / 0.999 | 0.729 / 0.999 |
| LDO, count | 248 | 1.000 / 1.000 | 0.875 / 0.983 | 0.737 / 0.990 | 0.648 / 0.953 | 0.483 / 0.931 | 0.749 / 0.971 |
| LDO, count+sparse | 248 | 1.000 / 1.000 | 0.896 / 0.980 | 0.778 / 0.987 | 0.686 / 0.944 | 0.518 / 0.914 | 0.776 / 0.965 |
| LDO, count+sparse, age test (--count-age 1.0) | 248 | 1.000 / 1.000 | 0.904 / 0.960 | 0.818 / 0.939 | 0.724 / 0.869 | 0.551 / 0.801 | 0.799 / 0.914 |
| OrthoFinder, pooled per replicate | 248 | 1.000 / 0.643 | 0.948 / 0.797 | 0.878 / 0.729 | 0.811 / 0.672 | 0.660 / 0.710 | 0.859 / 0.710 |
| OrthoFinder, per family | 248 | 1.000 / 0.000 | 0.949 / 0.350 | 0.858 / 0.388 | 0.803 / 0.487 | 0.633 / 0.573 | 0.849 / 0.359 |
| OrthoFinder, per family, patched normalisation | 248 | 1.000 / 0.643 | 0.950 / 0.806 | 0.883 / 0.707 | 0.829 / 0.632 | 0.662 / 0.609 | 0.865 / 0.679 |

## Clusters (LDO trees)

Cluster precision (exact + polytomy) / recall (reference clusters other than nested duplications), as cluster_analysis.py, pooled.

| Repair | dl0 | dl1 | dl2 | dl3 | dl5 | All |
|---|---|---|---|---|---|---|
| none | 1.000 / 1.000 | 0.772 / 0.746 | 0.645 / 0.591 | 0.571 / 0.518 | 0.503 / 0.432 | 0.664 / 0.608 |
| count | 1.000 / 1.000 | 0.803 / 0.787 | 0.672 / 0.636 | 0.599 / 0.564 | 0.521 / 0.475 | 0.682 / 0.646 |
| count+sparse | 1.000 / 1.000 | 0.824 / 0.817 | 0.685 / 0.673 | 0.615 / 0.604 | 0.526 / 0.516 | 0.689 / 0.679 |
| count+sparse+age | 1.000 / 1.000 | 0.765 / 0.760 | 0.600 / 0.604 | 0.537 / 0.548 | 0.423 / 0.446 | 0.612 / 0.623 |
