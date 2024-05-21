# Lundager-2024
Analysis to assess enrichment of genes in each cell type. Raw count matrices of human pancreatic cells were sourced from Baron et. al 2016 [(GSE84133)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133). The data was log1p-normalized with a scale factor of 10000, and averaged across each available cell type. A non-parametric one-sided Mann-Whitney U test was employed to determine if eight genes of interest (ACSL1, CENPU, TENT5C, G6PC2, GRB10, SSTR1, SSTR2, UBE2E2) were enriched for each cell type. To account for multiple comparisons and control the false discovery rate, Benjamini-Hochberg correction was applied to the resulting p-values.
