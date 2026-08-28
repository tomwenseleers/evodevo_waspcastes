Final GO and KEGG annotation

Direct wasp GO combines Galaxy InterProScan, Galaxy EggNOG-mapper and
EXCON EggNOG-mapper by set union at gene level. For each N13 HOG, direct
terms are retained through the established Polistes-Vespula conserved-GO
rule and combined with experimentally supported GO terms transferred from
the mapped Drosophila melanogaster orthologue.

The final GO table is N13_HOG_GO_final_long.tsv.gz. It uses EXP, IDA, IPI,
IMP, IGI, IEP and their high-throughput equivalents from FlyBase FB2025_03,
excluding NOT and ND annotations. Alternative GO sensitivity tables are
not exported by this final script. Source provenance remains auditable.

For KEGG pathway enrichment with clusterProfiler::enricher(), use one of
the KEGG_pathway_*_TERM2GENE.tsv files plus KEGG_pathway_TERM2NAME.tsv.
The union table maximises coverage; the conserved table requires the same
pathway annotation in both wasp species represented by an N13 HOG.
InterProScan contributes GO annotations only in this build; its separate
MetaCyc/Reactome pathway field is not treated as KEGG annotation.
InterPro rows with 13 rather than 15 fields lack optional GO/pathway fields;
the parser audit found no shifted or overfilled records.
Reload all key tables with load('final_GO_KEGG_analysis_workspace.RData').
