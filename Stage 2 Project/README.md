The vault of knowledge is inexhaustible, this is the report of my progress at the on going @hackbio  internship, specially for stage 2!

**Overview**

This project focuses on analysing a gene expression dataset pertaining to glioblastoma, aiming to visualize differential gene expression patterns and perform functional enrichment analysis. The insights garnered from this analysis are critical for understanding the underlying biological mechanisms of glioblastoma and potential therapeutic targets.

**Sample Information**

A sample information data frame was created to classify the tissue types, enabling a clear distinction between "Solid Tissue Normal," "Primary Tumor," and other sample types.

𝗛𝗲𝗮𝘁𝗺𝗮𝗽 𝗩𝗶𝘀𝘂𝗮𝗹𝗶𝘇𝗮𝘁𝗶𝗼𝗻: Generated sequential and divergent color heatmaps using the 𝗵𝗲𝗮𝘁𝗺𝗮𝗽.𝟮 package from 𝗴𝗽𝗹𝗼𝘁𝘀 to identify gene expression patterns for heatmap by Genes, samples and both

**Clustering Approaches**

- **Gene Clustering**: Revealed groups of co-expressed genes.
- **Sample Clustering**: Showed similarities among samples based on expression profiles.
- **Combined Clustering**: Offered insights into how genes and samples interact.

𝗗𝗶𝗳𝗳𝗲𝗿𝗲𝗻𝘁𝗶𝗮𝗹 𝗘𝘅𝗽𝗿𝗲𝘀𝘀𝗶𝗼𝗻 𝗔𝗻𝗮𝗹𝘆𝘀𝗶𝘀: We conducted differential expression analysis using the DESeq2 package, setting thresholds of adjusted p-value < 0.05 and |log2FoldChange| > 2 to identify significant genes. This process allowed us to isolate up regulated and down regulated genes.

**Significant Findings**

- **Upregulated Genes**: A set of genes with increased expression levels in tumor samples was identified.
- **Downregulated Genes**: Similarly, genes with decreased expression levels were also noted.

**Functional Enrichment Analysis**

We utilized the `clusterProfiler` package to perform functional enrichment analysis on the identified DEGs. The enrichment results were visualized using a lollipop plot to depict the top 5 KEGG pathways, where the size of each point indicates the number of genes associated with the pathway, and color intensity reflects the significance of each pathway (based on p-value).

**Conclusion**

This analysis demonstrates the intricate patterns of gene expression in glioblastoma and identifies potential pathways for further investigation. The visualizations and enrichment results provide a solid foundation for understanding the biological processes involved in glioblastoma progression.
