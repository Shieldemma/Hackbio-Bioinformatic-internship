Contributors: Ajisegiri Sewedo B, Balqees Mansour, Salaam Ridwan, Adekoya Adepeju, Melvin Khakabo, Elizabeth Abodunrin, Otuekong Emmanuel, Samuel Aladegbaiye

## Identifying and Validating Biomarkers for Breast Invasive Carcinoma through Differential Gene Expression and Functional Enrichment Techniques
## Introduction
Breast Invasive Carcinoma (BRCA) is a prevalent and aggressive form of breast cancer and the second most common cancer among women globally (De Schepper et al., 2022). Molecular techniques have identified biomarkers that can predict tumor behavior and improve treatment responses in various cancers (Neves Rebello et al., 2023). Similarly, this study seeks to use a combination of differential gene expression and functional enrichment analysis to identify important biomarkers for breast cancer.
## Methodology
In this study, the BRCA transcriptome dataset was retrieved from the TCGA database using the R studio. The dataset, which includes the gene sets associated with multiple samples, was prepared using the GDCprepare function. Validation checks were conducted and a DESeqDataSet object was created from the count matrix and sample information.  A model design based on tissue type (tumor vs. normal tissue) was specified. Genes with low counts across samples were filtered out to improve statistical power. The reference level for tissue type was set to "normal" to facilitate comparisons against this baseline. The DESeq function and eKegg function were executed to perform differential expression analysis and functional enrichment analysis, respectively. 
## Results
A volcano plot was used to visualize the upregulated and downregulated genes together. Additionally, boxplot charts were used to show the top 20 upregulated genes and top 20 downregulated genes. Similarly, the 20 functionally enriched pathways were shown in a lollipop chart. 
## Discussion
This study showed that Neuroactive ligand-receptor interaction is the most functionally enriched pathway. It modulates the tumor microenvironment, influencing cancer cell behavior and immune response (Yang et al., 2023). This is followed by PI3K Akt signaling pathway, which drives cell proliferation, survival, and metastasis (He et al., 2021). Cytokine-cytokine receptor interaction regulates immune responses; can promote tumor growth and metastasis through inflammation (Esquivel-Velázquez et al., 2015). Meanwhile, the PPAR signaling pathway is involved in lipid metabolism and inflammation (Michalik and Wahli, 2008). The fifth pathway on the list is the cell cycle, which controls cell division. The dysregulation leads to uncontrolled proliferation and tumor development (Williams and Stoeber, 2012). 
## Conclusion
This study identified several upregulated and downregulated genes as potential biomarkers for Breast Invasive Carcinoma (BIC), highlighting critical pathways such as the PI3K-Akt signaling pathway. Future research should focus on validating these biomarkers in clinical settings and exploring their roles in tumor progression. 

## References
De Schepper, M., Vincent-Salomon, A., Christgen, M., Van Baelen, K., Richard, F., Tsuda, H., Kurozumi, S., Jose Brito, M., Cserni, G., Schnitt, S., Larsimont, D., Kulka, J., Luis Fernandez, P., Rodríguez-Martínez, P., Aula Olivar, A., Melendez, C., Van Bockstal, M., Kovacs, A., Varga, Z., & Wesseling, J. (2022). Results of a worldwide survey on the currently used histopathological diagnostic criteria for invasive lobular breast cancer. Modern Pathology, 35(12), 1812–1820. https://doi.org/10.1038/s41379-022-01135-2

Esquivel-Velázquez, M., Ostoa-Saloma, P., Palacios-Arreola, M. I., Nava-Castro, K. E., Castro, J. I., and Morales-Montor, J. (2015). The role of cytokines in breast cancer development and progression. Journal of Interferon & Cytokine Research, 35(1), 1-16.https://doi.org/10.1089/jir.2014.0026. 

He, Y., Sun, M. M., Zhang, G. G., Yang, J., Chen, K. S., Xu, W. W., and Li, B. (2021). Targeting PI3K/Akt signal transduction for cancer therapy. Signal transduction and targeted therapy, 6(1), 425. https://doi.org/10.1038/s41392-021-00828-5.

Michalik, L., and Wahli, W. (2008). PPARs mediate lipid signaling in inflammation and cancer. PPAR research, 2008(1), 134059.https://doi.org/10.1155/2008/134059. 

Neves Rebello Alves, L., Dummer Meira, D., Poppe Merigueti, L., Correia Casotti, M., do Prado Ventorim, D., Ferreira Figueiredo Almeida, J., Pereira de Sousa, V., Cindra Sant'Ana, M., Gonçalves Coutinho da Cruz, R., Santos Louro, L., Mendonça Santana, G., Erik Santos Louro, T., Evangelista Salazar, R., Ribeiro Campos da Silva, D., Stefani Siqueira Zetum, A., Silva Dos Reis Trabach, R., Imbroisi Valle Errera, F., de Paula, F., de Vargas Wolfgramm Dos Santos, E., Fagundes de Carvalho, E., … Drumond Louro, I. (2023). Biomarkers in Breast Cancer: An Old Story with a New End. Genes, 14(7), 1364. https://doi.org/10.3390/genes14071364.

Williams, G. H., and Stoeber, K. (2012). The cell cycle and cancer. The Journal of pathology, 226(2), 352-364. https://doi.org/10.1073/pnas.94.7.2776. 

Yang, Y., Li, J., Jing, C., Zhai, Y., Bai, Z., Yang, Y., and Deng, W. (2023). Inhibition of neuroactive ligand–receptor interaction pathway can enhance immunotherapy response in colon cancer: an in-silico study. Expert Review of Anticancer Therapy, 23(11), 1205-1215.https://doi.org/10.1080/14737140.2023.2245567.


