# ZoomDeflate: Imputing scRNAseq Dropouts Through Matrix Reconstruction

**Authors**: [Jeremy P D'Silva](https://github.com/jpdsilva), [Jaeyoon Kim](https://github.com/jkim5209), [Nikhil Shankar](https://github.com/nshankar)

Single-cell RNA sequencing (scRNAseq) experiments are used to measure the gene expression of many cells at single-cell resolution. Plate-based scRNAseq methods are susceptible to a type of error called zero-inflation or dropouts, wherein transcripts that are present in the cell are not detected by scRNAseq, leading to artificial zeros (dropouts) in the expression data.

This library tests three methods to reconstruct zero-inflated single-cell data: ClusterMean (a novel algorithm), ZoomDeflate (independent rediscovery of [Mongia et al](https://github.com/aanchalMongia/McImpute_scRNAseq)), and ALRA (the low-rank approximation method developed in [Linderman et al](https://github.com/KlugerLab/ALRA)). The accompanying paper **`ZoomDeflate.pdf`** details our methodology and results.