# Challenge
*Coffea arabica* challenge for RUST stress

**_Coffea arabica_ challenge adressed by RNAseq under different types of mutagenesis at different levels of RUST infection**

## CODE

Using DESeq2 normalization

Unsupervised plots: (PCA, clustering/heatmaps) → They don’t use the DESeq2 design. They depend only on the expression matrix we feed. Los análisis no supervisados (PCA, clustering/heatmaps) no dependen del design de DESeq2; operan sobre la matriz de expresión preprocesada (VST / z-score) 

### Packages

      library(kohonen);	# This is the library for the SOM (Self-Organizing maps)
      library(ggplot2);	# This library is for transparency in the colors
      library(gplots);	# Easy heatmaps también pheat más facile
      library(VennDiagram);	# To creat Venn diagrams
      library(pheatmap);	# Pretty heatmaps
      library(dendsort);	# Sorting dendrograms, (ayudar a rotar los dendogramas, to make the color differences match better)
      library(DESeq2);	# Normalization and everything related to that

### INPUT

      url_rna_counts <- "https://raw.githubusercontent.com/jmvillalobos/RNAseq_curso2025/refs/heads/main/Cara_RNAc90_counts.txt"
      # url_rna_counts <- "Cara_RNAc90_counts.txt" # In case something happened with the link

      
      # MAIN RAW DATA reading
      rust_counts <- read.table(url_rna_counts, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

      
      head(rust_counts)

### Sample short names

      # file into rust_counts
      counts <- rust_counts # It just changed the name of the varibale, can be cleaned

      
      ## 1) Build simple names (strip(quitar) everything after the sample tag)
      raw_names <- colnames(counts) 
      simple    <- sub("^([TH])(\\d+).*", "\\1\\2", raw_names)   # E.g. "H10_trim_Cara_c90" -> "H10"
            
            # sub() → Functions as sub(pattern, replacement, x)
            # PATTERN
            # ^ → anchor the start of the string
            # ([TH]) → Encapsulate "group 1", that is a single character, either T or H
            # (\\d+) → Encapsulate "group 2", one ore more digits(numeric), \\d = digit
            # .* → any character and any lenght
            # REPLACEMENT
            # \\1 → Call encapsulated group 1
            # \\2 → Call encapsulated group 2
            # Nothing else since .# do not have a replacement, then it's dropped
            
      
      ## 2) Apply and check
      colnames(counts) <- simple             # Set the new names on counts
      if (any(duplicated(colnames(counts)))) stop("Duplicate sample names after renaming.")
            # A checkpoint, where if  "duplicated() = TRUE" we "stop" and message appears
      print(colnames(counts))
      head(counts)

      
      ## 3) Reorder to T1..T8 then H9..H16 (nice to have)
      pref <- substr(colnames(counts), 1, 1)                     # "T" or "H"
      num  <- as.numeric(sub("^[TH]", "", colnames(counts)))     # numeric index
            # ^ → anchor to the string start
            # [TH] → either T or H
            # "" → Delete [T or H]
            # as.numeric(..) → convert the remaining characters to numbers
      ord  <- order(factor(pref, levels=c("T","H")), num)
            # factor(pref, levels=c("T", "H")) → Forces T to sort before H, my decision not alphabetical
            # order(..., num) → sorts(clasificar) by a number from 1 to ...
      counts <- counts[, ord] 
            #Reindexes the coumps with that order
      print(colnames(counts))
      
      
      ## Ensure we have exactly the expected set name-info
      expected <- c(paste0("T", 1:8), paste0("H", 9:16))
            # paste0("T", 1:8) → "T1","T2",...,"T8"
            # paste0("H", 9:16) → "H9","H10",...,"H16"
      missing  <- setdiff(expected, colnames(counts))
      extra    <- setdiff(colnames(counts), expected)
            # setdiff(A, B) → gives items in A not in B
            # missing → samples you expected but don’t have
            # extra → unexpected columns present in counts
      if (length(missing)) stop("Missing: ", paste(missing, collapse=", "))
      if (length(extra))   stop("Extra: ",   paste(extra, collapse=", "))

### Zero-count filter

      # counts: your data frame with columns T1..T8, H9..H16
      
      # Drop genes with zero counts across all samples
      # Removes genes that are "completely absent" across "all samples"
      counts <- counts[rowSums(counts) > 0, , drop = FALSE]
      cat("After zero-row filter:", nrow(counts), "genes and", ncol(counts), "samples\n")
            # cat() → prints a quick summary

### Grouping into T0 & H24

When we want to se how RAW DATA behavior alone, say plot of the CPM/FPM of each gene in T1 to H16, grouping by condition T0 or H24 is not introducing BIAS, it is just labeling:
      
→ We aren't filtering based on the difference

→ We are "organizing" the samples

   → It become BIAS when we decide based on these groups which genes to keep

      ## Build colData for grouping (T0 vs H24) 
      samples <- colnames(counts)
      time <- ifelse(grepl("^T", samples), "T0", "H24")
            # grepl("^T", samples) → TRUE if the sample name starts with "T"
            # it TRUE the it assigns T0, else H24
      time <- factor(time, levels = c("T0", "H24"))
            # factor() → categorical variable
      coldata <- data.frame(time = time, row.names = samples)
            # data.frame() → Build up the metadata "table"

### FPM screen

In DESeq2 package, FPM = Fragments Per Million 

|                  CPM = Counts Per Million in edgeR & limma

<img width="300" height="150" alt="image" src="https://github.com/user-attachments/assets/ca0c1cfd-87b4-4ed5-a3f2-e6ae7d382c02" />

      ##  FPM >1 in ≥ half of reps per group -------------
      # HALF → We may have genes that are “off” in T0 but “on” in H24 (or vice versa)

      ## 1) Create DESeq 
      # DESeq2 needs integer counts; round if your file has decimals
      dds_raw <- DESeqDataSetFromMatrix(countData = as.matrix(round(counts)),
                                        colData = coldata,      # Provides metadata
                                        design = ~ time)        # Tells experimental design

      
      ## 2) Compute raw FPM/CPM
      # Raw FPM/CPM-like (before size-factor normalization) for *presence* filtering
      cpm_raw <- fpm(dds_raw, robust = FALSE)   # keep FALSE to match our pipeline
            # Here, it’s before normalization, quick way to check presence (expression >1)

            
      ## 3) Keep a gene if FPM>1 in at least half of the replicates in ANY group (T0 or H24)
      grp <- coldata$time
      mat_eval <- cpm_raw > 1                  # TRUE if FPM > 1, else FALSE.
      keep <- rep(FALSE, nrow(mat_eval))       
      cutoff <- ceiling(table(grp) / 2)        # half of replicates per group 


      ## 4) Group wise filtering loop 
      for (g in levels(grp)) {
        idx <- which(grp == g)
        keep <- keep | (rowSums(mat_eval[, idx, drop = FALSE]) >= cutoff[g])      keep = TRUE if the gene passed the rule in at least one group
      }
            # DO NOT introduce BIAS → since he OR (|) combines them: if a gene is expressed in ≥ half of T0 samples OR ≥ half of H24 samples, it’s kept.
            # Also the rule is symmetric in T0 & T24
            # NOT → keep genes that are higher in H24 than T0 = BIASED
            # YES → keep genes that are expressed in a reasonably consistent way in "at least one group"


      ## 5) Prints how many genes survived and the % of total.
      cat("Genes kept by FPM > 1 prevalence rule:", sum(keep), "out of", length(keep),
          sprintf("meaining a %.1f%%\n kept", 100 * sum(keep) / length(keep)))


      ## 6) Subset your objects for downstream steps
      counts_filt <- counts[keep, , drop = FALSE]
      dds_1 <- dds_raw[keep, ]  # we’ll use this for size-factor normalization next
            # Keeps only the filtered genes in both raw counts and DESeq2 object.

                  
      ## 7) Quick sanity check
      stopifnot(identical(colnames(counts_filt), rownames(coldata)))
            #Ensures sample names still match perfectly between counts_filt and coldata

### DESeq2 normalization

Size factor nomalization

→ Corrects for differences in sequencing depth / library size between samples

E.g. one sample may have 15M reads, another 25M → without correction, all counts in the 25M sample would look “higher” just due to depth 

      ## Size-factor normalization (DESeq2)
      dds_1 <- estimateSizeFactors(dds_1)      # median-of-ratios
      sizeFactors(dds_1)[1:5]                  # quick peek(vistazo)

<img width="309" height="30" alt="image" src="https://github.com/user-attachments/assets/b0426a2a-1de0-47e0-b9bc-40bd172d43a4" />

So:

→ T3 had the largest depth (1.21)

→ T5 is almost average (1.00)

### VST - Variance Stabilizing Transform

<a name="vst_1"></a>
vst → applies DESeq2’s variance-stabilizing transformation

"Raw counts" follow a negative binomial distribution with variance dependent on mean, what is twhy that makes PCA and heatmaps messy

→ Transform this "Raw counts" to the log2-like scale but also stabilizes variance so that genes with high counts don’t dominate

      ## Variance Stabilizing Transform (VST) ---------
      vst_1   <- vst(dds_1)
      vst_mat <- assay(vst_1)                  # genes x samples = matrix

<img width="383" height="100" alt="image" src="https://github.com/user-attachments/assets/10d6a18f-c3af-45c6-a31f-7b5a525b9453" />


###  Z-score standardization

Standardization across samples, With this each gene now has mean=0, sd=1 across samples

→ Makes it easier to compare patterns of expression across genes

      ## 1)Gene-wise z-score -----------
      # (center 0, sd 1 across samples)
      zmat <- t(scale(t(vst_mat)))             # rows = genes, cols = samples
            # scale() standardizes each column of a matrix by default (mean=0, sd=1)
            # t(vst_mat) → transpose (samples × genes) and vice versa

        
      ## 2) Sanity: mean≈0, sd≈1 per gene
      stopifnot(all(is.finite(zmat)))            # Ensures no NA, NaN, or Inf in the matrix
      
      ## 3) mean should be very close to 0
      cat("Mean of first 5 genes after z-score:\n"); print(rowMeans(zmat[1:5, ]))
      cat("SD of first 5 genes after z-score:\n"); print(apply(zmat[1:5, ], 1, sd))

<img width="434" height="85" alt="image" src="https://github.com/user-attachments/assets/45c3156e-244c-4341-9346-0faf1d7edea7" />

### PCA → Principal Component Analysis

prcomp() → is R’s function for PCA

Varirables:
<a name="pca"></a>
- pca     
<a name="var_expl"></a>
- var_expl

      ## PCA on samples (using prcomp) --------
      
      # We already z-scored by gene, so no extra scaling here
      # Input must be samples × features
      pca <- prcomp(t(zmat), center = FALSE, scale. = FALSE)        # observations = samples
            # center = FALSE → we don’t need to mean-center features, because z-score already did that
            # scale. = FALSE → we don’t need to scale to sd = 1, because z-score already did that too.
      var_expl <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
            # Each PC has a variance = pca$sdev^2.
            # Divide by total variance to get % explained by that PC.
            # Multiply by 100 → nice percentages.

<img width="383" height="55" alt="image" src="https://github.com/user-attachments/assets/8b72e887-f9f7-45ee-b64a-b70c6d927a43" />

### PCA colors

      library(RColorBrewer)   #Color
      display.brewer.all()

<img width="220" height="330" alt="image" src="https://github.com/user-attachments/assets/339f577b-d3e8-4261-824f-8a557d9b51e0" />

We will use Dark2

      # Colors and legend (robust) 
      pal_time  <- c(T0 = "#79cbb8", H24 = "#500472")                  # Hex fortmat
      #pal_time <- setNames(brewer.pal(2, "Dark2"), c("T0", "H24"))    # RColorBrewer
      col_time  <- pal_time[as.character(coldata$time)]
            # coldata$time → was the factor with values "T0" or "H24"
      leg_labs  <- unique(as.character(coldata$time))
            # unique(...) → collects unique condition labels (E.g. "T0", "H24")
      leg_labs  <- leg_labs[!is.na(leg_labs)]
            # removes any possible NA

### PCA sample plots

Principal Component Analysis (PCA) is a dimensionality reduction technique

In RNA-seq:

→ Each sample (T1…H16) has thousands of gene expression values (a high-dimensional space)

→ PCA finds new axes (principal components) that capture the largest sources of variation in the data

- PC1 = axis that explains the largest % of total variance.

- PC2 = second largest, orthogonal to PC1.

→ Plotting samples on PC1–PC2 shows which samples are more similar or different

**Variables for DIMENSIONS**

<a name="var_1"></a>
- var_1
<a name="var_2"></a>
- var_2

      ## PCA plots -----------
      
      ## 1) Pick the PCs you want to plot
      var_1 <- 1
      var_2 <- 3
      
      ## 2) Build a small data frame for ggplot
      df_pca <- data.frame(
        x     = pca$x[, var_1],
        y     = pca$x[, var_2],
        name  = rownames(pca$x),         # "T1".."H16"
        time  = factor(coldata$time, levels = c("T0","H24"))  # ensure order
      )
                  # pca$x → coordinates of each sample in PCA space
                  # coldata$time → group each sample into condition/time ("T0"/"H24")
      
      ## 3) PLOT
      ggplot(df_pca, aes(x, y, color = time, label = name)) +
        geom_point(size = 2.6) +            # Geometry and size dot
        ggrepel::geom_text_repel(show.legend = FALSE, max.overlaps = Inf, size = 3) +     # Size of each dot
              # geom_text_repel (from ggrepel) pushes text labels away from each other and from the points
        scale_color_manual(
          values = pal_time,                # <- uses your exact palette
          limits = c("T0","H24"),           # <- fixes legend order & ensures mapping
          name   = "Time"                   # <- legend title to match base plot
        ) +
        labs(
          title = sprintf("PCA – Samples → PC%d vs PC%d", var_1, var_2),
          x = sprintf("PC%d (%.1f%%)", var_1, var_expl[var_1]),
          y = sprintf("PC%d (%.1f%%)", var_2, var_expl[var_2])
        ) +
        theme_classic(base_size = 12) +     # Relación gráfico vs títulos 
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5),      # Moves tittle
                  # theme_classic() → removes grey background, gives a cleaner look
          legend.position = "right"
        )
      
      ## 4) save PCA matrices
      # write.csv(pca$x, file = "PCA_samples_scores.csv")      # samples in PC space
      # write.csv(pca$rotation, file = "PCA_gene_loadings.csv") # genes loadings

<img width="410" height="350" alt="image" src="https://github.com/user-attachments/assets/98d097fd-54e4-4358-851a-d6b934426d4f" />

### PCA gene plot

      ## PCA gene plot  ------
      ## 1) extract loadings (genes × PCs)
      loadings <- pca$rotation   # each row = a gene, each column = a PC

            
      ## 2) build a dataframe for ggplot
      df_load <- data.frame(
        gene = rownames(loadings),
        PCx  = loadings[, var_1],
        PCy  = loadings[, var_2]
      )
      
      ## 3) ggplot
      ggplot(df_load, aes(x = PCx, y = PCy)) +
        geom_point(alpha = 0.4, size = 0.6, color = "#AF9AB2") +
        labs(
          title = sprintf("PCA – Genes → PC%d vs PC%d", var_1, var_2),
          x = sprintf("PC%d loadings (%.1f%%)", var_1, var_expl[var_1]),
          y = sprintf("PC%d loadings (%.1f%%)", var_2, var_expl[var_2])
        ) +
        theme_classic(base_size = 12) +     # Relación gráfico vs títulos 
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5),      # Moves tittle
          legend.position = "right"
        )

<img width="410" height="350" alt="image" src="https://github.com/user-attachments/assets/3ca56c97-c4d1-47e8-88a2-6c53cd76b39b" />

### PCA color gene plots

*Missing the aesthetics and tittles

**INCOMPLETE**

      ggplot(df_load, aes(x = PCx, y = PCy, color = abs(PCx))) +
        geom_point(alpha = 0.6, size = 0.6) +
        scale_color_viridis_c() +
        theme_classic()

<img width="350" height="350" alt="image" src="https://github.com/user-attachments/assets/ecb5dbc0-73ec-48bd-aee7-691fba9ed119" />


### PCA for TOP genes

Variable that we already defined:

- [pca](#pca)
- [var_expl](#var_expl)
- [var_1](#var_1)
- [var_2](#var_2)

→ We're no longer plotting samples (PCA scores, pca$x) → These are the coordinates of each sample on the new PC axes

→ Ww’re plotting **genes using their PCA loadings** (pca$rotation) → the weights that define each PC axis → tell us which genes contribute the most on the separation we saw in the sample plot

→ Everything we did (zero-count filter → FPM > 1 prevalence → sizefactor normalization → VST → z-score per gene) was applied to the gene expression matrix (genes × samples) **That’s the input for PCA**

Meaning: Both samples and genes "see" the same preprocessed data → counts have been normalized, variance-stabilized, and z-scored across samples

> [!TIP]
> SAMPLES (scores) projection of normalized expression values for each sample → **Do T0 and H24 cluster apart?**
> 
> GENES (loadings) coefficients that define how much each gene “pushes” a sample along each PC axis → **Which genes made T0 vs H24 separate on PC1 (positive vs negative side)?**

*scores = where each sample lies in PC space

*loadings = how much each gene contributes to each PC


      ## PCA TOP genes ----
      
      # Inputs I already have
      # pca, var_expl, var_1, var_2
      gene_loadings <- pca$rotation  # safer name than 'loadings()'
      
      ## 1) Top K genes by sign (per side)
      K <- 100
      pos_genes <- names(sort(gene_loadings[, var_1], decreasing = TRUE))[1:K]
      neg_genes <- names(sort(gene_loadings[, var_1], decreasing = FALSE))[1:K]
      
      ## 2) Build df's with a 3-level flag
      df_load <- data.frame(
        gene = rownames(gene_loadings),
        PCx  = gene_loadings[, var_1],            # For each gene, store its loading on PC var_1 (x)
        PCy  = gene_loadings[, var_2],            # For each gene, store its loading on PC var_2 (y)
        group = "other"
      )
      ## 2.1) Tag genes as "top_pos", "top_neg", or "other" for coloring.
      df_load$group[df_load$gene %in% pos_genes] <- "top_pos"
      df_load$group[df_load$gene %in% neg_genes] <- "top_neg"
      df_load$group <- factor(df_load$group, levels = c("other","top_pos","top_neg"))
      
      ## 3) Colors (hex ok)
      col_other <- "#AF9AB2"  
      col_pos   <- "#f9665e"  # red
      col_neg   <- "#799fcb"  # blue  
      
      ## 4) Plotting
      ggplot(df_load, aes(PCx, PCy, color = group)) +
        geom_point(alpha = 0.65, size = 0.8) +
        scale_color_manual(
          values = c(other = col_other, top_pos = col_pos, top_neg = col_neg), # Color of dots
          labels = c(
            other   = "Other genes",
            top_pos = sprintf("Top %d positive on PC%d", K, var_1),
            top_neg = sprintf("Top %d negative on PC%d", K, var_1)
          ),
          name = "Highlight"
        ) +
        labs(
          title = sprintf("PCA – TOP Genes → PC%d vs PC%d", var_1, var_2),
          x = sprintf("PC%d loadings (%.1f%%)", var_1, var_expl[var_1]),
          y = sprintf("PC%d loadings (%.1f%%)", var_2, var_expl[var_2])
        ) +
        theme_classic(base_size = 12) +     # Relación gráfico vs títulos 
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5),      # Moves tittle
          legend.position = "top"
        )
      
      ## 5) The genes on CVS form
      write.csv(data.frame(gene = pos_genes,
                           loading = gene_loadings[pos_genes, var_1],
                           direction = "positive",
                           PC = paste0("PC", var_1)),
                sprintf("top_%d_positive_genes_PC%d.csv", K, var_1), row.names = FALSE)
      
      write.csv(data.frame(gene = neg_genes,
                           loading = gene_loadings[neg_genes, var_1],
                           direction = "negative",
                           PC = paste0("PC", var_1)),
                sprintf("top_%d_negative_genes_PC%d.csv", K, var_1), row.names = FALSE)

<img width="370" height="370" alt="image" src="https://github.com/user-attachments/assets/882f4eb8-7e2f-491a-9c42-a1fa06291b7b" />

> [!WARNING]
> 
> Sign is arbitrary
>
> With this plot **WE CAN'T SAY** “positive loadings = upregulated genes, negative loadings = repressed genes”
>
> Whether a gene appears “positive” or “negative” on PC1 can flip if the algorithm decides the opposite orientation
>
> What matters is **relative separation** positive vs negative genes, and how that aligns with your samples on plot, more to the right or to the left
>
> Genes with large +PC1 loadings are more associated with H24
>
>Genes with large –PC1 loadings are more associated with T0

### Heatmap Sample x Sample

Variable that we already defined:

- [vst_1](#vst_1)

      ## UNSUPERVISED CLUSTERING & HEATMAPS ----

      
      ## 1) Feedback & Fallback (only if vst_1/vst_mat don't exist) →  Ensuring VST was made
      if (!exists("vst_1")) {
        message("vst_1 not found, computing VST from dds_1...")
        stopifnot(exists("dds_1"))
        vst_1 <- vst(dds_1)
      }
      if (!exists("vst_mat")) vst_mat <- assay(vst_1)   # genes x samples
         # assay(vst1) → numeric matrix with rows = genes, cols = samples

      ## 2) SAMPLE–SAMPLE CLUSTERING -----------
      # Use Pearson correlation between samples (robust for RNA-seq after VST)
      cor_mat <- cor(vst_mat, method = "pearson")           # samples x samples
      ann_col <- data.frame(time = coldata$time)
          # Prepares sample annotations (here just time), to color the heatmap margins
      rownames(ann_col) <- rownames(coldata)
      stopifnot(identical(colnames(cor_mat), rownames(ann_col)))
          # Ensures the column order of the matrix matches the row order of the annotations
      
      
      ## 3) Heatmap of correlations sample -sample
      pheatmap(cor_mat,
               color            = colorRampPalette(rev(brewer.pal(12, "BrBG")))(255),  # (Darkness, "Color")
               border_color     = NA,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
                        # Clustering → Ward’s method on Euclidean distance applied to the correlation matrix
               clustering_method = "ward.D2",
               display_numbers  = FALSE,
               annotation_col   = ann_col,
               annotation_row   = ann_col,
               annotation_colors = list(time = pal_time),
               main = "Sample–sample correlation")
  
<img width="450" height="370" alt="image" src="https://github.com/user-attachments/assets/40b537d2-c99d-487b-a291-64cf9d89c1ad" />

### Dendogram

Hierarchical clustering of samples (or genes), where the height of the branches reflects the distance (or dissimilarity)

      ## Dendogram ------

      ## 1) Compute distance 
      dist_samples <- dist(t(vst_mat), method = "euclidean")                  # correlation distance
      # dist_samples <- dist(t(vst_mat[,c(-7,-8, -9, -15, -16)]), method = "euclidean") 
      
      ## 2) Hierarchical clustering
      hc_samples <- hclust(dist_samples, method = "ward.D2")       # method that minimizes variance within cluster
      
      ## 3) Plot dendrogram
      plot(hc_samples, hang = -1, main = "Sample clustering (Euclidean + Ward.D2)")
            # hang = -1 → forces all sample labels to sit at
      
      ## 4)  helper: height that produces exactly k clusters
      k <- 8      # Here we change the number of clusters we want
      cut_height_for_k <- function(hc, k) {
        n <- length(hc$height)
        stopifnot(k >= 2, k <= n)  # k = 1 means no cut
              # mid-point between the two merges that separate into k clusters
        (hc$height[n - k + 1] + ifelse(n - k + 2 <= n, hc$height[n - k + 2], Inf)) / 2
      }            # hc$height → stores the merge heights at which clusters combine

      # It figures out where to draw the red horizontal line to yield exactly k clusters
                  # Trick: for k clusters, we cut between the (k-1)-th and k-th largest merges
                        # (k-1)-th → largest merge height = the step where the tree goes from k clusters to k-1
                        #  k-th → largest merge height = the step just before that.

      n <- length(hc_samples$height)
      h_k <- mean(hc_samples$height[c(n - k + 1, n - k + 2)])
                  # The height computed for chosen k
      
      
      # 5) Add a horizontal cut line (example height = 150, adjust as needed)
      abline(h = h_k + -1, col="red", lty=2, lwd=2)       # We added a "-1" just to do not cover the other lines
      
      # 6) Alternatively, cut into a fixed number of clusters (e.g., k = 2)
      rect.hclust(hc_samples, k = k, border = "blue")  # draw cluster boxes
      
      ## 7) cluster summary tables membership
      grp <- cutree(hc_samples, k = k)
      print(grp)               # sample -> cluster id
      table(grp)               # cluster sizes

<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/a1fc760a-bf8f-466a-9d6c-bae385b8ab3e" />

### Heatmap Gene x Sample

- VST (Variance Stabilizing Transform) → removes dependence of variance on mean counts, making genes comparable across samples

- Row-z (row scaling) → each gene’s expression values are z-scored (mean = 0, sd = 1 across samples) so patterns are relative per gene, not absolute counts.

- Correlation distance (1 – correlation) → similarity between samples is measured by **how similar their expression profiles are** (ignores absolute magnitude, focuses on shape).

<a name="Ward.D2"></a>
- Ward.D2 clustering method → groups samples/genes by minimizing variance within clusters, producing balanced dendrogram branches.

- Top variable genes → selects the genes with the **highest variance across** samples, which are most informative for clustering patterns.

      ## HEATMAP – TOP VARIABLE GENES -------------
      
      ## 1) Select top N most variable genes across samples (after VST)
      topN <- 1000
      gene_var <- rowVars(vst_mat)
      ord <- order(gene_var, decreasing = TRUE)
      top_idx <- ord[seq_len(min(topN, length(ord)))]
      mat_top <- vst_mat[top_idx, , drop = FALSE]
      
      ## 2) Heatmap (scale by row = z-score per gene)
      pheatmap(mat_top,
               scale            = "row",                  # z-score per gene
               color            = colorRampPalette(brewer.pal(12, "RdBu"))(255),
               border_color     = NA,
               show_rownames    = FALSE,
               show_colnames    = TRUE,
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",
               clustering_method = "ward.D2",
               annotation_col   = ann_col,
               annotation_colors = list(time = pal_time),
               main = sprintf("Top %d variable genes (VST, row-z)", nrow(mat_top)))
      
      ## 3)  SAVE TO FILES -------------------
      # pdf("sample_correlation_heatmap.pdf", width = 7, height = 6); <repeat heatmap A>; dev.off()
      # pdf("top_variable_genes_heatmap.pdf", width = 8, height = 8); <repeat heatmap B>; dev.off()

<img width="450" height="350" alt="image" src="https://github.com/user-attachments/assets/180be73a-431a-40dd-a7d9-a4b44ad708c4" />





## DEGs in TIME CONTRAST T0 vs T24 CODE

### ↑↓ DEGs in time TEST

DESeq2 to find **differentially expressed genes** (DEGs)

> [!NOTE]
> Use of:
> 
> **LRT** → Likelihood Ratio Test (LRT), First step
>
> **FDR** → False Discovery Rate, Fourth step

Uses a **Likelihood Ratio Test (LRT)** to test the **effect of time (T0 vs H24)**. I’ll go line-by-line and add the “why”

The LRT asks: does the full model (~ time (LEVELS)) explain the counts significantly better than the null (~ 1)?

→ In other words: is there any effect of time on gene expression?

the LRT p-value is essentially testing whether the time coefficient is non-zero. With > 2 levels, LRT is preferred because it tests the overall effect across all levels

→ **HOW FDR works?**

When you test thousands of genes, many will look “significant” just by chance: 

- A raw p-value < 0.05 would mean ~5% false positives per test.

- With 20,000 genes, that’s ~1,000 false positives!

FDR controls the expected proportion of false positives among the genes you call significant

E.g. FDR < 0.01 means: among the genes you flag as DEGs, on average only ~1% are expected to be false discoveries

> [!CAUTION]
> FIRST IMPORTANT DOCUMENT
> 
> Differential expression results → Where we identify which genes have differential expression
>
> **"This gene is up-regulated at 24h relative to baseline T0"**

- log2FC > 0 → gene is globally up-regulated at H24 compared to T0

- log2FC < 0 → gene is globally down-regulated at H24 compared to T0

      ## DEGs with LRT (time effect) -----------

      ## 1) Fit the full (~ time) vs reduced (~ 1) model with Likelihood Ratio Test
      dds_1 <- DESeq(dds_1, test = "LRT", reduced = ~ 1)
             # dds_1 → filtered DESeqDataSet with design ~ time (LEVELS T0 & T24)
             # test = "LRT" → runs nbinomLRT (likelihood ratio test)
             # reduced = ~ 1 means the reduced (null) model has no covariates (intercept only)
             
      ## 2) Pull results (by default, coefficient is H24 vs T0)
      res_lrt <- results(dds_1)
            # Returns a DESeqResults table with one row per gene → log2FoldChange, lfcSE, stat, pvalue, **padj (FDR-adjusted p, FDR-corrected using Benjamini–Hochberg)**

      ## 3)  Replace NA padj with 1 (non-significant) 
      res_lrt$padj[is.na(res_lrt$padj)] <- 1
            # Treat NAs as non-significant

      ## 4) Call significant genes by FDR threshold
      alpha <- 0.01
            # alpha → FDR cutoff (1% here, false positives among the genes you call significant)
      deg_mask <- res_lrt$padj < alpha      # boolean vector marking significant genes.
      cat("DEGs at FDR <", alpha, ":", sum(deg_mask), "genes\n")       # Prints how many DEGs we have at that FDR
      
      ## 5) Quick summary table and top rows
      print(summary(res_lrt))
      head_res <- as.data.frame(res_lrt[order(res_lrt$padj), ])[1:10, ]
      cat("\nTop 10 by FDR:\n"); print(head_res)
        
      ## 6) Save complete LRT table 
      # write.csv(as.data.frame(res_lrt), file = "DEG_LRT_time_vs_null.csv")
            # This is the "Differential expression results"
                      # Each row means a gene, and it have:
                                  # log2 fold change (H24 vs T0, by default)
                                     # p-value and adjusted p-value (FDR)
                                     # This is where WE identify which genes are significant.

<img width="550" height="100" alt="image" src="https://github.com/user-attachments/assets/29870254-36aa-4e6b-885a-160ea3c7217f" />

- 1 → Started with 26,083 genes that passed the earlier zero-count filter

- 2 → Sets the default FDR threshold (padj < 0.1) **this is only for reporting here — in the final analysis, say: set alpha <- 0.01 (stricter)**

- 3 → Out of the 26,083 genes, 5427 genes (21%) **had a positive log2FoldChange (higher in H24 relative to T0) and were SIGNIFICANT** (alpha <- 0.01)

- 4 → Out of the 26,083 genes, 4881 genes (19%) **had a negative log2FoldChange (lower in H24 relative to T0) and were SIGNIFICANT** (alpha <- 0.01)

- 5 → DESeq2 flags genes as outliers if a single sample drives extreme counts (using Cook’s distance)

- 6 → Genes with very low counts across all samples are often filtered automatically

- 7 → This is the threshold used for "low counts"

- 8 → control outlier filtering via cooksCutoff.

- 9 → tweak filtering of low-information genes with independentFiltering

<img width="470" height="100" alt="image" src="https://github.com/user-attachments/assets/b70ec4aa-d552-4abc-a913-b9bdb6f57750" />

### DEGs for pattern discovery

This allows us to prepare data so that genes are comparable across samples and patterns **(clusters, SOM, heatmaps)** highlight relative up/down regulation, not absolute counts

→ We normalized, log2-transformed, and z-scored those DEGs (T0 vs 24) only.

> [!IMPORTANT]
>
> **→ BIOLOGICAL MEANING**
>
> FPM → corrects for sequencing depth and library size = compare gene X in sample T3 vs H10 fairly, because differences aren’t just due to different read totals
>
> LOG2 +1 → Expression changes become symmetric and interpretable. E.g. A jump from 10 → 20 counts = log2 fold-change of +1 & A drop from 20 → 10 counts = log2 fold-change of –1
>
> Gene-wise z-score → Each gene is centered to its own baseline → mean = 0, sd = 1. E.g. Positive = gene is up-regulated in that sample compared to its own average; Negative = gene is down-regulated
>
> Now clustering highlights **patterns of co-regulation** → Genes that rise/fall together across samples cluster into groups. Samples with similar overall expression profiles cluster together


> [!CAUTION]
> SECOND IMPORTANT DOCUMENT
> 
> Same DEGs as in the LRT document, but transformed into standardized values (s-scores) for pattern discovery (is the input for PCA, clustering, SOM, heatmaps, etc.)
      
      ## Transform DEGs for pattern discovery -------
      
      ## 0) make sure we actually have DEGs at FDR < 0.01
      if (sum(deg_mask) == 0) {
                # deg_mask is the logical vector from past step where "TRUE" if gene is significant at FDR < 0.01
        stop("No DEGs at FDR < 0.01. Revisit step ↑↓ DEGs in time TEST or use a different alpha temporarily")
      }       # If no genes are significant, it stops the script
        
      
      ## 1) Get normalized FPM (uses the size factors we estimated earlier)
      fpm_norm <- fpm(dds_1)   # normalized CPM-like values
      DEG_cpm  <- fpm_norm[deg_mask, , drop = FALSE]
            # Then we subset only DEGs (rows where deg_mask == TRUE)
      cat("DEG_cpm matrix:", nrow(DEG_cpm), "genes x", ncol(DEG_cpm), "samples\n")
            # We obtain a matrix → rows = significant genes, columns = sample
      
      ## 2) Log2(+1)
      log2_deg <- log2(DEG_cpm + 1)
            # log2 → stabilizes variance and makes fold-changes symmetric
            # The +1 avoids log(0)
      
      ## 3) Gene-wise z-score (center 0, sd 1 across samples)
      cs.log2 <- t(scale(t(log2_deg)))
            # scale() → standardizes each row to mean = 0, sd = 1, y lo transloca (wrap), change rows with columns
      
      ## 4) Quick sanity checks
      stopifnot(all(is.finite(cs.log2)))      # Ensures no missing/infinite values
      cat("z-score range:", sprintf("[%.2f, %.2f]\n", min(cs.log2), max(cs.log2)))
                                              # Prints the z-score range (should be roughly [-3, +3])
      stopifnot(identical(colnames(cs.log2), rownames(coldata)))
                                              # Double-checks that samples in matrix = metadata samples
      
      ## 5) save intermediates
      # write.csv(DEG_cpm,  "DEGs_FPM_normalized.csv")
      # write.csv(cs.log2,  "DEGs_log2p1_zscore.csv")
            # a processed expression matrix of only the DEGs that passed your FDR cutoff
            # Rows = significant genes
            # Columns = samples
            # Values → z-scores (relative expression; up = positive, down = negative, per gene).

### SOMs

Self Organizing Maps

→ We trained the self-organizing map using only those (TO vs T24) DEGs

Here we move from PCA-style linear analysis to a **nonlinear clustering method** that groups genes with similar expression dynamics into **"neurons"** on a grid

Empty neurons aren't bad, simply mean “no genes fit this prototype"; and also prevents forcing genes into artificial clusters that don’t match them

      ## Self-Organizing Map (SOM) --------
      
      ## 1) Sanity check
      stopifnot(exists("cs.log2"), nrow(cs.log2) > 0)
          # Ensures we have the input matrix cs.log2 (DEGs log2-normalized + z-scored)
      
      
      ## 2) Randomness fixing
      set.seed(42)
          # SOM training involves random initialization, so this fixes randomness
      
      
      ## 3) SOM grid definition
      som_grid_axes = 5       # The lenght in neurons as axis
      som_grid <- kohonen::somgrid(xdim = som_grid_axes, ydim = som_grid_axes, topo = "hexagonal", toroidal = TRUE)
                              # xdim = 5, ydim = 5 → a 5 × 5 grid = 25 neurons
                              # topo = "hexagonal" → neurons arranged in a honeycomb pattern (each neuron has 6 neighbors)
                              # toroidal = TRUE → edges wrap around (like donut), so no neuron is isolated on the edge
      
      ## 4) Training of the SOM
      som_fit <- kohonen::som(
        as.matrix(cs.log2),        # data: genes × samples, z-scored
        grid      = som_grid,      # Step 3
        rlen      = 200,           # training iterations
        keep.data = TRUE           # keep the original data mapped
      )
          # Input → each gene’s profile across samples
          # Output → a map of neurons where similar genes are grouped
          # Each gene gets mapped to a "best matching unit" (BMU = a neuron)

      
      ## 5) Gene classification & codebook vectors
      classif <- som_fit$unit.classif       # assigns each gene → neuron ID (which neuron the gene belong)
      names(classif) <- rownames(cs.log2)   # 
      codes <- som_fit$codes[[1]]           # matrix of prototype expression patterns (25 neurons × samples)
                  # It's basically each neuron’s "average gene behavior"
      
      ## 6) Ploting
      
      # CODEBOOK VECTORS → Shows the prototype expression patterns of each neuron
      plot(som_fit, main = "SOM – Codebook vectors")
      
      # MAPPING → Each gene is a dot placed on the neuron it belongs
      plot(som_fit, type = "mapping", pch = ".", main = "SOM – Mapping")
      
      # TRAINING PROGRESS → Learning curve of SOM, error reduction over iterations
      plot(som_fit, type = "changes", pch = ".", main = "SOM – Training")
      
      # COUNTS → Similar to an Histogram → How many genes were asigned at each neuron
      plot(som_fit, type = "counts", pch = ".", main = "SOM – Genes per neuron")
      
      # U-MATRIX → Distance between neurons, like cluster of clusters
      plot(som_fit, type = "dist.neighbours", pch = ".", main = "SOM – U-matrix")
      
      # QUANTIZATION ERROR → low = neuron well represents its genes, high error = more heterogeneous
      plot(som_fit, type = "quality", pch = ".", main = "SOM – Quantization error")
  
<img width="250" height="300" alt="image" src="https://github.com/user-attachments/assets/cd109ace-445a-4a69-ad77-1b6c76c5a475" />
<img width="250" height="300" alt="image" src="https://github.com/user-attachments/assets/88e8dba5-2f2b-41b4-800f-974331f1414e" />
<img width="250" height="300" alt="image" src="https://github.com/user-attachments/assets/e7129bf2-13a9-4c38-9c82-25516268bf5f" />
<img width="250" height="300" alt="image" src="https://github.com/user-attachments/assets/86effcac-b3e2-40e5-8cc2-6a3bcac8809c" />
<img width="250" height="300" alt="image" src="https://github.com/user-attachments/assets/21a14f14-74d2-4959-8781-9956ebabdcfd" />
<img width="250" height="300" alt="image" src="https://github.com/user-attachments/assets/472e22b6-ab5c-457a-a78f-9d837e78e3c1" />



### Heatmaps of SOM

We are clusterizing and visualizing the SOM neurons/codebooks that summarize DEG (T0 vs T24) patterns

Concepts:

[Ward.D2](#Ward.D2) →  busca las dos agrupaciones que, al unirse, aumenten lo menos posible la varianza interna

**DISTANCES**

- **Euclidean** → For samples (Samples which have similar globally transcriptome)

Measure the straight-line distance  

Sensitive to "magnitude" and overall level of the vectors

Doesn't matter if few genes change too much, the important is the GLOBAL TENDENCY


- **Manhattan** → For neurons and SOM codeboooks (Pattern of expression of many genes)

Measure the sum of absolute differences

We want to compare "patterns" not magnitude

Each neuron is a “prototype” of gene behavior

      ## Cluster codebooks + heatmaps --------
        
      ## 1) Info calling; codes → neurons x samples (from Step 9)
      codes <- som_fit$codes[[1]]   # Same matrix of prototype expression patterns (25 neurons × samples)
      stopifnot(identical(colnames(codes), rownames(coldata)))  # Chekpoint to see if samples align with coldata
                
      
      ## 2) Cluster samples (dendograms) (rows in the heatmap) and neurons (cols) 
      clust.sample <- hclust(dist(t(codes), method = "euclidean"), method = "ward.D2")
      clust.neuron <- hclust(dist(codes,   method = "manhattan"), method = "ward.D2")
      
      
      ## 3) Row annotations (samples) + colors
      ann_row <- data.frame(time = coldata$time)
      rownames(ann_row) <- rownames(coldata)
      ann_colors <- list(time = c(T0 = "#b89d47", H24 = "#fe9179"))  # named palette
      
      
      ## 4) Ploting → Heatmap 1: SOM codebooks (samples × neurons) 
      pheatmap(t(codes),
               border_color  = "grey60",
               scale         = "column",          # standardize each neuron across samples by z-score
               show_rownames = TRUE,
               show_colnames = FALSE,
               cluster_rows  = clust.sample,      # use our dendrogram for samples
               cutree_rows   = 2,
               annotation_row = ann_row,
               cluster_cols  = clust.neuron,      # dendrogram for neurons
               cutree_cols   = 3,
               annotation_colors = ann_colors,
               main = "SOM codebooks — samples × neurons")
      
      ## 5) Neuron clusters (at k clusters)
      k_neuron  <- 4                              # Aquí se cámbia 
      clust.pat <- cutree(clust.neuron, k = k_neuron)        # vector: nombre_neurona -> 1..k
      clust.aux <- paste0("C", clust.pat)                    # "C1","C2",...,"Ck"
      names(clust.aux) <- names(clust.pat)
      
      
      ## 6) Column annotation (neurons)
      ann_col <- data.frame(neuron_cluster = clust.aux)
      rownames(ann_col) <- names(clust.aux)
      
      
      ## 7) Color palette 
      pal_nc <- colorRampPalette(brewer.pal(8, "Set2"))(k_neuron)  # generates k colors
      names(pal_nc) <- paste0("C", seq_len(k_neuron))
      
      
      ## 8) Colores OF 'time' 
      ann_colors <- list(
        time = c(T0 = "#b89d47", H24 = "#fe9179"),
        neuron_cluster = pal_nc
      )
      
      
      ## 9) VERY IMPORTANT → forcing the factor to have the exact levels of palette
      ann_col$neuron_cluster <- factor(ann_col$neuron_cluster,
                                       levels = names(ann_colors$neuron_cluster))
      
      
      ## 10) Heatmap 2: with neuron-cluster colors on columns
      pheatmap(t(codes),
               border_color  = "grey60",
               scale         = "column",
               show_rownames = TRUE,
               show_colnames = FALSE,
               cluster_rows  = clust.sample,      cutree_rows = 2,
               annotation_row = ann_row,
               cluster_cols  = clust.neuron,      cutree_cols = k_neuron,
               annotation_col = ann_col,
               annotation_colors = ann_colors,
               main = "SOM codebooks — with neuron clusters")
      
      
      ## 11) Map each gene to a neuron cluster (pattern type)
      type.pattern <- clust.pat[classif]          # classif → Which genes belongs to which neuron
      names(type.pattern) <- names(classif)
      cat("Genes per neuron cluster:\n"); print(table(type.pattern))
      
      # (Optional) save outputs
      # write.csv(codes, file = "SOM_codebooks_neuron_by_sample.csv")
      # write.csv(data.frame(gene = names(type.pattern),
      #                      neuron_cluster = paste0("C", type.pattern)),
      #           file = "genes_to_SOM_clusters.csv", row.names = FALSE)

<img width="250" height="200" alt="image" src="https://github.com/user-attachments/assets/e89137da-948f-4053-a3ce-9742481c772a" />

<img width="250" height="200" alt="image" src="https://github.com/user-attachments/assets/6c89348b-0447-4000-9ff9-5d3e321140e1" />

<img width="300" height="30" alt="image" src="https://github.com/user-attachments/assets/6ca1fabb-39a6-4f88-9f36-4b0cb4904d98" />

## Annotation & Enrichment 0 vs 24       



*Este script fue obtenido y modificado del código realizado por el Dr. Carlos, se detallan las modificaciones realizadas y el why of the changes/adaptations:

- counts  → ya teníamos counts leído, filtrado por ceros, y lo más importante es que el objeto dds_1 con el diseño, los size factors, etc

> [!CAUTION]
> THIRD IMPORTANT DOCUMENT
> 
> ALL ENRICHMENT MODULE, PATHWAYS (ko or map) 41 documents

Full anotation document ()

      ## === 1) Universo y lista de DEGs coherentes con nuestra tubería ===
      universe_genes <- rownames(vst_mat)              # o rownames(dds_1)
      deg_genes      <- rownames(DEG_cpm)              # los significativos (FDR<0.01)
      
      ## === 2) Lee anotaciones ===
      ann <- read.delim("fullAnnotation.tsv.txt",
                        stringsAsFactors = FALSE,
                        row.names = 1)
      
      ## Helper para tokenizar columnas multivalor separadas por coma
      tokenize <- function(x) {
        x <- x[!is.na(x)]
        toks <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
        toks <- trimws(toks)
        toks[toks != ""]
      }
      
      ## === 3) Función ORA (Over-Representation Analysis) por columna ===
      run_ora <- function(category_col, universe, hits, ann_df) {
        # universo y hits presentes en anotación
        bg_vals  <- ann_df[intersect(universe, rownames(ann_df)), category_col]
        hit_vals <- ann_df[intersect(hits,     rownames(ann_df)), category_col]
        
        pop      <- table(tokenize(bg_vals))          # población por término
        aux      <- table(tokenize(hit_vals))         # hits por término
        if (length(pop) == 0 || length(aux) == 0) {
          return(list(p=structure(numeric(0)), fdr=structure(numeric(0))))
        }
        
        # vector de hits alineado a pop
        hit <- setNames(integer(length(pop)), names(pop))
        hit[names(aux)] <- as.integer(aux)
        
        # parámetros hipergeométrica
        M <- sum(pop)                                 # total etiquetas en universo
        n <- sum(hit)                                 # total etiquetas en hits
        
        pvals <- setNames(rep(1, length(pop)), names(pop))
        for (k in names(pop)) {
          # P(X >= hit[k]) con X~Hiper(M, pop[k], n)
          pvals[k] <- phyper(hit[k]-1, pop[k], M - pop[k], n, lower.tail = FALSE)
        }
        fdr <- p.adjust(pvals, method = "fdr")
        list(p = pvals, fdr = fdr, pop = pop, hit = hit)
      }
      
      ## === 4) Ejecuta para KEGG_Pathway y KEGG_Module (si existen) ===
      cats <- intersect(c("KEGG_Pathway", "KEGG_Module"), colnames(ann))
      ora_res <- lapply(cats, run_ora,
                        universe = universe_genes,
                        hits     = deg_genes,
                        ann_df   = ann)
      names(ora_res) <- cats
      
      ## === 5) Escribe resultados y genes por término significativo ===
      for (cat in names(ora_res)) {
        res <- ora_res[[cat]]
        if (!length(res$fdr)) next
        out_tab <- data.frame(term = names(res$fdr),
                              pval = unname(res$p[names(res$fdr)]),
                              FDR  = unname(res$fdr),
                              pop_count = unname(res$pop[names(res$fdr)]),
                              hit_count = unname(res$hit[names(res$fdr)]))
        out_tab <- out_tab[order(out_tab$FDR), ]
        write.table(out_tab,
                    file = sprintf("Enrichment_FDR_%s.tsv", cat),
                    sep  = "\t", quote = FALSE, row.names = FALSE)
        
        sig_terms <- out_tab$term[out_tab$FDR < 0.05]
        for (term in sig_terms) {
          gene_list <- rownames(ann)[grepl(term, ann[[cat]], fixed = TRUE)]
          gene_list <- intersect(gene_list, deg_genes)
          if (!length(gene_list)) next
          gene_out <- cbind(DEG_cpm[gene_list, , drop = FALSE],
                            ann[gene_list, grep("KEGG", colnames(ann), fixed = TRUE), drop = FALSE])
          write.table(gene_out,
                      file = sprintf("Enrichment_FDR_%s_%s.tsv", cat, term),
                      sep  = "\t", quote = FALSE)
        }
      }
**Enrichment_FDR_KEGG_Pathway**

Archivo del profe: 
<img width="933" height="681" alt="image" src="https://github.com/user-attachments/assets/0ba334f0-7149-47e4-9c58-720578dcf48e" />

Nuestro archivo: 
<img width="1176" height="976" alt="image" src="https://github.com/user-attachments/assets/17c36284-1d2f-46f9-91c4-41b0d1f11384" />


### Heatmaps of DEGs 0 vs 24

      ## Heatmaps de DEGs
      
      
      ## 1) Heatmap clásico de DEGs (genes × muestras)
      
      # ---- LIBRERÍAS ----
      library(pheatmap)
      library(RColorBrewer)
      library(matrixStats)
      
      # ---- ENTRADAS QUE YA TIENES ----
      # cs.log2  : matriz DEGs z-score (genes x muestras)
      # coldata  : data.frame con rownames = nombres de muestra y al menos la columna `time` (T0 / H24)
      
      stopifnot(identical(colnames(cs.log2), rownames(coldata)))
      
      # Paleta para condición/tiempo (ajústala si ya definiste una)
      pal_time <- c(T0 = "#79cbb8", H24 = "#500472")
      
      # Anotación de columnas (muestras)
      ann_col <- data.frame(time = coldata$time)
      rownames(ann_col) <- rownames(coldata)
      ann_colors <- list(time = pal_time)
      
      # (Opcional) si hay muchos DEGs, mostrar los más variables para que sea legible
      topN_1 <- min(2000, nrow(cs.log2))
      ord  <- order(rowVars(cs.log2), decreasing = TRUE)
      matH <- cs.log2[ord[seq_len(topN_1)], , drop = FALSE]
      
      # Heatmap
      pheatmap(
        matH,
        scale                  = "none",                             # ya está en z-score
        color                  = colorRampPalette(rev(brewer.pal(11, "RdBu")))(256),
        clustering_distance_rows = "correlation",
        clustering_distance_cols = "correlation",
        clustering_method      = "ward.D2",
        show_rownames          = FALSE,
        show_colnames          = TRUE,
        annotation_col         = ann_col,
        annotation_colors      = ann_colors,
        border_color           = NA,
        main                   = sprintf("DEGs T0 vs H24 — %d genes (z-score por gen)", nrow(matH))
      )


<img width="716" height="606" alt="image" src="https://github.com/user-attachments/assets/687922f8-eff7-4656-9206-cb3b50e5a5af" />

      
      ## 2) Heatmap de DEGs con anotaciones ricas (tiempo, severidad, variedad, tratamiento)
      
      
      # Si no tenemos nuestro archivo, bosquejo rápido:
      sample_info <- data.frame(
        time        = coldata$time,                                    # ya lo tenemos
        severidad   = factor(c("baja","baja","baja","baja","alta","alta","alta","baja",
                               "baja","baja","baja","baja","alta","alta","alta","baja"), 
                             levels = c("baja","alta")),
        variedad    = factor(c("catuai_mutante","catuai_mutante","catuai_mutante","catuai_mutante",
                               "catuai_mutante","catuai_mutante","catuai_wildtype","catuai_wildtype",
                               "catuai_mutante","catuai_mutante","catuai_mutante","catuai_mutante",
                               "catuai_mutante","catuai_mutante","catuai_wildtype","catuai_wildtype")),
        tratamiento = factor(c("EMS","azida de sodio","EMS","azida de sodio","EMS","azida de sodio","sin mutación","sin mutación",
                               "EMS","azida de sodio","EMS","azida de sodio","EMS","azida de sodio","sin mutación","sin mutación"))
      )
      rownames(sample_info) <- rownames(coldata)
      
      # Chequeo
      stopifnot(identical(rownames(sample_info), colnames(cs.log2)))
      
      # Paletas para cada anotación
      pal_time      <- c(T0 = "#79cbb8", H24 = "#500472")
      pal_severidad <- c(baja = "#88c057", alta = "#f06355")
      pal_variedad  <- c(catuai_mutante = "#6fa8dc", catuai_wildtype = "#ffd966")
      pal_trat      <- setNames(RColorBrewer::brewer.pal(6, "Set2")[1:nlevels(sample_info$tratamiento)],
                                levels(sample_info$tratamiento))
      
      ann_colors <- list(
        time        = pal_time,
        severidad   = pal_severidad,
        variedad    = pal_variedad,
        tratamiento = pal_trat
      )
      
      # Selección de genes (reutilizamos matH del Heatmap DEGs genes*samples)
      if (!exists("matH")) {
        topN <- min(2000, nrow(cs.log2))
        ord  <- order(matrixStats::rowVars(cs.log2), decreasing = TRUE)
        matH <- cs.log2[ord[seq_len(topN)], , drop = FALSE]
      }
      
      # Heatmap con múltiples anotaciones
      pheatmap(
        matH,
        scale                     = "none",
        color                     = colorRampPalette(rev(brewer.pal(11, "RdBu")))(256),
        clustering_distance_rows  = "correlation",
        clustering_distance_cols  = "correlation",
        clustering_method         = "ward.D2",
        show_rownames             = FALSE,
        show_colnames             = TRUE,
        annotation_col            = sample_info[, c("time","severidad","variedad","tratamiento"), drop = FALSE],
        annotation_colors         = ann_colors,
        border_color              = NA,
        main                      = "DEGs T0 vs H24 — con anotaciones"
      )

<img width="716" height="606" alt="image" src="https://github.com/user-attachments/assets/326e0f71-c91f-4468-8797-b38ee8955d77" />
      
      # Heatmap de enrequeciemiento
      
      library(pheatmap)
      library(dplyr)
      library(readr)
      
      # 1) Filtra/colapsa:
      kegg_clean <- kegg %>%
        mutate(core = sub("^(ko|map)", "", term)) %>%
        group_by(core) %>% slice_min(FDR, with_ties = FALSE) %>% ungroup() %>%
        mutate(score = -log10(FDR))
      if (!"term_name" %in% names(kegg_clean)) kegg_clean$term_name <- kegg_clean$term
      
      # 2) (Opcional) filtra categorías no vegetales (enfermedades humanas)
      drop_pat <- grepl("^map05", kegg_clean$term) |
        grepl("diabetes|cancer|cardio|Alzheimer|Parkinson|disease",
              kegg_clean$term_name, ignore.case = TRUE)
      kegg_clean <- kegg_clean[!drop_pat, ]
      
      # 3) Top N ordenado por FDR
      topN_2 <- 30
      #kegg_top <- kegg_clean %>% arrange(FDR) %>% slice_head(n = min(topN_2, n())) # arroja error
      kegg_top <- kegg_clean %>% arrange(FDR) %>% slice_head(n = min(topN_2, nrow(.)))
      # kegg_top <- head(arrange(kegg, FDR), topN) # No indica los nombres
      
      # 4) Matriz pathway × contraste
      mat_path <- matrix(kegg_top$score, ncol = 1)
      rownames(mat_path) <- kegg_top$term_name
      colnames(mat_path) <- "All_DEGs (T0 vs H24)"
      
      # 5) Heatmap (sin dendrograma de filas)
      pheatmap(
        mat_path,
        cluster_rows  = FALSE,   # <- sin dendrograma
        cluster_cols  = FALSE,
        color         = colorRampPalette(c("#ffffff","#ffe082","#f57c00","#b71c1c"))(256),
        show_rownames = TRUE,
        main          = "KEGG pathway enrichment (−log10 FDR)"
      )

<img width="716" height="606" alt="image" src="https://github.com/user-attachments/assets/5d32242e-ec84-4bf8-98b7-8caa459df3be" />

### Venn diagram of DEGs 0 vs 24


# Venn Diagrams -------

      library(VennDiagram)
      
      deg_up   <- rownames(res_lrt)[res_lrt$padj < 0.01 & res_lrt$log2FoldChange > 0]
      deg_down <- rownames(res_lrt)[res_lrt$padj < 0.01 & res_lrt$log2FoldChange < 0]
      
      venn.plot <- venn.diagram(
        x = list(
          "Upregulated (H24)"   = deg_up,
          "Downregulated (H24)" = deg_down
        ),
        filename = NULL,
        fill = c("#66c2a5", "#fc8d62"),
        alpha = 0.5,
        cex = 1.5,
        cat.cex = 1.5,
        cat.pos = 0
      )
      grid::grid.draw(venn.plot)

<img width="716" height="606" alt="image" src="https://github.com/user-attachments/assets/b4506a52-bba1-4394-b91a-c98a2a33b2be" />



## MAIN CONTRAST CODE
      
           
      # Coffea arabica challenge adressed by RNAseq under different types of mutagenesis at different levels of RUST infection
      
      # Packages ------
      
      library(kohonen);	      #This is the library for the SOM
      library(ggplot2);     	#This library is for transparency in the colors
      library(ggrepel)
      library(gplots);	      #Easy heatmaps también pheat más facile
      library(VennDiagram);	  #self explanatory
      library(pheatmap);	    #pretty heatmaps
      library(dendsort);    	#sorting dendrograms ayudar a rotar los dendogramas a
      library(DESeq2);	      #Normalization and everything related to that
      library(RColorBrewer)   #Color
      library(matrixStats)   # rowVars
      
      
      # INPUT -----
      url_rna_counts <- "https://raw.githubusercontent.com/jmvillalobos/RNAseq_curso2025/refs/heads/main/Cara_RNAc90_counts.txt"
      #url_rna_counts <- "Cara_RNAc90_counts.txt"
      
      # MAIN RAW DATA 
      rust_counts <- read.table(url_rna_counts, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
      
      head(rust_counts)
      
      
      # DATA TREATMENT -----------
      
      
      ### SHORT NAMES ----
      counts <- rust_counts
      
      ## 1) Build simple names 
      raw_names <- colnames(counts)
      simple    <- sub("^([TH])(\\d+).*", "\\1\\2", raw_names)
      
      ## 2) Apply and check
      colnames(counts) <- simple
      if (any(duplicated(colnames(counts)))) stop("Duplicate sample names after renaming")
      print(colnames(counts))
      head(counts)
      
      ## 3) Reorder from 1 to 16
      pref <- substr(colnames(counts), 1, 1)                     # "T" or "H"
      num  <- as.numeric(sub("^[TH]", "", colnames(counts)))     # numeric index
      ord  <- order(factor(pref, levels=c("T","H")), num)
      counts <- counts[, ord]
      print(colnames(counts))
      
      ## Ensure we have exactly the expected set name-info
      expected <- c(paste0("T", 1:8), paste0("H", 9:16))
      missing  <- setdiff(expected, colnames(counts))
      extra    <- setdiff(colnames(counts), expected)
      if (length(missing)) stop("Missing: ", paste(missing, collapse=", "))
      if (length(extra))   stop("Extra: ",   paste(extra, collapse=", "))
      
      
      
      
      ### Zero-count filter -----
      
      # Drop genes with zero counts across all samples
      counts <- counts[rowSums(counts) > 0, , drop = FALSE]
      cat("After zero-row filter:", nrow(counts), "genes and", ncol(counts), "samples\n")
      
      
      
      ### Grouping -------- 
      # grupos A/B/Control
      
      # Nos quedamos solo con H9–H15 (excluye H16)
      keep_samples <- c("H9","H10","H11","H12","H13","H14","H15")
      stopifnot(all(keep_samples %in% colnames(counts)))
      counts <- counts[, keep_samples, drop = FALSE]                 #Último movimiento a COUNTS
      
      
      # Etiquetado de grupos:
      #  A = alta susceptibilidad → H13, H14
      #  B = resistencia/baja severidad → H10, H12, H11, H9
      #  Control = H15 (Catuaí WT)
      grupo <- setNames(rep(NA_character_, length(keep_samples)), keep_samples)
      grupo[c("H13","H14")]                 <- "A"
      grupo[c("H10","H12","H11","H9")]      <- "B"
      grupo["H15"]                          <- "Control"
      
      coldata <- data.frame(
        grupo = factor(grupo, levels = c("A","B","Control")),
        row.names = names(grupo)
      )
      
      # chequeo rápido
      stopifnot(!any(is.na(coldata$grupo)))
      print(coldata)
      
      
      
      ### FPM screen ----
      
      # DESeq2 necesita counts enteros
      dds_raw <- DESeqDataSetFromMatrix(
        countData = as.matrix(round(counts)),
        colData   = coldata,
        design    = ~ grupo
      )
      
      # FPM-like (antes de normalizar por size factors)
      cpm_raw <- fpm(dds_raw, robust = FALSE)
      
      # Regla de prevalencia: FPM > 1 en al menos la mitad de réplicas del grupo
      grp     <- coldata$grupo
      mat_eval <- cpm_raw > 1
      # como nuestro Control = 1 muestra, exigir ≥ 1 ya es lo máximo posible. Si no, perdemos genes expresados solo ahí
      keep <- rep(FALSE, nrow(mat_eval))
      cutoff <- ceiling(table(grp) / 2)  # A=1, B=2, Control=1
            # Para A (2 muestras: H13, H14) → cutoff = 1 (al menos 1 de las 2 debe pasar).
            # Para B (4 muestras: H10, H12, H11, H9) → cutoff = 2 (al menos 2 de 4 deben pasar).
            # Para Control (1 muestra: H15) → cutoff = 1 (basta con que en esa única muestra FPM > 1).
      
      
      for (g in levels(grp)) {
        idx <- which(grp == g)
        keep <- keep | (rowSums(mat_eval[, idx, drop = FALSE]) >= cutoff[g])
      }
      
      cat("Genes kept by FPM > 1 prevalence rule:", sum(keep), "of", length(keep),
          sprintf(" (%.1f%%)\n", 100 * sum(keep) / length(keep)))
      
      # Subset
      counts_filt <- counts[keep, , drop = FALSE]
      dds_1 <- dds_raw[keep, ]
      
      # Quick sanity check
      stopifnot(identical(colnames(counts_filt), rownames(coldata)))
      
      
      
      ### Size-factor normalization (DESeq2) ------
      dds_1 <- estimateSizeFactors(dds_1)      # median-of-ratios
      sizeFactors(dds_1)
      sizeFactors(dds_1)[1:5]     # Quick peek
      
      
      
      ### VST →Variance Stabilizing Transform (VST) -------
      vst_1   <- vst(dds_1)
      vst_mat <- assay(vst_1)                  # genes x samples
      head(vst_mat)
      
      
      ### FPM genes*samples ---------  
      fpm_norm <- fpm(dds_1)   # genes x muestras
      head(fpm_norm)
      
      
      
      ### Gene-wise z-score ----
      zmat <- t(scale(t(vst_mat)))
      
      # Now our columns are H9 to H15
      stopifnot(all(is.finite(zmat)))
      cat("Mean of first 5 genes after z-score:\n"); print(rowMeans(zmat[1:5, ]))
      cat("SD of first 5 genes after z-score:\n"); print(apply(zmat[1:5, ], 1, sd))
      
      
      
      
      
      # PCA'S -----------------
      
      ### pre-PCA on samples (using prcomp) --------
      
      # We already z-scored by gene, so no extra scaling here.
      # zmat viene de VST + z-score (genes x muestras)
      pca <- prcomp(t(zmat), center = FALSE, scale. = FALSE)  # observations = samples
      var_expl <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
      
      
      
      ### Colors and legend (robust) ----
      pal_grupo <- c(A = "#e41a1c", B = "#377eb8", Control = "#4daf4a")
      col_grupo <- pal_grupo[as.character(coldata$grupo)]
      leg_labs  <- levels(coldata$grupo)
      
      
      
      ### PCA plots -----------
      
      ## 1) Pick the PCs you want to plot
      var_1 <- 1
      var_2 <- 2
      
      ## 2) Build a small data frame for ggplot
      df_pca <- data.frame(
        x     = pca$x[, var_1],
        y     = pca$x[, var_2],
        name  = rownames(pca$x),        # "H9".."H16"
        grupo = factor(coldata$grupo, levels = levels(coldata$grupo)) # ensure order
      )
      # pca$x → coordinates of each sample in PCA space
      # coldata$time → group each sample into condition/time ("T0"/"H24")
      
      ## 3) PLOT
      ggplot(df_pca, aes(x, y, color = grupo, label = name)) +
        geom_point(size = 2.6) +
        ggrepel::geom_text_repel(show.legend = FALSE, max.overlaps = Inf, size = 3) +
        scale_color_manual(values = pal_grupo, name = "Grupo") +
        labs(
          title = sprintf("PCA – Samples → PC%d vs PC%d", var_1, var_2),
          x = sprintf("PC%d (%.1f%%)", var_1, var_expl[var_1]),
          y = sprintf("PC%d (%.1f%%)", var_2, var_expl[var_2])
        ) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      
      ## 4) save PCA matrices
      # write.csv(pca$x, file = "PCA_samples_scores.csv")      # samples in PC space
      # write.csv(pca$rotation, file = "PCA_gene_loadings.csv") # genes loadings
      
      
      
      
      
      
      ## PCA gene plot  ------
      # Step 1: extract loadings (genes × PCs)
      loadings <- pca$rotation   # each row = a gene, each column = a PC
      
      
      # Step 3: build a dataframe for ggplot
      df_load <- data.frame(
        gene = rownames(loadings),
        PCx  = loadings[, var_1],
        PCy  = loadings[, var_2]
      )
      
      # Step 4: ggplot
      ggplot(df_load, aes(x = PCx, y = PCy, color = abs(PCx))) +
        geom_point(alpha = 0.4, size = 0.6, color = "#AF9AB2") +
        scale_color_viridis_c() +
        labs(
          title = sprintf("PCA – Genes ", var_1),
          x = sprintf("PC%d loadings (%.1f%%)", var_1, var_expl[var_1]),
          y = sprintf("PC%d loadings (%.1f%%)", var_2, var_expl[var_2]),
          color = sprintf("|loading PC%d|", var_1)
        ) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      
      
      # 5) COLORFULL PLOT PCA
      ggplot(df_load, aes(x = PCx, y = PCy, color = abs(PCx))) +
        geom_point(alpha = 0.6, size = 0.6) +
        scale_color_viridis_c() +
        labs(
          title = sprintf("PCA – Genes (|loading PC%d|)", var_1),
          x = sprintf("PC%d loadings (%.1f%%)", var_1, var_expl[var_1]),
          y = sprintf("PC%d loadings (%.1f%%)", var_2, var_expl[var_2]),
          color = sprintf("|loading PC%d|", var_1)
        ) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      
      
      
      
      
      
      
      
      
      
      
      
      ## PCA TOP genes ----
      
      # Inputs I already have
      # pca, var_expl, var_1, var_2
      gene_loadings <- pca$rotation  # safer name than 'loadings()'
      
      # Top K by sign (per side)
      K <- 100  # cambia a gusto
      ord_pos <- order(gene_loadings[, var_1], decreasing = TRUE)
      ord_neg <- order(gene_loadings[, var_1], decreasing = FALSE)
      pos_genes <- rownames(gene_loadings)[ord_pos][1:K]
      neg_genes <- rownames(gene_loadings)[ord_neg][1:K]
      
      # Build df with a 3-level flag
      df_load_top <- transform(df_load, group = "other")
      df_load_top$group[df_load_top$gene %in% pos_genes] <- "top_pos"
      df_load_top$group[df_load_top$gene %in% neg_genes] <- "top_neg"
      df_load_top$group <- factor(df_load_top$group, levels = c("other","top_pos","top_neg"))
      
      # Colors (hex ok)
      col_other <- "#AF9AB2"  
      col_pos   <- "#f9665e"  # red
      col_neg   <- "#799fcb"  # blue  
      
      
      ggplot(df_load_top, aes(PCx, PCy, color = group)) +
        geom_point(alpha = 0.65, size = 0.8) +
        scale_color_manual(
          values = c(other = col_other, top_pos = col_pos, top_neg = col_neg),
          labels = c(
            other   = "Other genes",
            top_pos = sprintf("Top %d positive on PC%d", K, var_1),
            top_neg = sprintf("Top %d negative on PC%d", K, var_1)
          ),
          name = "Highlight"
        ) +
        labs(
          title = sprintf("PCA – TOP Genes → PC%d vs PC%d", var_1, var_2),
          x = sprintf("PC%d loadings (%.1f%%)", var_1, var_expl[var_1]),
          y = sprintf("PC%d loadings (%.1f%%)", var_2, var_expl[var_2])
        ) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      
      # (Opcional) exportar listas
      write.csv(data.frame(gene = pos_genes,
                           loading = gene_loadings[pos_genes, var_1],
                           direction = "positive",
                           PC = paste0("PC", var_1)),
                sprintf("top_%d_positive_genes_PC%d.csv", K, var_1), row.names = FALSE)
      
      write.csv(data.frame(gene = neg_genes,
                           loading = gene_loadings[neg_genes, var_1],
                           direction = "negative",
                           PC = paste0("PC", var_1)),
                sprintf("top_%d_negative_genes_PC%d.csv", K, var_1), row.names = FALSE)
      
      
      
      
      
      
      
      
      #CLUSTERING & HEATMAPS -------
      
      ## 1) Fallback (only if vst_1/vst_mat don't exist) 
      if (!exists("vst_1")) {
        message("vst_1 not found, computing VST from dds_1...")
        stopifnot(exists("dds_1"))
        vst_1 <- vst(dds_1)
      }
      if (!exists("vst_mat")) vst_mat <- assay(vst_1)   # genes x samples
      # assay(vst1) → numeric matrix with rows = genes, cols = samples
      
      
      ### Heatmap sample sample ---------
      
      ## 2) SAMPLE–SAMPLE CLUSTERING → matriz de correlaciones entre muestras
      
      # Use Pearson correlation between samples (robust for RNA-seq after VST)
      cors <- cor(zmat, method = "pearson")     # samples x samples
      ann <- data.frame(Grupo = coldata$grupo)
      # Prepares sample annotations (here just time), to color the heatmap margins
      rownames(ann) <- rownames(coldata)
      
      
      ## 3) Heatmap of correlations sample -sample
      pheatmap(cors,
               annotation_col = ann,
               annotation_row = ann,
               color = colorRampPalette(c("blue", "white", "red"))(100),
               clustering_distance_rows = "correlation",
               clustering_distance_cols = "correlation",             #ASK WHY NOT EUCLIDEAN
               clustering_method = "average",         #ASK WHY NOT WARD.D2
               display_numbers  = FALSE,
               main = "Sample-to-sample correlation ")
      
      
      
      ## Heatmap top genes -------------
      
      # Select top N most variable genes across samples (after VST)
      rv <- rowVars(zmat)  # varianza por gen
      ntop <- 500          # número de genes a graficar
      sel_genes <- order(rv, decreasing = TRUE)[1:ntop]
      zmat_var <- zmat[sel_genes, ]
      
      # Heatmap (scale by row = z-score per gene)
      pheatmap(zmat_var,
               scale = "row",  # normaliza por gen (z-score)
               show_rownames = FALSE,
               show_colnames = TRUE,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "correlation",
               clustering_method = "ward.D2",
               annotation_col = ann,
               color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(255),
               main = sprintf("Top %d variable genes (unsupervised clustering)", ntop))
      
      ## SAVE TO FILES -
      # pdf("sample_correlation_heatmap.pdf", width = 7, height = 6); <repeat heatmap A>; dev.off()
      # pdf("top_variable_genes_heatmap.pdf", width = 8, height = 8); <repeat heatmap B>; dev.off()
      
      
      
      
      
      
      
      
      # DEGs & DOCUMENTS ----------
      
      ## 1) Counts se encuentra en ### Zero-count filter 
      
      
      ### EDIT 1) Parámetros del par y subconjunto de muestras
      
      AAA <- "A"   # numerador del log2FC (A vs B)
      BBB   <- "B"   # referencia
      CCC   <- "C"  # control
      #### AQUI SE CAMBIA ####
      contr_1 <- AAA
      contr_2 <- CCC
      
      
      samples_A <- c("H13","H14")
      samples_B <- c("H9","H10","H11","H12")
      samples_C <- c("H15")
      keep_samples_1 <- c(samples_A, samples_C) #### AQUI SE CAMBIA ####
      keep_samples_1
      
      # Subset columnas de acuerdo al contraste que queremos
      stopifnot(all(c(samples_B,  samples_A, samples_C) %in% colnames(counts))) # Corroboramos que nuestro documento tenga las mismas muestras que nuestros clusters
      counts_1 <- counts 
      colnames(counts) # Quick peek
      colnames(counts_1) # Quick peek
      stopifnot(all(keep_samples_1 %in% colnames(counts_1))) # Que sean las mismas
      counts_1 <- counts_1[, keep_samples_1, drop = FALSE]  # THE IMPORTATN COUNTS
      keep_samples_1 # Quick peek
      colnames(counts_1) # Quick peek
      
      
      # coldata_1 factor y nivel
      ### EDIT 2) colData con factor Group (A/B) y niveles (contr_1, contr_2)
      grp_vec <- setNames(rep(NA_character_, length(keep_samples_1)), keep_samples_1)
      #### AQUI SE CAMBIA ####
      grp_vec[samples_A] <- "A"       # Comenta el nivel no utilizado
      #grp_vec[samples_B] <- "B"
      grp_vec[samples_C] <- "C"
      
      grp_vec
      
      
      coldata_1 <- data.frame(Group = factor(grp_vec, levels = c(contr_1, contr_2)))
      rownames(coldata_1) <- names(grp_vec)
      
      # coldata (Dr. carlos) = coldata_1
      print(coldata_1)
      
      
      
      
      #### dds en contrastes --------
      
      ### EDIT 3) dds con diseño ~ Group (dos niveles)
      dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(counts_1)),
                                    colData = coldata_1,
                                    design = ~ Group)
      head(dds)
      
      
      # FPM-like (para referencia/consistencia con el flujo del profe)
      cpm <- fpm(dds, robust = FALSE)
      head(cpm)
      
      # VST “blind” 
      vst <- vst(dds, blind = TRUE)
      head(vst)
      
      
      
      # Matrices derivadas 
      std.mat <- t(scale(t(log(fpm(dds) + 1, 2))))  # matriz z-score desde log2(FPM+1)
      head(std.mat)
      mat.vst <- assay(vst)
      head(mat.vst)
      mat.cpm <- fpm(dds)
      head(mat.cpm)
      

### DEGs & DOCUMENTS

      # DEGs & DOCUMENTS ----------
      
      ## 1) Counts se encuentra en ### Zero-count filter 
      
      
      ### EDIT 1) Parámetros del par y subconjunto de muestras
      
      AAA <- "A"   # numerador del log2FC (A vs B)
      BBB   <- "B"   # referencia
      CCC   <- "C"  # control
      #### AQUI SE CAMBIA ####
      contr_1 <- CCC
      contr_2 <- AAA
      
      
      samples_A <- c("H13","H14")
      samples_B <- c("H9","H10","H11","H12")
      samples_C <- c("H15")
      keep_samples_1 <- c(samples_C, samples_A) #### AQUI SE CAMBIA ####
      keep_samples_1
      
      # Subset columnas de acuerdo al contraste que queremos
      stopifnot(all(c(samples_B,  samples_A, samples_C) %in% colnames(counts))) # Corroboramos que nuestro documento tenga las mismas muestras que nuestros clusters
      counts_1 <- counts 
      colnames(counts) # Quick peek
      colnames(counts_1) # Quick peek
      
      counts_1 <- counts_1[, keep_samples_1, drop = FALSE]  # THE IMPORTATN COUNTS
      stopifnot(all(keep_samples_1 %in% colnames(counts_1))) # Que sean las mismas
      keep_samples_1 # Quick peek
      colnames(counts_1) # Quick peek
      
      
      # coldata_1 factor y nivel
      ### EDIT 2) colData con factor Group (A/B) y niveles (contr_1, contr_2)
      grp_vec <- setNames(rep(NA_character_, length(keep_samples_1)), keep_samples_1)
      #### AQUI SE CAMBIA ####
      grp_vec[samples_A] <- "A"       # Comenta el nivel no utilizado
      #grp_vec[samples_B] <- "B"
      grp_vec[samples_C] <- "C"
      
      grp_vec
      
      
      coldata_1 <- data.frame(Group = factor(grp_vec, levels = c(contr_1, contr_2)))
      rownames(coldata_1) <- names(grp_vec)
      
      # coldata (Dr. carlos) = coldata_1
      print(coldata_1)
      
      
      
      
      #### dds en contrastes --------
      
      ### EDIT 3) dds con diseño ~ Group (dos niveles)
      dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(counts_1)),
                                    colData = coldata_1,
                                    design = ~ Group)
      head(dds)
      
      
      # FPM-like (para referencia/consistencia con el flujo del profe)
      cpm <- fpm(dds, robust = FALSE)
      head(cpm)
      
      # VST “blind” 
      vst <- vst(dds, blind = TRUE)
      head(vst)
      
      
      
      # Matrices derivadas 
      std.mat <- t(scale(t(log(fpm(dds) + 1, 2))))  # matriz z-score desde log2(FPM+1)
      head(std.mat)
      mat.vst <- assay(vst)
      head(mat.vst)
      mat.cpm <- fpm(dds)
      head(mat.cpm)
      
      
      
      ### LRT TEST ------
      
      ### EDIT 4) LRT (full ~ Group vs reduced ~ 1)
      
      # TEST LRT, comentar si WALD
      #dds <- DESeq(dds, test = "LRT", reduced = ~ 1)
      #res <- results(dds, alpha = 0.05)  # p-ajus y log2FC = contr_2 vs contr_1 (por niveles de Group)
      
      # TEST WALD, comentar si LRT
      dds <- DESeq(dds)  # Wald por defecto
      res <- results(dds, contrast = c("Group", contr_2, contr_1), alpha = 0.05)
      
      # --
        
      
      
      p.adj <- res$padj
      names(p.adj) <- rownames(mat.vst)
      p.adj[is.na(p.adj)] <- 1
      
      sum(p.adj < 0.05)
      
      # [1] n_degs (se imprimirá el conteo real)
      
      
      
      # Selección de DEGs como en el flujo original (desde mat.cpm)
      DEG.cpm <- mat.cpm[(p.adj < 0.05), , drop = FALSE]
      
      log2  <- log(DEG.cpm + 1, 2)
      cs.log2 <- t(scale(t(log2)))   # z-score por gen (para heatmaps/SOM, etc.)
      
      
      
      ### Print res as a table -----
      res_df    <- as.data.frame(res)
      head(res_df)
      head(res)
      
      # contrastes están en contr_1 y contr_2
      out_file <- paste0("DEG_gene_results_", contr_1, "_vs_", contr_2, ".csv")
      
      write.csv(res_df, file = out_file, row.names = TRUE )
      
      
      getwd()
      
      ### Enriquecimiento  --------
      annotation <- read.delim("fullAnnotation.tsv.txt", stringsAsFactors = FALSE, row.names = 1)
      
      enriched.FDR <- enriched.p <- vector("list", 2)
      names(enriched.p) <- names(enriched.FDR) <- c("KEGG_Pathway", "KEGG_Module")
      
      for (i in 1:2) {
        population <- table(unlist(strsplit(annotation[rownames(annotation) %in% rownames(std.mat),
                                                       names(enriched.p)[[i]]], split = ",")))
        results.p <- array(1, length(population)); names(results.p) <- names(population)
        results.FDR <- results.p
        hit <- population * 0
        
        aux <- table(unlist(strsplit(annotation[rownames(annotation) %in% rownames(DEG.cpm),
                                                names(enriched.p)[[i]]], split = ",")))
        hit[names(aux)] <- aux
        
        for (k in names(aux)) {
          # Over-representation: phyper(hitInSample-1, hitInPopulation, failInPopulation, sampleSize, lower.tail=FALSE)
          results.p[k] <- phyper(hit[k] - 1, population[k], sum(population) - population[k], sum(hit), lower.tail = FALSE)
        }
        
        results.FDR <- p.adjust(results.p, method = "fdr")
        enriched.p[[i]] <- results.p
        enriched.FDR[[i]] <- results.FDR
        
        write.table(t(t(results.FDR)),
                    file = paste0("Enrichment_FDR_", names(enriched.p)[[i]], ".tsv"),
                    sep = "\t")
        
        for (j in names(results.FDR)[results.FDR < 0.05]) {
          enriched.current <- intersect(
            rownames(annotation)[grep(j, annotation[, names(enriched.p)[[i]]])],
            rownames(DEG.cpm)
          )
          write.table(
            cbind(DEG.cpm[enriched.current, ],
                  annotation[enriched.current, grep("KEGG", colnames(annotation))]),
            file = paste0("Enrichment_FDR_", names(enriched.p)[[i]], "_", j, ".tsv"),
            sep  = "\t"
          )
        }
      }
      
      results.FDR[[1]]
      
      
      
      
      # DEGs for pattern discovery -----
      
      ## Transform DEGs for pattern discovery  ----
      contr_1
      contr_2
      tag <- paste0(contr_2, "_vs_", contr_1)  # E.j. A_vs_C
      
      # 0) máscara (contraste) de DEGs por FDR
      alpha <- 0.05
      padj <- res$padj
      padj[is.na(padj)] <- 1              # <- reemplaza NAs por 1 (no significativo)
      deg_mask <- padj < alpha
      
      cat("N NAs en res$padj:", sum(is.na(res$padj)), "\n")
      cat("DEGs @FDR<", alpha, ": ", sum(deg_mask), "\n", sep="")
      
      if (!any(deg_mask)) stop("No DEGs @ FDR<", alpha, " en ", tag)
      
      
      # 1) FPM normalizados con los size factors de ESTE dds
      fpm_norm <- fpm(dds)                      # genes x muestras
      DEG_cpm  <- fpm_norm[deg_mask, , drop=FALSE]
      cat("[", tag, "] DEG_cpm:", nrow(DEG_cpm), "genes x", ncol(DEG_cpm), "muestras\n")
      
      # 2) log2(+1)
      log2_deg <- log2(DEG_cpm + 1)
      
      # 3) z-score por gen
      cs.log2 <- t(scale(t(log2_deg)))
      
      # 4) checks
      stopifnot(all(is.finite(cs.log2)))
      stopifnot(identical(colnames(cs.log2), rownames(coldata_1)))
      cat("z-score range: [", round(min(cs.log2),2), ", ", round(max(cs.log2),2), "]\n", sep="")
      
      # 5) guardar (opcional)
      # write.csv(as.data.frame(DEG_cpm), file=paste0("DEGs_FPMnorm_", tag, ".csv"))
      # write.csv(as.data.frame(cs.log2), file=paste0("DEGs_log2p1_zscore_", tag, ".csv"))



Analicemos pathways
Los csv no analicemos module archives "es priorizar"

<img width="1919" height="1199" alt="image" src="https://github.com/user-attachments/assets/026de977-6b07-42a5-851b-ba6f25eba874" />

eggnote es estadistica aplicada → Cualquier persona que use estadistica aplicada le puede llamar AI


####################################3
#################################

# Evidencia de que mi counts es el mismo counts del profesor
<img width="1217" height="799" alt="image" src="https://github.com/user-attachments/assets/922679a8-4851-448b-8d06-7a8143b08b11" />



Volcano plot con los padj 



Diagrama de venn




> [!NOTE]
> Useful information that users should know, even when skimming content.

> [!TIP]
> Helpful advice for doing things better or more easily.

> [!IMPORTANT]
> Key information users need to know to achieve their goal.

> [!WARNING]
> Urgent info that needs immediate user attention to avoid problems.

> [!CAUTION]
> Advises about risks or negative outcomes of certain actions.
