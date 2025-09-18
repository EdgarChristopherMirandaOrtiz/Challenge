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
      library(dplyr)
      library(tibble)
      library(stringr)
      library(readr)
      library(ggvenn)
      
      
      
      
      # INPUT -----
      #url_rna_counts <- "https://raw.githubusercontent.com/jmvillalobos/RNAseq_curso2025/refs/heads/main/Cara_RNAc90_counts.txt"
      url_rna_counts <- "Cara_RNAc90_counts.txt"
      
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
      
      counts_good <- counts
      
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
      pal_grupo <- c(A = "#fdac4f", B = "#a2593d", C = "#3d6863")
      pall_grupo <- c(A = "#fdac4f", B = "#a2593d", Control = "#3d6863")
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
        geom_point(alpha = 0.4, size = 0.6, color = "gray80") +
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
      col_other <- "gray80"  
      col_pos   <- "#75621d"  # red
      col_neg   <- "#76244c"  # blue  
      
      
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
               color = colorRampPalette(c("#76244c", "#f7efe6", "#75621d"))(100),
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
      
      
      
      
      # --DEGs & DOCUMENTS ------------------------------------------------
      
      # ENRIQ... log2FC (B/A) -------------------------------------------
      # Corre "annotation"
      
      #  Parámetros y subconjunto 
      AAA <- "A"   
      BBB   <- "B"   
      CCC   <- "C"  # control
      contr_1a <- AAA       # referencia (B)
      contr_2a <- BBB       # numerador (A)
      samples_A <- c("H13","H14")
      samples_B <- c("H9","H10","H11","H12")
      samples_C <- c("H15")
      keep_samples_a <- c(samples_A, samples_B)
      stopifnot(all(keep_samples_a %in% colnames(counts)))
      
      
      counts_a <- counts[, keep_samples_a, drop = FALSE]
      
      #  colData (niveles en orden ref→num) 
      grp_vec_a <- setNames(rep(NA_character_, length(keep_samples_a)), keep_samples_a)
      grp_vec_a[samples_A] <- "A"
      grp_vec_a[samples_B] <- "B"
      
      coldata_a <- data.frame(Group = factor(grp_vec_a, levels = c(contr_1a, contr_2a)))
      rownames(coldata_a) <- names(grp_vec_a)
      
      ### DESeq2: WALD -------
      dds_a <- DESeqDataSetFromMatrix(countData = as.matrix(round(counts_a)),
                                      colData   = coldata_a,
                                      design    = ~ Group)
      dds_a <- dds_a[rowSums(counts(dds_a)) > 0, ]
      dds_a <- DESeq(dds_a)
      res_a <- results(dds_a, contrast = c("Group", contr_2a, contr_1a), alpha = 0.05)
      resdf_a <- as.data.frame(res_a)
      
      # matrices derivadas (para enriquecimiento / heatmaps) 
      fpm_a     <- fpm(dds_a)
      vst_a     <- vst(dds_a, blind = TRUE)
      std.mat_a <- t(scale(t(log2(fpm(dds_a) + 1))))   # z-score desde log2(FPM+1)
      
      # --- DEGs por padj 
      padj_a <- res_a$padj; names(padj_a) <- rownames(res_a); padj_a[is.na(padj_a)] <- 1
      DEG.cpm_a <- fpm_a[padj_a < 0.05, , drop = FALSE]
      
      
      ### Exporta tabla ----------
      res_sig_a <- subset(resdf_a, !is.na(padj) & padj < 0.05)
      write.csv(res_sig_a, file = paste0("DEG_gene_results_", contr_2a, "_vs_", contr_1a, ".csv"),
                row.names = TRUE)
      
      ### enriquecimiento ---------
      # ENRIQUECIMIENTO (KEGG_Pathway / KEGG_Module), 
      annotation <- read.delim("fullAnnotation.tsv.txt", stringsAsFactors = FALSE, row.names = 1)
      
      
      enriched.FDR_a <- enriched.p_a <- vector("list", 2)
      names(enriched.p_a) <- names(enriched.FDR_a) <- c("KEGG_Pathway", "KEGG_Module")
      
      for (i in 1:2) {
        fam <- names(enriched.p_a)[[i]]
        
        population <- table(unlist(strsplit(annotation[rownames(annotation) %in% rownames(std.mat_a),
                                                       fam], split = ",")))
        results.p  <- array(1, length(population)); names(results.p) <- names(population)
        results.FDR<- results.p
        hit        <- population * 0
        
        aux <- table(unlist(strsplit(annotation[rownames(annotation) %in% rownames(DEG.cpm_a),
                                                fam], split = ",")))
        hit[names(aux)] <- aux
        
        for (k in names(aux)) {
          results.p[k] <- phyper(hit[k]-1, population[k], sum(population)-population[k],
                                 sum(hit), lower.tail = FALSE)
        }
        
        results.FDR         <- p.adjust(results.p, method = "fdr")
        enriched.p_a[[i]]   <- results.p
        enriched.FDR_a[[i]] <- results.FDR
        
        # archivo resumen por familia (Pathway/Module) con sufijo a
        write.table(t(t(results.FDR)),
                    file = paste0("Enrichment_FDR_", fam, "_", contr_2a, "_vs_", contr_1a, ".tsv"),
                    sep  = "\t")
        # export por cada término significativo (sus genes DEG + anotación KEGG)
        for (j in names(results.FDR)[results.FDR < 0.05]) {
          enriched.current <- intersect(
            rownames(annotation)[grep(j, annotation[, fam])],
            rownames(DEG.cpm_a)
          )
          if (length(enriched.current)) {
            write.table(
              cbind(DEG.cpm_a[enriched.current, ],
                    annotation[enriched.current, grep("KEGG", colnames(annotation))]),
              file = paste0("Enrichment_FDR_", fam, "_", j, "_", contr_2a, "_vs_", contr_1a, ".tsv"),
              sep  = "\t"
            )
          }
        }
      }
      
      
      # Corre "annotation"
      # ENRIQ... log2FC (A/C) -------------------------------------------
      
      # Parámetros y subconjunto 
      contr_1b <- CCC   # referencia = A
      contr_2b <- AAA   # numerador = C
      keep_samples_b <- c(samples_C, samples_A)
      stopifnot(all(keep_samples_b %in% colnames(counts)))
      
      counts_b <- counts[, keep_samples_b, drop = FALSE]
      
      # colData (niveles en orden ref→num) 
      grp_vec_b <- setNames(rep(NA_character_, length(keep_samples_b)), keep_samples_b)
      grp_vec_b[samples_C] <- "C"
      grp_vec_b[samples_A] <- "A"
      
      coldata_b <- data.frame(Group = factor(grp_vec_b, levels = c(contr_1b, contr_2b)))
      rownames(coldata_b) <- names(grp_vec_b)
      
      ### DESeq2: WALD ---------
      dds_b <- DESeqDataSetFromMatrix(countData = as.matrix(round(counts_b)),
                                      colData   = coldata_b,
                                      design    = ~ Group)
      dds_b <- dds_b[rowSums(counts(dds_b)) > 0, ]
      dds_b <- DESeq(dds_b)
      res_b  <- results(dds_b, contrast = c("Group", contr_2b, contr_1b), alpha = 0.05)
      resdf_b <- as.data.frame(res_b)
      
      #  matrices derivadas (para enriquecimiento / heatmaps) 
      fpm_b     <- fpm(dds_b)
      vst_b     <- vst(dds_b, blind = TRUE)
      std.mat_b <- t(scale(t(log2(fpm(dds_b) + 1))))   # z-score desde log2(FPM+1)
      
      # DEGs por padj 
      padj_b <- res_b$padj; names(padj_b) <- rownames(res_b); padj_b[is.na(padj_b)] <- 1
      DEG.cpm_b <- fpm_b[padj_b < 0.05, , drop = FALSE]
      
      ### Exporta tabla --------
      res_sig_b <- subset(resdf_b, !is.na(padj) & padj < 0.05)
      write.csv(res_sig_b,
                file = paste0("DEG_gene_results_", contr_2b, "_vs_", contr_1b, ".csv"),
                row.names = TRUE)
      
      ### enriquecimiento -------------
      # ENRIQUECIMIENTO (KEGG_Pathway / KEGG_Module) 
      if (!exists("annotation")) {
        annotation <- read.delim("fullAnnotation.tsv.txt", stringsAsFactors = FALSE, row.names = 1)
      }
      
      enriched.FDR_b <- enriched.p_b <- vector("list", 2)
      names(enriched.p_b) <- names(enriched.FDR_b) <- c("KEGG_Pathway", "KEGG_Module")
      
      for (i in 1:2) {
        fam <- names(enriched.p_b)[[i]]
        
        population <- table(unlist(strsplit(
          annotation[rownames(annotation) %in% rownames(std.mat_b), fam], split = ","
        )))
        results.p   <- array(1, length(population)); names(results.p) <- names(population)
        results.FDR <- results.p
        hit         <- population * 0
        
        aux <- table(unlist(strsplit(
          annotation[rownames(annotation) %in% rownames(DEG.cpm_b), fam], split = ","
        )))
        hit[names(aux)] <- aux
        
        for (k in names(aux)) {
          results.p[k] <- phyper(hit[k]-1, population[k], sum(population)-population[k],
                                 sum(hit), lower.tail = FALSE)
        }
        
        results.FDR          <- p.adjust(results.p, method = "fdr")
        enriched.p_b[[i]]    <- results.p
        enriched.FDR_b[[i]]  <- results.FDR
        
        write.table(t(t(results.FDR)),
                    file = paste0("Enrichment_FDR_", fam, "_", contr_2b, "_vs_", contr_1b, ".tsv"),
                    sep = "\t")
        
        for (j in names(results.FDR)[results.FDR < 0.05]) {
          enriched.current <- intersect(
            rownames(annotation)[grep(j, annotation[, fam])],
            rownames(DEG.cpm_b)
          )
          if (length(enriched.current)) {
            write.table(
              cbind(DEG.cpm_b[enriched.current, ],
                    annotation[enriched.current, grep("KEGG", colnames(annotation))]),
              file = paste0("Enrichment_FDR_", fam, "_", j, "_", contr_2b, "_vs_", contr_1b, ".tsv"),
              sep  = "\t"
            )
          }
        }
      }
      
      
      # Corre "annotation"
      # ENRIQ... log2FC (B/C) -------------------------------------------
      
      #  Parámetros y subconjunto 
      contr_1c <- CCC   # referencia = C
      contr_2c <- BBB   # numerador = B
      keep_samples_c <- c(samples_C, samples_B)
      stopifnot(all(keep_samples_c %in% colnames(counts)))
      
      counts_c <- counts[, keep_samples_c, drop = FALSE]
      
      # colData (niveles en orden ref→num) 
      grp_vec_c <- setNames(rep(NA_character_, length(keep_samples_c)), keep_samples_c)
      grp_vec_c[samples_C] <- "C"
      grp_vec_c[samples_B] <- "B"
      
      coldata_c <- data.frame(Group = factor(grp_vec_c, levels = c(contr_1c, contr_2c)))
      rownames(coldata_c) <- names(grp_vec_c)
      
      ### DESeq2: WALD ---------------------
      dds_c <- DESeqDataSetFromMatrix(countData = as.matrix(round(counts_c)),
                                      colData   = coldata_c,
                                      design    = ~ Group)
      dds_c <- dds_c[rowSums(counts(dds_c)) > 0, ]
      dds_c <- DESeq(dds_c)
      res_c  <- results(dds_c, contrast = c("Group", contr_2c, contr_1c), alpha = 0.05)
      resdf_c <- as.data.frame(res_c)
      
      # matrices derivadas (para enriquecimiento / heatmaps) 
      fpm_c     <- fpm(dds_c)
      vst_c     <- vst(dds_c, blind = TRUE)
      std.mat_c <- t(scale(t(log2(fpm(dds_c) + 1))))   # z-score desde log2(FPM+1)
      
      # DEGs por padj 
      padj_c <- res_c$padj; names(padj_c) <- rownames(res_c); padj_c[is.na(padj_c)] <- 1
      DEG.cpm_c <- fpm_c[padj_c < 0.05, , drop = FALSE]
      
      ### Exporta tabla ---------
      res_sig_c <- subset(resdf_c, !is.na(padj) & padj < 0.05)
      write.csv(res_sig_c,
                file = paste0("DEG_gene_results_", contr_2c, "_vs_", contr_1c, ".csv"),
                row.names = TRUE)
      
      ### enriquecimiento -------
      # (KEGG_Pathway / KEGG_Module) 
      if (!exists("annotation")) {
        annotation <- read.delim("fullAnnotation.tsv.txt", stringsAsFactors = FALSE, row.names = 1)
      }
      
      enriched.FDR_c <- enriched.p_c <- vector("list", 2)
      names(enriched.p_c) <- names(enriched.FDR_c) <- c("KEGG_Pathway", "KEGG_Module")
      
      for (i in 1:2) {
        fam <- names(enriched.p_c)[[i]]
        
        population <- table(unlist(strsplit(
          annotation[rownames(annotation) %in% rownames(std.mat_c), fam], split = ","
        )))
        results.p   <- array(1, length(population)); names(results.p) <- names(population)
        results.FDR <- results.p
        hit         <- population * 0
        
        aux <- table(unlist(strsplit(
          annotation[rownames(annotation) %in% rownames(DEG.cpm_c), fam], split = ","
        )))
        hit[names(aux)] <- aux
        
        for (k in names(aux)) {
          results.p[k] <- phyper(hit[k]-1, population[k], sum(population)-population[k],
                                 sum(hit), lower.tail = FALSE)
        }
        
        results.FDR          <- p.adjust(results.p, method = "fdr")
        enriched.p_c[[i]]    <- results.p
        enriched.FDR_c[[i]]  <- results.FDR
        
        write.table(t(t(results.FDR)),
                    file = paste0("Enrichment_FDR_", fam, "_", contr_2c, "_vs_", contr_1c, ".tsv"),
                    sep = "\t")
        
        for (j in names(results.FDR)[results.FDR < 0.05]) {
          enriched.current <- intersect(
            rownames(annotation)[grep(j, annotation[, fam])],
            rownames(DEG.cpm_c)
          )
          if (length(enriched.current)) {
            write.table(
              cbind(DEG.cpm_c[enriched.current, ],
                    annotation[enriched.current, grep("KEGG", colnames(annotation))]),
              file = paste0("Enrichment_FDR_", fam, "_", j, "_", contr_2c, "_vs_", contr_1c, ".tsv"),
              sep  = "\t"
            )
          }
        }
      }
      
      
      
      
      
      
      
      # PLOTING ---------
      
      ### Preparativos (correr una sola vez) -----
      
      
      stopifnot(exists("counts"), is.matrix(counts) || is.data.frame(counts))
      counts <- as.matrix(counts)
      
      # Define nombres de grupos (solo etiquetas, para los niveles del factor)
      AAA <- "A"   # numerador del log2FC (A vs B)      
      BBB <- "B"   # referencia
      CCC <- "C"   # control
      
      # Define los vectores de muestras por grupo (reutilizados en los 3 contrastes)
      samples_A <- c("H13","H14")
      samples_B <- c("H9","H10","H11","H12")
      samples_C <- c("H15")
      
      stopifnot(all(c(samples_A, samples_B, samples_C) %in% colnames(counts_good)))
      
      # Umbral para volcano y heatmap (ajústalo si quieres)
      alpha_padj <- 0.05
      lfc_thr    <- 1
      n_top_hm   <- 30
      
      
      
      
      
      
      
      
      
      ### Contraste 1: A vs B  (log2FC = A/B) ----------
      
      # log2FC = contr_2 / contr_1
      
      # log2FC > 0 ⇒ el gen está sobreexpresado (upregulated) en contr_2 respecto a contr_1
      # log2FC < 0 ⇒ el gen está subexpresado (downregulated) en contr_2 respecto a contr_1
      # Puedes decir “sobre/subrepresentado” si hablas de abundancia/CPM, pero para DEGs es más claro up/downregulated
      
      
      # 1) Parámetros del par y subconjunto de muestras
      contr_1a <- AAA   # referencia (denominador)
      contr_2a <- BBB   # numerador
      keep_samples_a <- c(samples_B, samples_A)  # B primero, A después
      stopifnot(all(keep_samples_a %in% colnames(counts_good)))
      
      counts_1a <- counts_good[, keep_samples_a, drop = FALSE]
      
      # 2) colData con factor Group y niveles (contr_1a, contr_2a)
      grp_vec_a <- setNames(rep(NA_character_, length(keep_samples_a)), keep_samples_a)
      grp_vec_a[samples_A] <- "A"
      grp_vec_a[samples_B] <- "B"
      
      coldata_a <- data.frame(Group = factor(grp_vec_a, levels = c(contr_1a, contr_2a)))
      rownames(coldata_a) <- names(grp_vec_a)
      
      # 2.5) Chequeos y redondeo a enteros (recomendado para DESeq2)
      stopifnot(identical(colnames(counts_1a), rownames(coldata_a)))     # mismas muestras
      if (any(is.na(counts_1a))) stop("Hay NAs en counts_1a")
      
      bad_a <- which(abs(counts_1a - round(counts_1a)) > 1e-8, arr.ind = TRUE)
      cat("Celdas no enteras en counts_1a:", nrow(bad_a), "\n")
      counts_1a_int <- as.matrix(round(counts_1a))
      
      # 3) DESeq2 (sin shrink)
      dds_a <- DESeqDataSetFromMatrix(countData = counts_1a_int, colData = coldata_a, design = ~ Group)
      dds_a <- dds_a[rowSums(counts(dds_a)) > 0, ]   # por si queda algún gen en cero
      dds_a <- DESeq(dds_a)
      res_a <- results(dds_a, contrast = c("Group", contr_2a, contr_1a))  # A vs B
      resdf_a <- as.data.frame(res_a)
      resdf_a$gene <- rownames(resdf_a)
      
      # Etiquetas para barplot (dirección y significancia)
      resdf_a$direction <- ifelse(resdf_a$log2FoldChange > 0, paste0("up in ", contr_2a),
                                  ifelse(resdf_a$log2FoldChange < 0, paste0("down in ", contr_2a), "no change"))
      resdf_a$sig <- ifelse(!is.na(resdf_a$padj) & resdf_a$padj < alpha_padj, "significant", "ns")
      
      # 4) Normalización para gráficos (VST)
      vst_a <- vst(dds_a, blind = TRUE)
      normcounts_a <- assay(vst_a)  # para PCA/Heatmap
      
      #### PCA (A vs B) --------------------------------------------------------------
      # Separación global por grupos; puntos = muestras
      # PCA (Contraste A vs B) — usa vst_a y sufijos 'a'
      pca_a <- DESeq2::plotPCA(vst_a, intgroup = "Group", returnData = TRUE)
      
      # solo niveles del contraste, en orden: B (ref) y A
      pca_a$Group <- factor(as.character(pca_a$Group), levels = c(contr_1a, contr_2a))
      
      # colores fijos para B y A
      col_ab <- setNames(c("#fdac4f", "#a2593d"), c(contr_1a, contr_2a))
      
      ggplot(pca_a, aes(PC1, PC2, color = Group, label = name)) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(show.legend = FALSE, size = 3) +
        scale_color_manual(values = col_ab, breaks = c(contr_1a, contr_2a), name = "Group") +
        labs(
          title    = paste0("PCA (", contr_2a, " vs ", contr_1a, ")"),
          subtitle = "Transformación VST",
          x = paste0("PC1: ", round(100 * attr(pca_a, "percentVar")[1], 1), "%"),
          y = paste0("PC2: ", round(100 * attr(pca_a, "percentVar")[2], 1), "%")
        ) +
        theme_classic(base_size = 10) +
        theme(
          plot.title    = element_text(hjust = 0.5, face = "bold"),  # centrado
          plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30"),  # a la derecha
          axis.title    = element_text(size = 12),
          axis.text     = element_text(size = 12),
          legend.title  = element_text(size = 10),
          legend.text   = element_text(size = 10)
        )
      
      
      
      #### MA-plot (A vs B) ----------------------------------------------------------
      # log2FC vs abundancia media; DESeq2 resalta significativos
      # --- preparar datos MA (A = log10 baseMean, M = log2FC) ---
      ma_a <- transform(as.data.frame(res_a),
                        A = log10(baseMean + 1),
                        M = log2FoldChange)
      
      # reglas de significancia y dirección (reusa tus umbrales)
      if (!exists("alpha_padj")) alpha_padj <- 0.05
      if (!exists("lfc_thr"))    lfc_thr    <- 1
      
      ma_a$flag <- with(ma_a,
                        ifelse(!is.na(padj) & padj < alpha_padj & abs(M) >= lfc_thr,
                               ifelse(M > 0, paste0("up in ", contr_2a),
                                      paste0("down in ", contr_2a)),
                               "ns"))
      
      # paleta consistente con el resto
      cols_a <- c("ns" = "grey80",
                  setNames("#a5923d", paste0("up in ",   contr_2a)),
                  setNames("#aa4d6f", paste0("down in ", contr_2a)))
      
      # Límites que quieras (ajústalos a gusto)
      x_lim_a <- c(0, 6)      # eje X en log10(baseMean+1)
      y_lim_a <- c(-12, 12)   # eje Y en log2FC
      
      # --- GGPlot con el “look” del PCA ---
      p_ma_a <- ggplot(ma_a, aes(A, M)) +
        geom_point(aes(color = flag), size = 1.2, alpha = 0.7) +
        geom_hline(yintercept = 0, colour = "grey40") +
        geom_hline(yintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
        scale_color_manual(values = cols_a, name = "") +
        labs(
          title    = sprintf("MA-plot (%s vs %s)", contr_2a, contr_1a),
          subtitle = sprintf("DESeq2 • (padj < %.2f) & (log2FC ≥ %g)", alpha_padj, lfc_thr),
          x = "log10(baseMean + 1)",
          y = "log2 Fold Change"
        ) +
        coord_cartesian(xlim = x_lim_a, ylim = y_lim_a) +
        theme_classic(base_size = 12) +
        theme(
          plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30"),
          axis.title    = element_text(size = 12),
          axis.text     = element_text(size = 12),
          legend.title  = element_text(size = 10),
          legend.text   = element_text(size = 10)
        )
      
      p_ma_a
      
      
      #### Volcano (A vs B) ----------------------------------------------------------
      
      # 1) clasificar puntos (up/down/ns) y evitar Inf cuando padj = 0
      volcano_a <- transform(
        resdf_a,
        padj_plot = pmin(padj, 1),
        dir = ifelse(is.na(padj) | padj >= alpha_padj | abs(log2FoldChange) < lfc_thr,
                     "ns",
                     ifelse(log2FoldChange > 0,
                            paste0("up in ",   contr_2a),
                            paste0("down in ", contr_2a)))
      )
      
      # 2) paleta como en MA-plot
      cols_a <- setNames(  c("#aa4d6f", "grey75", "#a5923d"),
        c(paste0("down in ", contr_2a), "ns", paste0("up in ", contr_2a))
      )
      
      
      # 3) gráfico
      library(ggplot2)
      ggplot(volcano_a, aes(x = log2FoldChange, y = -log10(padj_plot), color = dir)) +
        geom_point(alpha = 0.75, size = 1.2) +
        geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", colour = "grey40") +
        geom_hline(yintercept = -log10(alpha_padj),  linetype = "dashed", colour = "grey40") +
        scale_color_manual(values = cols_a, name = NULL) +
        labs(
          title    = "Volcano (A vs B)",
          subtitle = "DESeq2 • (padj < 0.05) & (log2FC ≥ 1)",
          x = paste0("log2FC (", contr_2a, " / ", contr_1a, ")"),
          y = "-log10(padj)"
        ) +
        theme_classic(base_size = 12) +
        theme(
          plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30"),
          axis.title    = element_text(size = 12),
          axis.text     = element_text(size = 12),
          legend.title  = element_text(size = 10),
          legend.text   = element_text(size = 10)
        )
      
      
      
      
      #### Heatmap (A vs B) ----------------------------------
      # top genes por padj (A vs B)
      
      # 1) Selección de genes más significativos por padj
      n_top_hm <- 150                                    # cambiar
      ord_a    <- order(resdf_a$padj, na.last = NA)
      top_a    <- head(rownames(resdf_a)[ord_a], n_top_hm)
      
      # 2) Matriz VST normalizada + centrado por gen (z-score por fila sin escalar sd)
      mat_a0 <- normcounts_a[top_a, , drop = FALSE]
      mat_a  <- mat_a0 - rowMeans(mat_a0)
      
      # 3) Anotaciones de columnas (Group) con niveles en orden (contr_1a = referencia, contr_2a = numerador)
      ann_a <- data.frame(Group = coldata_a$Group,
                          row.names = rownames(coldata_a))
      ann_a$Group <- factor(as.character(ann_a$Group),
                            levels = c(contr_1a, contr_2a))
      
      # 1) Colores extremos (ajústalos si quieres)
      col_neg_a <- "#aa4d6f"  # morado (menor que la media de su gen)
      col_pos_a <- "#a5923d"  # verde (mayor que la media de su gen)
      
      # Color group
      ann_colors_a <- list(
        Group = setNames(
          c("#fdac4f", "#a2593d"),            # B = café, A = amarillo (ajusta si quieres)
          c(contr_1a,  contr_2a)
        )
      )
      
      # 2) Degradado azul → blanco → rojo (blanco en 0 para que el centro sea neutro)
      nbreaks_a <- 256
      lim_a     <- max(abs(mat_a))  # simétrico alrededor de 0
      breaks_a  <- seq(-lim_a, lim_a, length.out = nbreaks_a)
      hm_cols_a <- colorRampPalette(c(col_neg_a, "#f7efe6", col_pos_a))(nbreaks_a - 1)
      
      # 3) Ticks de la barra de color (leyenda de la escala)
      ticks_a <- pretty(c(-lim_a, lim_a), n = 5)
      
      # 4) Heatmap
      pheatmap(mat_a,
               color = hm_cols_a, breaks = breaks_a,
               scale = "none",                         # ya centramos por fila
               border_color = NA,
               show_rownames = FALSE,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "correlation",
               clustering_method = "ward.D2",
               annotation_col    = ann_a,
               annotation_colors = ann_colors_a,
               legend_breaks = ticks_a,
               legend_labels = sprintf("%.1f", ticks_a),
               main = sprintf("Heatmap top %d en %s vs %s\n",             # \n → Sirve para dar enter
                              n_top_hm, contr_2a, contr_1a),
               fontsize_col = 10, treeheight_row = 30, treeheight_col = 30)
      
      
      #### Barplot (A vs B):  --------------------------
      # número de DEGs por dirección
      # 0) Umbrales (usa los que ya tengas; si no existen, pongo defaults)
      if (!exists("alpha_padj")) alpha_padj <- 0.05
      if (!exists("lfc_thr"))    lfc_thr    <- 1
      
      # 1) Etiquetas dinámicas respecto a contr_2a (numerador del log2FC)
      lab_up_a   <- paste0("up in ",   contr_2a)
      lab_down_a <- paste0("down in ", contr_2a)
      
      # 2) Tabla de conteos SOLO con genes significativos
      tmp_a <- within(resdf_a, {
        sig <- !is.na(padj) & (padj < alpha_padj) & (abs(log2FoldChange) >= lfc_thr)
        direction <- ifelse(log2FoldChange > 0, lab_up_a,
                            ifelse(log2FoldChange < 0, lab_down_a, NA))
      })
      bar_df_a <- as.data.frame(table(tmp_a$direction[tmp_a$sig]))
      names(bar_df_a) <- c("dir", "n")
      bar_df_a$dir <- factor(bar_df_a$dir, levels = c(lab_up_a, lab_down_a))
      
      # 3) Colores (mismo esquema que MA/volcano: rojo = up, azul = down)
      cols_a <- setNames(c("#a5923d", "#aa4d6f"), c(lab_up_a, lab_down_a))
      
      # 4) Gráfico
      if (nrow(bar_df_a) > 0) {
        ggplot(bar_df_a, aes(x = dir, y = n, fill = dir)) +
          geom_col(width = 0.65) +
          geom_text(aes(label = n), vjust = -0.35, size = 4) +
          scale_fill_manual(values = cols_a, guide = "none") +
          labs(
            title    = paste0("DEGs por dirección (", contr_2a, " vs ", contr_1a, ")"),
            subtitle = paste0("padj < ", alpha_padj, "  &  (log2FC) ≥ ", lfc_thr),
            x        = paste0("Dirección (respecto a ", contr_2a, ")"),
            y        = "Número de genes"
          ) +
          ylim(0, max(bar_df_a$n) * 1.15) +
          theme_classic(base_size = 12) +
          theme(
            plot.title    = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30"),
            axis.text     = element_text(size = 10),
            axis.title    = element_text(size = 12)
          )
      } else {
        message("No hay DEGs con los umbrales (padj & log2FC) para el barplot.")
      }
      
      
      
      
      
      
      ## Contraste 2: A vs C  (log2FC = A/C) ---------------------------------------
      
      # Numerador / denominador para results()
      contr_1b <- CCC   # C = Control (denominador)
      contr_2b <- AAA   # A = numerador
      
      # Muestras del contraste: primero C (ref), luego A
      keep_samples_b <- c(samples_C, samples_A)
      stopifnot(all(keep_samples_b %in% colnames(counts_good)))
      
      counts_1b <- counts_good[, keep_samples_b, drop = FALSE]
      
      ## 2) colData (mismos niveles que el contraste)
      grp_vec_b <- setNames(rep(NA_character_, length(keep_samples_b)), keep_samples_b)
      grp_vec_b[samples_A] <- "A"
      grp_vec_b[samples_C] <- "C"
      
      coldata_b <- data.frame(Group = factor(grp_vec_b, levels = c(contr_1b, contr_2b)))
      rownames(coldata_b) <- names(grp_vec_b)
      
      ## 2.5) Chequeos + redondeo (recomendado por DESeq2)
      stopifnot(identical(colnames(counts_1b), rownames(coldata_b)))
      if (any(is.na(counts_1b))) stop("Hay NAs en counts_1b")
      
      bad_b <- which(abs(counts_1b - round(counts_1b)) > 1e-8, arr.ind = TRUE)
      cat("Celdas no enteras en counts_1b:", nrow(bad_b), "\n")
      counts_1b_int <- as.matrix(round(counts_1b))
      
      ## 3) DESeq2 (sin shrink)
      dds_b <- DESeqDataSetFromMatrix(countData = counts_1b_int, colData = coldata_b, design = ~ Group)
      dds_b <- dds_b[rowSums(counts(dds_b)) > 0, ]
      dds_b <- DESeq(dds_b)
      res_b <- results(dds_b, contrast = c("Group", contr_2b, contr_1b))  # A vs C
      resdf_b <- as.data.frame(res_b)
      resdf_b$gene <- rownames(resdf_b)
      
      # Umbrales si no existen
      if (!exists("alpha_padj")) alpha_padj <- 0.05
      if (!exists("lfc_thr"))    lfc_thr    <- 1
      
      # Etiquetas (dirección/significancia)
      resdf_b$direction <- ifelse(resdf_b$log2FoldChange > 0, paste0("up in ", contr_2b),
                                  ifelse(resdf_b$log2FoldChange < 0, paste0("down in ", contr_2b), "no change"))
      resdf_b$sig <- ifelse(!is.na(resdf_b$padj) & resdf_b$padj < alpha_padj, "significant", "ns")
      
      ## 4) VST para gráficos
      vst_b <- vst(dds_b, blind = TRUE)
      normcounts_b <- assay(vst_b)
      
      #### PCA (A vs C) -------------------------------------------
      # plotPCA usa por defecto ~500 genes más variables (no solo DEGs)
      pca_b <- DESeq2::plotPCA(vst_b, intgroup = "Group", returnData = TRUE)
      pca_b$Group <- factor(as.character(pca_b$Group), levels = c(contr_1b, contr_2b))
      
      # Colores: Control (C) = #3d6863, A = mismo tono que usaste para A
      col_ac <- setNames(c("#3d6863", "#a2593d"), c(contr_1b, contr_2b))
      
      ggplot(pca_b, aes(PC1, PC2, color = Group, label = name)) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(show.legend = FALSE, size = 3) +
        scale_color_manual(values = col_ac, breaks = c(contr_1b, contr_2b), name = "Group") +
        labs(
          title    = paste0("PCA (", contr_2b, " vs ", contr_1b, ")"),
          subtitle = "Transformación VST",
          x = paste0("PC1: ", round(100 * attr(pca_b, "percentVar")[1], 1), "%"),
          y = paste0("PC2: ", round(100 * attr(pca_b, "percentVar")[2], 1), "%")
        ) +
        theme_classic(base_size = 10) +
        theme(
          plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30"),
          axis.title    = element_text(size = 12),
          axis.text     = element_text(size = 12),
          legend.title  = element_text(size = 10),
          legend.text   = element_text(size = 10)
        )
      
      
      #### MA-plot (A vs C) --------------------------------------------------
      
      # 1) clasificar puntos (up/down/ns) y evitar Inf cuando padj = 0
      volcano_b <- transform(
        resdf_b,
        padj_plot = pmin(padj, 1),
        dir = ifelse(is.na(padj) | padj >= alpha_padj | abs(log2FoldChange) < lfc_thr,
                     "ns",
                     ifelse(log2FoldChange > 0,
                            paste0("up in ",   contr_2b),
                            paste0("down in ", contr_2b)))
      )
      
      # 2) paleta como en MA-plot
      cols_b <- setNames(  c("#aa4d6f", "grey75", "#a5923d"),
                           c(paste0("down in ", contr_2b), "ns", paste0("up in ", contr_2b))
      )
      
      
      # 3) gráfico
      library(ggplot2)
      ggplot(volcano_b, aes(x = log2FoldChange, y = -log10(padj_plot), color = dir)) +
        geom_point(alpha = 0.75, size = 1.2) +
        geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", colour = "grey40") +
        geom_hline(yintercept = -log10(alpha_padj),  linetype = "dashed", colour = "grey40") +
        scale_color_manual(values = cols_b, name = NULL) +
        labs(
          title    = "Volcano (A vs C)",
          subtitle = "DESeq2 • (padj < 0.05) & (log2FC ≥ 1)",
          x = paste0("log2FC (", contr_2b, " / ", contr_1b, ")"),
          y = "-log10(padj)"
        ) +
        theme_classic(base_size = 12) +
        theme(
          plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30"),
          axis.title    = element_text(size = 12),
          axis.text     = element_text(size = 12),
          legend.title  = element_text(size = 10),
          legend.text   = element_text(size = 10)
        )
      
      
      
      
      
      
      
      
      
      
      #### Volcano (A vs C) ------------------------------------------------
      ## MA-plot (A vs C) — idéntico al de A vs B
      ma_b <- transform(as.data.frame(res_b),
                        A = log10(baseMean + 1),
                        M = log2FoldChange)
      
      ma_b$flag <- with(ma_b,
                        ifelse(!is.na(padj) & padj < alpha_padj & abs(M) >= lfc_thr,
                               ifelse(M > 0, paste0("up in ",   contr_2b),
                                      paste0("down in ", contr_2b)),
                               "ns"))
      
      # MISMA paleta que usaste: ns gris, up verde, down morado
      cols_b <- setNames(c("#aa4d6f", "grey80", "#a5923d"),
                         c(paste0("down in ", contr_2b), "ns", paste0("up in ", contr_2b)))
      
      # MISMOS límites (ajústalos si quieres)
      x_lim_b <- c(0, 6)
      y_lim_b <- c(-12, 12)
      
      library(ggplot2)
      ggplot(ma_b, aes(A, M, color = flag)) +
        geom_point(size = 1.2, alpha = 0.7) +
        geom_hline(yintercept = 0,                  colour = "grey40") +
        geom_hline(yintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", colour = "grey40") +
        scale_color_manual(values = cols_b, name = NULL) +
        labs(title    = paste0("MA-plot (", contr_2b, " vs ", contr_1b, ")"),
             subtitle = sprintf("DESeq2 • (padj < %.2f) & |log2FC| ≥ %g", alpha_padj, lfc_thr),
             x = "log10(baseMean + 1)", y = "log2 Fold Change") +
        coord_cartesian(xlim = x_lim_b, ylim = y_lim_b) +
        theme_classic(base_size = 12) +
        theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30"),
              axis.title    = element_text(size = 12),
              axis.text     = element_text(size = 12),
              legend.title  = element_text(size = 10),
              legend.text   = element_text(size = 10))
      
      
      #### Heatmap (A vs C) ---------------------------------------------------
      # top genes por padj
      # Top N por padj (solo DEGs más significativos)
      n_top_hm <- 150
      ord_b    <- order(resdf_b$padj, na.last = NA)
      top_b    <- head(rownames(resdf_b)[ord_b], n_top_hm)
      
      mat_b0 <- normcounts_b[top_b, , drop = FALSE]
      mat_b  <- mat_b0 - rowMeans(mat_b0)  # centrado por gen
      
      ann_b <- data.frame(Group = coldata_b$Group, row.names = rownames(coldata_b))
      ann_b$Group <- factor(as.character(ann_b$Group), levels = c(contr_1b, contr_2b))
      
      # Colores de grupo: Control con #3d6863
      ann_colors_b <- list(Group = setNames(c("#3d6863", "#a2593d"), c(contr_1b, contr_2b)))
      
      # Gradiente simétrico alrededor de 0 (mismo esquema que A vs B)
      col_neg_b <- "#aa4d6f"; col_pos_b <- "#a5923d"
      nbreaks_b <- 256
      lim_b     <- max(abs(mat_b))
      breaks_b  <- seq(-lim_b, lim_b, length.out = nbreaks_b)
      hm_cols_b <- colorRampPalette(c(col_neg_b, "#f7efe6", col_pos_b))(nbreaks_b - 1)
      ticks_b   <- pretty(c(-lim_b, lim_b), n = 5)
      
      pheatmap(mat_b,
               color = hm_cols_b, breaks = breaks_b, scale = "none",
               border_color = NA, show_rownames = FALSE,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "correlation",
               clustering_method = "ward.D2",
               annotation_col    = ann_b,
               annotation_colors = ann_colors_b,
               legend_breaks = ticks_b,
               legend_labels = sprintf("%.1f", ticks_b),
               main = sprintf("Heatmap top %d en %s vs %s\n", n_top_hm, contr_2b, contr_1b),
               fontsize_col = 10, treeheight_row = 30, treeheight_col = 30)
      
      
      #### Barplot (A vs C) --------------------------------------------------
      lab_up_b   <- paste0("up in ",   contr_2b)
      lab_down_b <- paste0("down in ", contr_2b)
      
      tmp_b <- within(resdf_b, {
        sig <- !is.na(padj) & (padj < alpha_padj) & (abs(log2FoldChange) >= lfc_thr)
        direction <- ifelse(log2FoldChange > 0, lab_up_b,
                            ifelse(log2FoldChange < 0, lab_down_b, NA))
      })
      bar_df_b <- as.data.frame(table(tmp_b$direction[tmp_b$sig]))
      names(bar_df_b) <- c("dir", "n")
      bar_df_b$dir <- factor(bar_df_b$dir, levels = c(lab_up_b, lab_down_b))
      
      cols_bar_b <- setNames(c("#a5923d", "#aa4d6f"), c(lab_up_b, lab_down_b))
      
      if (nrow(bar_df_b) > 0) {
        ggplot(bar_df_b, aes(x = dir, y = n, fill = dir)) +
          geom_col(width = 0.65) +
          geom_text(aes(label = n), vjust = -0.35, size = 4) +
          scale_fill_manual(values = cols_bar_b, guide = "none") +
          labs(
            title    = paste0("DEGs por dirección (", contr_2b, " vs ", contr_1b, ")"),
            subtitle = paste0("(padj < ", alpha_padj, ")", "  &  (log2FC ≥ ", lfc_thr, ")"),
            x        = paste0("Dirección (respecto a ", contr_2b, ")"),
            y        = "Número de genes"
          ) +
          ylim(0, max(bar_df_b$n) * 1.15) +
          theme_classic(base_size = 12) +
          theme(
            plot.title    = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30")
          )
      } else {
        message("No hay DEGs con los umbrales (padj & log2FC) para el barplot (A vs C).")
      }
      
      
      
      
      
      ## Contraste 3: B vs C  (log2FC = B/C) ---------------------------------------
      # Qué se ve en cada gráfico (resumen):
      # - PCA: TODAS las muestras (vst) de los grupos del contraste.
      # - MA-plot: TODOS los genes analizados; se resaltan los DEGs según padj y |log2FC|.
      # - Volcano: TODOS los genes; se resaltan los DEGs.
      # - Heatmap: SOLO los top N DEGs (por padj).
      # - Barplot: SOLO cuenta DEGs (por dirección up/down).
      
      ## 1) Parámetros y subconjunto (B vs C)
      contr_1c <- CCC   # referencia (denominador, C)
      contr_2c <- BBB   # numerador  (B)
      keep_samples_c <- c(samples_C, samples_B)
      stopifnot(all(keep_samples_c %in% colnames(counts_good)))
      counts_1c <- counts_good[, keep_samples_c, drop = FALSE]
      
      ## 2) colData con niveles en orden (C, B)
      grp_vec_c <- setNames(rep(NA_character_, length(keep_samples_c)), keep_samples_c)
      grp_vec_c[samples_B] <- "B"
      grp_vec_c[samples_C] <- "C"
      coldata_c <- data.frame(Group = factor(grp_vec_c, levels = c(contr_1c, contr_2c)))
      rownames(coldata_c) <- names(grp_vec_c)
      
      ## 2.5) Chequeos y redondeo
      stopifnot(identical(colnames(counts_1c), rownames(coldata_c)))
      if (any(is.na(counts_1c))) stop("Hay NAs en counts_1c")
      bad_c <- which(abs(counts_1c - round(counts_1c)) > 1e-8, arr.ind = TRUE)
      cat("Celdas no enteras en counts_1c:", nrow(bad_c), "\n")
      counts_1c_int <- as.matrix(round(counts_1c))
      
      ## 3) DESeq2 (sin shrink)
      dds_c  <- DESeqDataSetFromMatrix(countData = counts_1c_int, colData = coldata_c, design = ~ Group)
      dds_c  <- dds_c[rowSums(counts(dds_c)) > 0, ]
      dds_c  <- DESeq(dds_c)
      res_c  <- results(dds_c, contrast = c("Group", contr_2c, contr_1c))  # B vs C
      resdf_c <- as.data.frame(res_c)
      resdf_c$gene <- rownames(resdf_c)
      
      ## Etiquetas up/down y significancia (para barplot/volcano)
      if (!exists("alpha_padj")) alpha_padj <- 0.05
      if (!exists("lfc_thr"))    lfc_thr    <- 1
      resdf_c$direction <- ifelse(resdf_c$log2FoldChange > 0, paste0("up in ", contr_2c),
                                  ifelse(resdf_c$log2FoldChange < 0, paste0("down in ", contr_2c), "no change"))
      resdf_c$sig <- ifelse(!is.na(resdf_c$padj) & resdf_c$padj < alpha_padj, "significant", "ns")
      
      ## 4) Normalización para gráficos (VST)
      vst_c        <- vst(dds_c, blind = TRUE)
      normcounts_c <- assay(vst_c)
      
      #### PCA (B vs C) --------------------------------------------------------------
      # Muestras (TODAS las del contraste) en espacio PCA; color por Group.
      pca_c <- DESeq2::plotPCA(vst_c, intgroup = "Group", returnData = TRUE)
      pca_c$Group <- factor(as.character(pca_c$Group), levels = c(contr_1c, contr_2c))
      # Paleta fija: Control = #3d6863, B = #fdac4f
      col_bc <- setNames(c("#3d6863", "#fdac4f"), c(contr_1c, contr_2c))
      
      ggplot(pca_c, aes(PC1, PC2, color = Group, label = name)) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(show.legend = FALSE, size = 3) +
        scale_color_manual(values = col_bc, breaks = c(contr_1c, contr_2c), name = "Group") +
        labs(
          title    = paste0("PCA (", contr_2c, " vs ", contr_1c, ")"),
          subtitle = "Transformación VST",
          x = paste0("PC1: ", round(100 * attr(pca_c, "percentVar")[1], 1), "%"),
          y = paste0("PC2: ", round(100 * attr(pca_c, "percentVar")[2], 1), "%")
        ) +
        theme_classic(base_size = 12) +
        theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30"),
              axis.title    = element_text(size = 12),
              axis.text     = element_text(size = 12),
              legend.title  = element_text(size = 10),
              legend.text   = element_text(size = 10))
      
      #### MA-plot (B vs C) ----------------------------------------------------------
      # TODOS los genes; resaltamos DEGs (padj & |log2FC|).
      ma_c <- transform(as.data.frame(res_c),
                        A = log10(baseMean + 1),
                        M = log2FoldChange)
      
      ma_c$flag <- with(ma_c,
                        ifelse(!is.na(padj) & padj < alpha_padj & abs(M) >= lfc_thr,
                               ifelse(M > 0, paste0("up in ",   contr_2c),
                                      paste0("down in ", contr_2c)),
                               "ns"))
      
      cols_c <- setNames(c("#aa4d6f", "grey80", "#a5923d"),
                         c(paste0("down in ", contr_2c), "ns", paste0("up in ", contr_2c)))
      x_lim_c <- c(0, 6)
      y_lim_c <- c(-12, 12)
      
      ggplot(ma_c, aes(A, M, color = flag)) +
        geom_point(size = 1.2, alpha = 0.7) +
        geom_hline(yintercept = 0, colour = "grey40") +
        geom_hline(yintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", colour = "grey40") +
        scale_color_manual(values = cols_c, name = NULL) +
        labs(title    = paste0("MA-plot (", contr_2c, " vs ", contr_1c, ")"),
             subtitle = sprintf("DESeq2 • (padj < %.2f) & (log2FC ≥ %g)", alpha_padj, lfc_thr),
             x = "log10(baseMean + 1)", y = "log2 Fold Change") +
        coord_cartesian(xlim = x_lim_c, ylim = y_lim_c) +
        theme_classic(base_size = 12) +
        theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30"),
              axis.title    = element_text(size = 12),
              axis.text     = element_text(size = 12),
              legend.title  = element_text(size = 10),
              legend.text   = element_text(size = 10))
      
      #### Volcano (B vs C) ----------------------------------------------------------
      # TODOS los genes; resaltamos DEGs.
      volcano_c <- transform(
        resdf_c,
        padj_plot = pmin(padj, 1),
        dir = ifelse(is.na(padj) | padj >= alpha_padj | abs(log2FoldChange) < lfc_thr,
                     "ns",
                     ifelse(log2FoldChange > 0,
                            paste0("up in ",   contr_2c),
                            paste0("down in ", contr_2c)))
      )
      cols_c <- setNames(c("#aa4d6f", "grey75", "#a5923d"),
                         c(paste0("down in ", contr_2c), "ns", paste0("up in ", contr_2c)))
      
      ggplot(volcano_c, aes(x = log2FoldChange, y = -log10(padj_plot), color = dir)) +
        geom_point(alpha = 0.75, size = 1.2) +
        geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", colour = "grey40") +
        geom_hline(yintercept = -log10(alpha_padj),  linetype = "dashed", colour = "grey40") +
        scale_color_manual(values = cols_c, name = NULL) +
        labs(title    = paste0("Volcano (", contr_2c, " vs ", contr_1c, ")"),
             subtitle = sprintf("DESeq2 • (padj < %.2f) & (log2FC ≥ %g)", alpha_padj, lfc_thr),
             x = paste0("log2FC (", contr_2c, " / ", contr_1c, ")"),
             y = "-log10(padj)") +
        theme_classic(base_size = 12) +
        theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30"),
              axis.title    = element_text(size = 12),
              axis.text     = element_text(size = 12),
              legend.title  = element_text(size = 10),
              legend.text   = element_text(size = 10))
      
      #### Heatmap (B vs C) -----------------------------------------------
      # top genes
      # SOLO los top N DEGs (por padj); z-score por gen (centrado a media 0).
      n_top_hm <- if (exists("n_top_hm")) n_top_hm else 150
      ord_c    <- order(resdf_c$padj, na.last = NA)
      top_c    <- head(rownames(resdf_c)[ord_c], n_top_hm)
      
      mat_c0 <- normcounts_c[top_c, , drop = FALSE]
      mat_c  <- mat_c0 - rowMeans(mat_c0)
      
      ann_c <- data.frame(Group = coldata_c$Group, row.names = rownames(coldata_c))
      ann_c$Group <- factor(as.character(ann_c$Group), levels = c(contr_1c, contr_2c))
      
      # Anotación de columnas (Group): Control (#3d6863), B (#fdac4f)
      ann_colors_c <- list(Group = setNames(c("#3d6863", "#fdac4f"), c(contr_1c, contr_2c)))
      
      # Gradiente simétrico (mismos colores que usaste antes)
      nbreaks_c <- 256
      lim_c     <- max(abs(mat_c))
      breaks_c  <- seq(-lim_c, lim_c, length.out = nbreaks_c)
      hm_cols_c <- colorRampPalette(c("#aa4d6f", "#f7efe6", "#a5923d"))(nbreaks_c - 1)
      ticks_c   <- pretty(c(-lim_c, lim_c), n = 5)
      
      pheatmap(mat_c,
               color = hm_cols_c, breaks = breaks_c,
               scale = "none", border_color = NA,
               show_rownames = FALSE,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "correlation",
               clustering_method = "ward.D2",
               annotation_col    = ann_c,
               annotation_colors = ann_colors_c,
               legend_breaks = ticks_c,
               legend_labels = sprintf("%.1f", ticks_c),
               main = sprintf("Heatmap top %d en %s vs %s\n", n_top_hm, contr_2c, contr_1c),
               fontsize_col = 10, treeheight_row = 30, treeheight_col = 30)
      
      #### Barplot (B vs C):  -----------------------------
      # número de DEGs por dirección
      lab_up_c   <- paste0("up in ",   contr_2c)
      lab_down_c <- paste0("down in ", contr_2c)
      
      tmp_c <- within(resdf_c, {
        sig <- !is.na(padj) & (padj < alpha_padj) & (abs(log2FoldChange) >= lfc_thr)
        direction <- ifelse(log2FoldChange > 0, lab_up_c,
                            ifelse(log2FoldChange < 0, lab_down_c, NA))
      })
      bar_df_c <- as.data.frame(table(tmp_c$direction[tmp_c$sig]))
      names(bar_df_c) <- c("dir", "n")
      bar_df_c$dir <- factor(bar_df_c$dir, levels = c(lab_up_c, lab_down_c))
      
      cols_bar_c <- setNames(c("#a5923d", "#aa4d6f"), c(lab_up_c, lab_down_c))
      
      if (nrow(bar_df_c) > 0) {
        ggplot(bar_df_c, aes(x = dir, y = n, fill = dir)) +
          geom_col(width = 0.65) +
          geom_text(aes(label = n), vjust = -0.35, size = 4) +
          scale_fill_manual(values = cols_bar_c, guide = "none") +
          labs(title    = paste0("DEGs por dirección (", contr_2c, " vs ", contr_1c, ")"),
               subtitle = paste0("padj < ", alpha_padj, "  &  (log2FC ≥ ", lfc_thr, ")"),
               x = paste0("Dirección (respecto a ", contr_2c, ")"),
               y = "Número de genes") +
          ylim(0, max(bar_df_c$n) * 1.15) +
          theme_classic(base_size = 12) +
          theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
                plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30"),
                axis.text     = element_text(size = 10),
                axis.title    = element_text(size = 12))
      } else {
        message("No hay DEGs con los umbrales (padj & log2FC) para el barplot (B vs C).")
      }
      
      
      
      
      
      # NEURONS B/A-------------
      ##  Etiquetas del contraste (B/A) 
      # bloque de “ENRIQ… B/A” ya definimos:
      # contr_1a <- AAA  # A (denominador)
      # contr_2a <- BBB  # B (numerador)
      tag_a <- paste0(contr_2a, "_vs_", contr_1a)   # "B_vs_A"
      
      ##  0) Máscara de DEGs por FDR 
      alpha_a   <- 0.05
      padj_a    <- res_a$padj
      padj_a[is.na(padj_a)] <- 1
      deg_mask_a <- padj_a < alpha_a
      
      cat("N NAs en res_a$padj:", sum(is.na(res_a$padj)), "\n")
      cat(tag_a, "DEGs at FDR <", alpha_a, ": ", sum(deg_mask_a), "\n")
      
      if (!any(deg_mask_a)) stop("No DEGs at FDR <", alpha_a, " en ", tag_a)
      
      ## 1) FPM normalizados de ESTE dds (a) 
      fpm_norm_a <- fpm(dds_a)                                   # genes x muestras
      DEG_cpm_a  <- fpm_norm_a[deg_mask_a, , drop = FALSE]
      cat("[", tag_a, "] DEG_cpm_a:", nrow(DEG_cpm_a), "genes x", ncol(DEG_cpm_a), "muestras\n")
      
      ## 2) log2(+1) 
      log2_deg_a <- log2(DEG_cpm_a + 1)
      
      ## 3) z-score por gen 
      cs.log2_a <- t(scale(t(log2_deg_a)))
      
      ## 4) checks 
      stopifnot(all(is.finite(cs.log2_a)))
      stopifnot(identical(colnames(cs.log2_a), rownames(coldata_a)))
      cat("z-score range (", tag_a, "): [",
          round(min(cs.log2_a),2), ", ", round(max(cs.log2_a),2), "]\n", sep="")
      
      ## (opcional) guardado
      # write.csv(as.data.frame(DEG_cpm_a), file=paste0("DEGs_FPMnorm_", tag_a, ".csv"))
      # write.csv(as.data.frame(cs.log2_a), file=paste0("DEGs_log2p1_zscore_", tag_a, ".csv"))
      
      
      
      ### SOMs B/A ----------
      #Self-Organizing Map 
      
      
      ## 1) Sanity check
      stopifnot(exists("cs.log2_a"), nrow(cs.log2_a) > 0)
      head(cs.log2_a)
      
      ## 2) Fijar aleatoriedad
      set.seed(42)
      
      ## 3) Definición de la grilla
      som_grid_axes_a <- 3
      som_grid_a <- kohonen::somgrid(xdim = som_grid_axes_a,
                                     ydim = som_grid_axes_a,
                                     topo = "hexagonal",
                                     toroidal = TRUE)
      
      ## 4) Entrenamiento del SOM
      som_fit_a <- kohonen::som(
        as.matrix(cs.log2_a),
        grid      = som_grid_a,
        rlen      = 200,
        keep.data = TRUE
      )
      
      ## 5) Clasificación de genes y codebooks
      classif_a <- som_fit_a$unit.classif
      names(classif_a) <- rownames(cs.log2_a)
      codes_a <- som_fit_a$codes[[1]]   # neuronas × muestras
      
      ## 6) 
      ### Plots rápidos del SOM
      plot(som_fit_a,                          main = paste0("SOM – Codebook vectors at ",   tag_a))
      plot(som_fit_a, type = "mapping",        main = paste0("SOM – Mapping at ",            tag_a), pch = ".")
      plot(som_fit_a, type = "changes",        main = paste0("SOM – Training at ",           tag_a), pch = ".")
      plot(som_fit_a, type = "counts",         main = paste0("SOM – Genes per neuron at ",   tag_a), pch = ".")
      plot(som_fit_a, type = "dist.neighbours",main = paste0("SOM – U-matrix at ",           tag_a), pch = ".")
      plot(som_fit_a, type = "quality",        main = paste0("SOM – Quantization error at ", tag_a), pch = ".")
      
      
      
      ### Heatmaps SOMs & Cluster codebooks B/A ------
      
      
      ## 0) Alias de metadatos (para no pisar objetos)
      coldata_a2 <- coldata_a
      coldata_a2$Group <- droplevels(coldata_a2$Group)
      
      if (!exists("pal_grupo")) pal_grupo <- c(A="#fdac4f", B="#a2593d", C="#3d6863")
      pal_a <- pal_grupo[levels(coldata_a2$Group)]
      
      ## 1) Alinear columnas de codes con coldata_a2
      stopifnot(all(colnames(codes_a) %in% rownames(coldata_a2)))
      codes_a <- codes_a[, rownames(coldata_a2), drop = FALSE]
      stopifnot(identical(colnames(codes_a), rownames(coldata_a2)))
      
      ## 2) Clustering de muestras y neuronas
      clust.sample_a <- hclust(dist(t(codes_a), method = "euclidean"),  method = "ward.D2")
      clust.neuron_a <- hclust(dist(codes_a,   method = "manhattan"),   method = "ward.D2")
      
      ## 3) Anotación de filas (muestras) + colores
      ann_row_a <- data.frame(Group = coldata_a2$Group, row.names = rownames(coldata_a2))
      ann_colors_row_a <- list(Group = pal_a)
      
      
      # 1) Colores extremos (ajústalos si quieres)
      col_neg_a <- "#aa4d6f"  # morado (menor que la media de su gen)
      col_pos_a <- "#a5923d"  # verde (mayor que la media de su gen)
      
      # 2) Matriz VST normalizada + centrado por gen (z-score por fila sin escalar sd)
      mat_a0 <- normcounts_a[top_a, , drop = FALSE]
      mat_a  <- mat_a0 - rowMeans(mat_a0)
      
      # 2) Degradado (blanco en 0 para que el centro sea neutro)
      nbreaks_a <- 256
      lim_a     <- max(abs(mat_a))  # simétrico alrededor de 0
      breaks_a  <- seq(-lim_a, lim_a, length.out = nbreaks_a)
      hm_cols_a <- colorRampPalette(c(col_neg_a, "#f7efe6", col_pos_a))(nbreaks_a - 1)
      
      # 3) Ticks de la barra de color (leyenda de la escala)
      ticks_a <- pretty(c(-lim_a, lim_a), n = 5)
      
      
      #### Heatmap 1 (samples × neurons) -------
      pheatmap(t(codes_a),
               color = hm_cols_a, breaks = breaks_a,
               border_color  = "grey60",
               scale         = "column",
               show_rownames = TRUE,
               show_colnames = FALSE,
               cluster_rows  = clust.sample_a,  cutree_rows = 2,
               annotation_row = ann_row_a,
               cluster_cols  = clust.neuron_a,  cutree_cols = 3,
               annotation_colors = ann_colors_row_a,
               legend_breaks = ticks_a,
               legend_labels = sprintf("%.1f", ticks_a),
               main = paste0("SOM codebooks — samples × neurons (", tag_a, ")"))
      
      ## 5) Clusters de neuronas
      k_neuron_a  <- 3
      clust.pat_a <- cutree(clust.neuron_a, k = k_neuron_a)
      clust.aux_a <- paste0("C", clust.pat_a); names(clust.aux_a) <- names(clust.pat_a)
      
      ## 6) Anotación de columnas (neurons)
      ann_col_a <- data.frame(Pattern = clust.aux_a, row.names = names(clust.aux_a))
      
      ## 7) Paleta para los clusters de neuronas
      pal_nc_a <- colorRampPalette(brewer.pal(8, "Accent"))(k_neuron_a)
      names(pal_nc_a) <- paste0("C", seq_len(k_neuron_a))
      
      ## 8) Colores combinados
      ann_colors_a <- list(Group = pal_a, Pattern = pal_nc_a)
      
      ## 9) Forzar niveles
      ann_col_a$Pattern <- factor(ann_col_a$Pattern, levels = names(ann_colors_a$Pattern))
      
      #### Heatmap 2 clusters de neuronas ----
      pheatmap(t(codes_a),
               color = hm_cols_a, breaks = breaks_a,
               border_color  = "grey60",
               scale         = "column",
               show_rownames = TRUE,
               show_colnames = FALSE,
               cluster_rows  = clust.sample_a,  cutree_rows = 2,
               annotation_row = ann_row_a,
               cluster_cols  = clust.neuron_a,  cutree_cols = k_neuron_a,
               annotation_col = ann_col_a,
               annotation_colors = ann_colors_a,
               legend_breaks = ticks_a,
               legend_labels = sprintf("%.1f", ticks_a),
               main = paste0("SOM codebooks — with neuron clusters (", tag_a, ")"))
      
      ## 11) Gene → cluster de neurona
      type.pattern_a <- clust.pat_a[classif_a]
      names(type.pattern_a) <- names(classif_a)
      cat("Genes per neuron cluster (", tag_a, "):\n", sep = ""); print(table(type.pattern_a))
      
      ## (opcional) export
      # write.csv(codes_a, file = paste0("SOM_codebooks_neuron_by_sample_", tag_a, ".csv"))
      # write.csv(data.frame(gene = names(type.pattern_a),
      #                      neuron_cluster = paste0("C", type.pattern_a)),
      #           file = paste0("genes_to_SOM_clusters_", tag_a, ".csv"), row.names = FALSE)
      
      
      
      
      
      
      # NEURONS A/C-------------
      ##  Etiquetas del contraste (A/C)
      # Deben existir de tu bloque “ENRIQ… A/C”:
      # contr_1b <- CCC  # C (denominador / Control)
      # contr_2b <- AAA  # A (numerador)
      tag_b <- paste0(contr_2b, "_vs_", contr_1b)   # "A_vs_C"
      
      ##  0) Máscara de DEGs por FDR
      alpha_b   <- 0.05
      padj_b    <- res_b$padj
      padj_b[is.na(padj_b)] <- 1
      deg_mask_b <- padj_b < alpha_b
      
      cat("N NAs en res_b$padj:", sum(is.na(res_b$padj)), "\n")
      cat(tag_b, "DEGs at FDR <", alpha_b, ": ", sum(deg_mask_b), "\n")
      if (!any(deg_mask_b)) stop("No DEGs at FDR <", alpha_b, " en ", tag_b)
      
      ## 1) FPM normalizados de ESTE dds (b)
      fpm_norm_b <- fpm(dds_b)                                   # genes x muestras
      DEG_cpm_b  <- fpm_norm_b[deg_mask_b, , drop = FALSE]
      cat("[", tag_b, "] DEG_cpm_b:", nrow(DEG_cpm_b), "genes x", ncol(DEG_cpm_b), "muestras\n")
      
      ## 2) log2(+1)
      log2_deg_b <- log2(DEG_cpm_b + 1)
      
      ## 3) z-score por gen
      cs.log2_b <- t(scale(t(log2_deg_b)))
      
      ## 4) checks
      stopifnot(all(is.finite(cs.log2_b)))
      stopifnot(identical(colnames(cs.log2_b), rownames(coldata_b)))
      cat("z-score range (", tag_b, "): [",
          round(min(cs.log2_b),2), ", ", round(max(cs.log2_b),2), "]\n", sep="")
      
      ## (opcional) guardado
      # write.csv(as.data.frame(DEG_cpm_b), file=paste0("DEGs_FPMnorm_", tag_b, ".csv"))
      # write.csv(as.data.frame(cs.log2_b), file=paste0("DEGs_log2p1_zscore_", tag_b, ".csv"))
      
      
      ### SOMs A/C ----------
      ## 1) Sanity check
      stopifnot(exists("cs.log2_b"), nrow(cs.log2_b) > 0)
      
      ## 2) Fijar aleatoriedad
      set.seed(42)
      
      ## 3) Definición de la grilla
      som_grid_axes_b <- 3
      som_grid_b <- kohonen::somgrid(xdim = som_grid_axes_b,
                                     ydim = som_grid_axes_b,
                                     topo = "hexagonal",
                                     toroidal = TRUE)
      
      ## 4) Entrenamiento del SOM
      som_fit_b <- kohonen::som(
        as.matrix(cs.log2_b),
        grid      = som_grid_b,
        rlen      = 200,
        keep.data = TRUE
      )
      
      ## 5) Clasificación de genes y codebooks
      classif_b <- som_fit_b$unit.classif
      names(classif_b) <- rownames(cs.log2_b)
      codes_b <- som_fit_b$codes[[1]]   # neuronas × muestras
      
      ## 6) Plots rápidos del SOM
      plot(som_fit_b,                           main = paste0("SOM – Codebook vectors at ",   tag_b))
      plot(som_fit_b, type = "mapping",         main = paste0("SOM – Mapping at ",            tag_b), pch = ".")
      plot(som_fit_b, type = "changes",         main = paste0("SOM – Training at ",           tag_b), pch = ".")
      plot(som_fit_b, type = "counts",          main = paste0("SOM – Genes per neuron at ",   tag_b), pch = ".")
      plot(som_fit_b, type = "dist.neighbours", main = paste0("SOM – U-matrix at ",           tag_b), pch = ".")
      plot(som_fit_b, type = "quality",         main = paste0("SOM – Quantization error at ", tag_b), pch = ".")
      
      
      ### Heatmaps SOMs & Cluster codebooks A/C ------
      ## 0) Alias de metadatos (para no pisar objetos)
      coldata_b2 <- coldata_b
      coldata_b2$Group <- droplevels(coldata_b2$Group)
      
      # Paleta global (incluye Control)
      if (!exists("pal_grupo")) pal_grupo <- c(A="#fdac4f", B="#a2593d", C="#3d6863")
      pal_b <- pal_grupo[levels(coldata_b2$Group)]
      
      ## 1) Alinear columnas de codes con coldata_b2
      stopifnot(all(colnames(codes_b) %in% rownames(coldata_b2)))
      codes_b <- codes_b[, rownames(coldata_b2), drop = FALSE]
      stopifnot(identical(colnames(codes_b), rownames(coldata_b2)))
      
      ## 2) Clustering de muestras y neuronas
      clust.sample_b <- hclust(dist(t(codes_b), method = "euclidean"),  method = "ward.D2")
      clust.neuron_b <- hclust(dist(codes_b,   method = "manhattan"),   method = "ward.D2")
      
      ## 3) Anotación de filas (muestras) + colores
      ann_row_b <- data.frame(Group = coldata_b2$Group, row.names = rownames(coldata_b2))
      ann_colors_row_b <- list(Group = pal_b)
      
      # Colores para gradiente (mismo estilo que ‘a’)
      col_neg_b <- "#aa4d6f"
      col_pos_b <- "#a5923d"
      nbreaks_b <- 256
      lim_b     <- 3                   # z-score ~ [-3, 3] tras 'scale="column"'
      breaks_b  <- seq(-lim_b, lim_b, length.out = nbreaks_b)
      hm_cols_b <- colorRampPalette(c(col_neg_b, "#f7efe6", col_pos_b))(nbreaks_b - 1)
      ticks_b   <- pretty(c(-lim_b, lim_b), n = 5)
      
      #### Heatmap 1 (samples × neurons) -------
      pheatmap(t(codes_b),
               color = hm_cols_b, breaks = breaks_b,
               border_color  = "grey60",
               scale         = "column",
               show_rownames = TRUE,
               show_colnames = FALSE,
               cluster_rows  = clust.sample_b,  cutree_rows = 2,
               annotation_row = ann_row_b,
               cluster_cols  = clust.neuron_b,  cutree_cols = 3,
               annotation_colors = ann_colors_row_b,
               legend_breaks = ticks_b,
               legend_labels = sprintf("%.1f", ticks_b),
               main = paste0("SOM codebooks — samples × neurons (", tag_b, ")"))
      
      ## 5) Clusters de neuronas
      k_neuron_b  <- 3
      clust.pat_b <- cutree(clust.neuron_b, k = k_neuron_b)
      clust.aux_b <- paste0("C", clust.pat_b); names(clust.aux_b) <- names(clust.pat_b)
      
      ## 6) Anotación de columnas (neurons)
      ann_col_b <- data.frame(Pattern = clust.aux_b, row.names = names(clust.aux_b))
      
      ## 7) Paleta para los clusters de neuronas
      pal_nc_b <- colorRampPalette(brewer.pal(8, "Accent"))(k_neuron_b)
      names(pal_nc_b) <- paste0("C", seq_len(k_neuron_b))
      
      ## 8) Colores combinados
      ann_colors_b <- list(Group = pal_b, Pattern = pal_nc_b)
      
      ## 9) Forzar niveles
      ann_col_b$Pattern <- factor(ann_col_b$Pattern, levels = names(ann_colors_b$Pattern))
      
      #### Heatmap 2 clusters de neuronas ----
      pheatmap(t(codes_b),
               color = hm_cols_b, breaks = breaks_b,
               border_color  = "grey60",
               scale         = "column",
               show_rownames = TRUE,
               show_colnames = FALSE,
               cluster_rows  = clust.sample_b,  cutree_rows = 2,
               annotation_row = ann_row_b,
               cluster_cols  = clust.neuron_b,  cutree_cols = k_neuron_b,
               annotation_col = ann_col_b,
               annotation_colors = ann_colors_b,
               legend_breaks = ticks_b,
               legend_labels = sprintf("%.1f", ticks_b),
               main = paste0("SOM codebooks — with neuron clusters (", tag_b, ")"))
      
      ## 11) Gene → cluster de neurona
      type.pattern_b <- clust.pat_b[classif_b]
      names(type.pattern_b) <- names(classif_b)
      cat("Genes per neuron cluster (", tag_b, "):\n", sep = ""); print(table(type.pattern_b))
      
      
      
      
      
      
      
      # NEURONS B/C-------------
      ##  Etiquetas del contraste (B/C)
      # Deben existir de tu bloque “ENRIQ… B/C”:
      # contr_1c <- CCC  # C (denominador / Control)
      # contr_2c <- BBB  # B (numerador)
      tag_c <- paste0(contr_2c, "_vs_", contr_1c)   # "B_vs_C"
      
      ##  0) Máscara de DEGs por FDR
      alpha_c   <- 0.05
      padj_c    <- res_c$padj
      padj_c[is.na(padj_c)] <- 1
      deg_mask_c <- padj_c < alpha_c
      
      cat("N NAs en res_c$padj:", sum(is.na(res_c$padj)), "\n")
      cat(tag_c, "DEGs at FDR <", alpha_c, ": ", sum(deg_mask_c), "\n")
      if (!any(deg_mask_c)) stop("No DEGs at FDR <", alpha_c, " en ", tag_c)
      
      ## 1) FPM normalizados de ESTE dds (c)
      fpm_norm_c <- fpm(dds_c)                                   # genes x muestras
      DEG_cpm_c  <- fpm_norm_c[deg_mask_c, , drop = FALSE]
      cat("[", tag_c, "] DEG_cpm_c:", nrow(DEG_cpm_c), "genes x", ncol(DEG_cpm_c), "muestras\n")
      
      ## 2) log2(+1)
      log2_deg_c <- log2(DEG_cpm_c + 1)
      
      ## 3) z-score por gen
      cs.log2_c <- t(scale(t(log2_deg_c)))
      
      ## 4) checks
      stopifnot(all(is.finite(cs.log2_c)))
      stopifnot(identical(colnames(cs.log2_c), rownames(coldata_c)))
      cat("z-score range (", tag_c, "): [",
          round(min(cs.log2_c),2), ", ", round(max(cs.log2_c),2), "]\n", sep="")
      
      ## (opcional) guardado
      # write.csv(as.data.frame(DEG_cpm_c), file=paste0("DEGs_FPMnorm_", tag_c, ".csv"))
      # write.csv(as.data.frame(cs.log2_c), file=paste0("DEGs_log2p1_zscore_", tag_c, ".csv"))
      
      
      ### SOMs B/C ----------
      ## 1) Sanity check
      stopifnot(exists("cs.log2_c"), nrow(cs.log2_c) > 0)
      
      ## 2) Fijar aleatoriedad
      set.seed(42)
      
      ## 3) Definición de la grilla
      som_grid_axes_c <- 3
      som_grid_c <- kohonen::somgrid(xdim = som_grid_axes_c,
                                     ydim = som_grid_axes_c,
                                     topo = "hexagonal",
                                     toroidal = TRUE)
      
      ## 4) Entrenamiento del SOM
      som_fit_c <- kohonen::som(
        as.matrix(cs.log2_c),
        grid      = som_grid_c,
        rlen      = 200,
        keep.data = TRUE
      )
      
      ## 5) Clasificación de genes y codebooks
      classif_c <- som_fit_c$unit.classif
      names(classif_c) <- rownames(cs.log2_c)
      codes_c <- som_fit_c$codes[[1]]   # neuronas × muestras
      
      ## 6) Plots rápidos del SOM
      plot(som_fit_c,                           main = paste0("SOM – Codebook vectors at ",   tag_c))
      plot(som_fit_c, type = "mapping",         main = paste0("SOM – Mapping at ",            tag_c), pch = ".")
      plot(som_fit_c, type = "changes",         main = paste0("SOM – Training at ",           tag_c), pch = ".")
      plot(som_fit_c, type = "counts",          main = paste0("SOM – Genes per neuron at ",   tag_c), pch = ".")
      plot(som_fit_c, type = "dist.neighbours", main = paste0("SOM – U-matrix at ",           tag_c), pch = ".")
      plot(som_fit_c, type = "quality",         main = paste0("SOM – Quantization error at ", tag_c), pch = ".")
      
      
      ### Heatmaps SOMs & Cluster codebooks B/C ------
      ## 0) Alias de metadatos (para no pisar objetos)
      coldata_c2 <- coldata_c
      coldata_c2$Group <- droplevels(coldata_c2$Group)
      
      # Paleta global (incluye Control)
      if (!exists("pal_grupo")) pal_grupo <- c(A="#fdac4f", B="#a2593d", C="#3d6863")
      pal_c <- pal_grupo[levels(coldata_c2$Group)]
      
      ## 1) Alinear columnas de codes con coldata_c2
      stopifnot(all(colnames(codes_c) %in% rownames(coldata_c2)))
      codes_c <- codes_c[, rownames(coldata_c2), drop = FALSE]
      stopifnot(identical(colnames(codes_c), rownames(coldata_c2)))
      
      ## 2) Clustering de muestras y neuronas
      clust.sample_c <- hclust(dist(t(codes_c), method = "euclidean"),  method = "ward.D2")
      clust.neuron_c <- hclust(dist(codes_c,   method = "manhattan"),   method = "ward.D2")
      
      ## 3) Anotación de filas (muestras) + colores
      ann_row_c <- data.frame(Group = coldata_c2$Group, row.names = rownames(coldata_c2))
      ann_colors_row_c <- list(Group = pal_c)
      
      # Colores para gradiente
      col_neg_c <- "#aa4d6f"
      col_pos_c <- "#a5923d"
      nbreaks_c <- 256
      lim_c     <- 3
      breaks_c  <- seq(-lim_c, lim_c, length.out = nbreaks_c)
      hm_cols_c <- colorRampPalette(c(col_neg_c, "#f7efe6", col_pos_c))(nbreaks_c - 1)
      ticks_c   <- pretty(c(-lim_c, lim_c), n = 5)
      
      #### Heatmap 1 (samples × neurons) -------
      pheatmap(t(codes_c),
               color = hm_cols_c, breaks = breaks_c,
               border_color  = "grey60",
               scale         = "column",
               show_rownames = TRUE,
               show_colnames = FALSE,
               cluster_rows  = clust.sample_c,  cutree_rows = 2,
               annotation_row = ann_row_c,
               cluster_cols  = clust.neuron_c,  cutree_cols = 3,
               annotation_colors = ann_colors_row_c,
               legend_breaks = ticks_c,
               legend_labels = sprintf("%.1f", ticks_c),
               main = paste0("SOM codebooks — samples × neurons (", tag_c, ")"))
      
      ## 5) Clusters de neuronas
      k_neuron_c  <- 3
      clust.pat_c <- cutree(clust.neuron_c, k = k_neuron_c)
      clust.aux_c <- paste0("C", clust.pat_c); names(clust.aux_c) <- names(clust.pat_c)
      
      ## 6) Anotación de columnas (neurons)
      ann_col_c <- data.frame(Pattern = clust.aux_c, row.names = names(clust.aux_c))
      
      ## 7) Paleta para los clusters de neuronas
      pal_nc_c <- colorRampPalette(brewer.pal(8, "Accent"))(k_neuron_c)
      names(pal_nc_c) <- paste0("C", seq_len(k_neuron_c))
      
      ## 8) Colores combinados
      ann_colors_c <- list(Group = pal_c, Pattern = pal_nc_c)
      
      ## 9) Forzar niveles
      ann_col_c$Pattern <- factor(ann_col_c$Pattern, levels = names(ann_colors_c$Pattern))
      
      #### Heatmap 2 clusters de neuronas ----
      pheatmap(t(codes_c),
               color = hm_cols_c, breaks = breaks_c,
               border_color  = "grey60",
               scale         = "column",
               show_rownames = TRUE,
               show_colnames = FALSE,
               cluster_rows  = clust.sample_c,  cutree_rows = 2,
               annotation_row = ann_row_c,
               cluster_cols  = clust.neuron_c,  cutree_cols = k_neuron_c,
               annotation_col = ann_col_c,
               annotation_colors = ann_colors_c,
               legend_breaks = ticks_c,
               legend_labels = sprintf("%.1f", ticks_c),
               main = paste0("SOM codebooks — with neuron clusters (", tag_c, ")"))
      
      ## 11) Gene → cluster de neurona
      type.pattern_c <- clust.pat_c[classif_c]
      names(type.pattern_c) <- names(classif_c)
      cat("Genes per neuron cluster (", tag_c, "):\n", sep = ""); print(table(type.pattern_c))
      
      
      
      
      
      
      # PATHWAYS --------------------------
      
      
      # (B vs A, log2FC = B/A) --------------------------------------
      # Requisitos ya creados:
      # contr_1a == "A" (denominador), contr_2a == "B" (numerador)
      # res_a (DESeqResults), vst_a (DESeqTransform), coldata_a (niveles A/B)
      # enriched.FDR_a[["KEGG_Pathway"]], annotation (con columna "KEGG_Pathway")
      
      if (!exists("alpha_padj")) alpha_padj <- 0.05
      if (!exists("lfc_thr"))    lfc_thr    <- 1
      
      pal_grupo <- if (exists("pal_grupo")) pal_grupo else c(A="#fdac4f", B="#a2593d", Control="#3d6863")
      cols_dir_a <- c("down"="#aa4d6f", "up"="#a5923d")  # colores por dirección
      col_neg_a  <- "#aa4d6f";  col_pos_a <- "#a5923d"   # extremos heatmap
      
      # Filtros opcionales de rutas (deja NULL / character(0) si no quieres filtrar)
      keep_paths_a <- NULL
      drop_paths_a <- character(0)
      
      # --- 1) DEGs y FDR por pathway
      resdf_a     <- as.data.frame(res_a)
      deg_mask_a  <- !is.na(resdf_a$padj) & resdf_a$padj < alpha_padj
      deg_genes_a <- rownames(resdf_a)[deg_mask_a]
      
      fdr_path_a <- enriched.FDR_a[["KEGG_Pathway"]]
      stopifnot(!is.null(fdr_path_a))
      
      sig_paths_all_a <- names(fdr_path_a)[fdr_path_a < alpha_padj]
      sig_paths_all_a <- sub("^ko", "map", sig_paths_all_a)
      sig_paths_all_a <- sig_paths_all_a[!duplicated(sig_paths_all_a)]
      
      sig_paths_a <- sig_paths_all_a
      if (!is.null(keep_paths_a)) sig_paths_a <- intersect(sig_paths_a, keep_paths_a)
      if (length(drop_paths_a))   sig_paths_a <- setdiff(sig_paths_a, drop_paths_a)
      
      if (length(sig_paths_a) == 0) {
        message("[B_vs_A] No hay pathways enriquecidos (o fueron filtrados).")
      } else {
        
        # --- 2) ruta -> genes DEG de esa ruta
        path2genes_a <- setNames(vector("list", length(sig_paths_a)), sig_paths_a)
        for (p in sig_paths_a) {
          idx <- grep(p, annotation[,"KEGG_Pathway"], fixed = TRUE)
          g   <- intersect(rownames(annotation)[idx], deg_genes_a)
          path2genes_a[[p]] <- g
        }
        
        # --- 3) resumen por ruta + dirección por mediana(LFC)
        lfc_a <- resdf_a$log2FoldChange; names(lfc_a) <- rownames(resdf_a)
        
        path_summary_a <- do.call(rbind, lapply(names(path2genes_a), function(p){
          g <- path2genes_a[[p]]
          if (length(g) == 0) return(NULL)
          lfc_g <- lfc_a[g]
          data.frame(
            pathway    = p,
            FDR        = unname(fdr_path_a[p]),
            n_genes    = length(g),
            n_up       = sum(lfc_g >=  lfc_thr, na.rm=TRUE),
            n_down     = sum(lfc_g <= -lfc_thr, na.rm=TRUE),
            median_lfc = stats::median(lfc_g, na.rm=TRUE),
            dir_lab    = ifelse(stats::median(lfc_g, na.rm=TRUE) > 0, "up", "down"),
            stringsAsFactors = FALSE
          )
        }))
        rownames(path_summary_a) <- NULL
        
        ### Dotplot ---------------------------------
        path_summary_a$base_id <- sub("^ko", "map", path_summary_a$pathway)
        pth_a <- path_summary_a[!duplicated(path_summary_a$base_id), ]
        pth_a$mlogFDR <- -log10(pth_a$FDR)
        
        dotp_a <- ggplot(pth_a,
                         aes(x = mlogFDR, y = reorder(base_id, mlogFDR),
                             size = n_genes, color = dir_lab)) +
          geom_point() +
          scale_color_manual(limits = c("down","up"),
                             values = cols_dir_a,
                             name = "Dirección (median LFC)") +
          scale_size_area(max_size = 9, name = "DEGs en ruta") +
          labs(title = "Pathways enriquecidos (B vs A)",
               subtitle = sprintf("FDR < %.2f  •  color = dirección por mediana(log2FC)", alpha_padj),
               x = expression(-log[10]("FDR")), y = NULL) +
          theme_classic(base_size = 12) +
          theme(plot.title = element_text(hjust = .5, face = "bold"),
                plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30"))
        dotp_a
        
        # Heatmaps por pathway ----------------------------
        # (DEGs centrados por gen)
        ann_a <- data.frame(Group = coldata_a$Group, row.names = rownames(coldata_a))
        ann_a$Group <- factor(as.character(ann_a$Group), levels = c(contr_1a, contr_2a))  # A, B
        ann_colors_a <- list(Group = pal_grupo[levels(ann_a$Group)])
        
        nbreaks_a <- 256
        make_hm_a <- function(p){
          g <- path2genes_a[[p]]
          if (length(g) < 2) { message("Ruta con <2 genes DEG: ", p); return(invisible(NULL)) }
          m0 <- assay(vst_a)[g, , drop=FALSE]
          m  <- m0 - rowMeans(m0)               # centrar por gen (patrón)
          lim <- max(abs(m))
          brk <- seq(-lim, lim, length.out = nbreaks_a)
          cols <- colorRampPalette(c(col_neg_a, "#f7efe6", col_pos_a))(nbreaks_a - 1)
          
          pheatmap::pheatmap(m, color = cols, breaks = brk,
                             scale="none", border_color=NA,
                             clustering_distance_rows="euclidean",
                             clustering_distance_cols="correlation",
                             clustering_method="ward.D2",
                             show_rownames=FALSE, show_colnames=TRUE,
                             annotation_col=ann_a, annotation_colors=ann_colors_a,
                             main=paste0("Heatmap (DEGs) — ", p, "  [B vs A]"),
                             fontsize_col=10, treeheight_row=30, treeheight_col=30)
        }
        
        sig_to_plot_a <- pth_a$pathway
        invisible(lapply(sig_to_plot_a, make_hm_a))
      }
      
      
      
      # (A vs C, log2FC = A/C) -----------
      
      # Requisitos ya creados en tus bloques previos:
      # contr_1b == "C" (denominador), contr_2b == "A" (numerador)
      # res_b (DESeqResults), vst_b (DESeqTransform), coldata_b (Group con niveles C/A)
      # enriched.FDR_b (lista con "KEGG_Pathway"), annotation (con columna "KEGG_Pathway")
      
      # 0) Parámetros y helpers
      if (!exists("alpha_padj")) alpha_padj <- 0.05
      if (!exists("lfc_thr"))    lfc_thr    <- 1
      
      # Paletas (consistentes con lo previo)
      pal_grupo <- if (exists("pal_grupo")) pal_grupo else c(A="#fdac4f", B="#a2593d", Control="#3d6863")
      cols_dir_b <- c("down"="#aa4d6f", "up"="#a5923d")  # para dirección de ruta
      col_neg_b  <- "#aa4d6f";  col_pos_b <- "#a5923d"   # extremos del heatmap
      
      # Vectores de filtro opcional (si no quieres filtrar, déjalos como están)
      keep_paths_b <- NULL                 # p.ej.: c("Phenylpropanoid biosynthesis", "MAPK signaling pathway - plant")
      drop_paths_b <- character(0)         # p.ej.: c("Quorum sensing","Biofilm formation")
      
      # --- 1) Conjunto de DEGs y objetos base
      resdf_b <- as.data.frame(res_b)
      deg_mask_b <- !is.na(resdf_b$padj) & resdf_b$padj < alpha_padj
      deg_genes_b <- rownames(resdf_b)[deg_mask_b]
      
      # FDR por pathway (del bloque de enriquecimiento)
      fdr_path_b <- enriched.FDR_b[["KEGG_Pathway"]]
      stopifnot(!is.null(fdr_path_b))
      
      sig_paths_all_b <- names(fdr_path_b)[fdr_path_b < alpha_padj]
      sig_paths_all_b <- sub("^ko", "map", sig_paths_all_b)
      sig_paths_all_b <- sig_paths_all_b[!duplicated(sig_paths_all_b)]
      
      # Aplicar filtros opcionales (no filtra si keep=NULL y drop=vacío)
      sig_paths_b <- sig_paths_all_b
      
      if (!is.null(keep_paths_b)) sig_paths_b <- intersect(sig_paths_b, keep_paths_b)
      if (length(drop_paths_b))   sig_paths_b <- setdiff(sig_paths_b, drop_paths_b)
      
      #Path_summary_b
      if (length(sig_paths_b) == 0) {
        message(sprintf("[A_vs_C] No hay pathways enriquecidos (o fueron filtrados)."))
      } else {
        
        # --- 2) Mapeo ruta -> genes DEG de esa ruta
        # Usamos coincidencia por texto como en tu enriquecimiento (grep sobre annotation$KEGG_Pathway)
        path2genes_b <- setNames(vector("list", length(sig_paths_b)), sig_paths_b)
        for (p in sig_paths_b) {
          idx <- grep(p, annotation[,"KEGG_Pathway"], fixed = TRUE)
          g   <- intersect(rownames(annotation)[idx], deg_genes_b)
          path2genes_b[[p]] <- g
        }
        
        # --- 3) Resumen por ruta: dirección, tamaños, FDR
        lfc_b <- resdf_b$log2FoldChange; names(lfc_b) <- rownames(resdf_b)
        
        path_summary_b <- do.call(rbind, lapply(names(path2genes_b), function(p){
          g <- path2genes_b[[p]]
          if (length(g) == 0) return(NULL)
          lfc_g <- lfc_b[g]
          n_up   <- sum(lfc_g >=  lfc_thr, na.rm=TRUE)
          n_down <- sum(lfc_g <= -lfc_thr, na.rm=TRUE)
          data.frame(
            pathway    = p,
            FDR        = unname(fdr_path_b[p]),
            n_genes    = length(g),
            n_up       = n_up,
            n_down     = n_down,
            median_lfc = stats::median(lfc_g, na.rm=TRUE),
            dir_lab    = ifelse(stats::median(lfc_g, na.rm=TRUE) > 0, "up", "down"),
            stringsAsFactors = FALSE
          )
        }))
        rownames(path_summary_b) <- NULL
        
        # --- 4) 
        ### Dotplot ---------------------------
        # (enriquecimiento con dirección
        
        path_summary_b$base_id <- sub("^ko", "map", path_summary_b$pathway)  # añade esta columna a tu tabla de rutas
        pth_b_nodup   <- path_summary_b[!duplicated(path_summary_b$base_id), ]  # una por ruta
        
      
        
        pth_b_nodup$mlogFDR <- -log10(pth_b_nodup$FDR)
        ggplot(pth_b_nodup,
               aes(x = mlogFDR, y = reorder(base_id, mlogFDR),
                   size = n_genes, color = dir_lab)) +
          geom_point() +
          scale_color_manual(limits = c("down","up"),  values = cols_dir_b, name = "Dirección (median LFC)") +
          scale_size_area(max_size = 9, name = "DEGs en ruta") +
          labs(
            title    = "Pathways enriquecidos (A vs C)",
            subtitle = sprintf("FDR < %.2f  •  color = dirección por mediana(log2FC)", alpha_padj),
            x = expression(-log[10]("FDR")), y = NULL
          ) +
          theme_classic(base_size = 12) +
          theme(
            plot.title    = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30")
          )
        
        
        if (is.null(pth_b_nodup) || nrow(pth_b_nodup)==0) {
          message("[A_vs_C] No hay genes DEG asociados a las rutas significativas.")
        } else {
          
      
          
          ## Dotplot  ---------------------------
          # (enriquecimiento con dirección)
          # Solo es para ver si se imprime bien el fot plot, por eso está doble
          path_summary_b$base_id <- sub("^ko", "map", path_summary_b$pathway)  # añade esta columna a tu tabla de rutas
          pth_b_nodup   <- path_summary_b[!duplicated(path_summary_b$base_id), ]  # una por ruta
          
          pathways_keeper <- c("map00052", "map00400", "map00592", "map01110", "map01230")
          
          pth_b_filtered <- pth_b_nodup %>%
            filter(pathway %in% pathways_keeper)
          
          pth_b_filtered$mlogFDR <- -log10(pth_b_filtered$FDR)
          dotp_b <- (ggplot(pth_b_filtered,
                 aes(x = mlogFDR, y = reorder(base_id, mlogFDR),
                     size = n_genes, color = dir_lab)) +
            geom_point() +
            scale_color_manual(limits = c("down","up"),  values = cols_dir_b, name = "Dirección (median LFC)") +
            scale_size_area(max_size = 9, name = "DEGs en ruta") +
            labs(
              title    = "Pathways enriquecidos (A vs C)",
              subtitle = sprintf("FDR < %.2f  •  color = dirección por mediana(log2FC)", alpha_padj),
              x = expression(-log[10]("FDR")), y = NULL
            ) +
            theme_classic(base_size = 12) +
            theme(
              plot.title    = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30")
            ))
        
          dotp_b
          ### Heatmaps por pathway -----------------------------
          # (DEGs de cada ruta, centrados por gen
          # Anotación de columnas (A y C) con tus colores
          ann_b <- data.frame(Group = coldata_b$Group, row.names = rownames(coldata_b))
          ann_b$Group <- factor(as.character(ann_b$Group), levels = c(contr_1b, contr_2b))  # C, A
          ann_colors_b <- list(Group = pal_grupo[levels(ann_b$Group)])
          
          # Paleta continua simétrica alrededor de 0
          nbreaks_b <- 256
          # pequeña función para hacer cada heatmap:
          make_hm_b <- function(p){
            g <- path2genes_b[[p]]
            if (length(g) < 2) { message("Ruta con <2 genes DEG: ", p); return(invisible(NULL)) }
            mat0 <- assay(vst_b)[g, , drop = FALSE]
            # centrar por gen (patrón, no nivel)
            mat  <- mat0 - rowMeans(mat0)
            # límites simétricos para una escala comparable
            lim  <- max(abs(mat))
            breaks_b <- seq(-lim, lim, length.out = nbreaks_b)
            hm_cols_b <- colorRampPalette(c(col_neg_b, "#f7efe6", col_pos_b))(nbreaks_b - 1)
            # clustering por patrón
            pheatmap::pheatmap(
              mat,
              color = hm_cols_b, breaks = breaks_b,
              scale = "none", border_color = NA,
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "correlation",
              clustering_method = "ward.D2",
              show_rownames = FALSE, show_colnames = TRUE,
              annotation_col = ann_b, annotation_colors = ann_colors_b,
              main = paste0("Heatmap (DEGs) — ", p, "  [A vs C]"),
              fontsize_col = 10, treeheight_row = 30, treeheight_col = 30
            )
          }
          
          # Si quieres limitar a las Top-N por evidencia:
          # pth_b_nodup <- pth_b_nodup[order(pth_b_nodup$FDR), ]
          # sig_to_plot <- head(pth_b_nodup$pathway, 8)
          sig_to_plot <- pth_b_nodup$pathway
          
          invisible(lapply(sig_to_plot, make_hm_b))
        } 
      } 
      
      dotp_b
      
      
      
      #  (B vs C, log2FC = B/C) ------------------------------
      # Requisitos ya creados:
      # contr_1c == "C" (denominador), contr_2c == "B" (numerador)
      # res_c (DESeqResults), vst_c (DESeqTransform), coldata_c (niveles C/B)
      # enriched.FDR_c[["KEGG_Pathway"]], annotation (con columna "KEGG_Pathway")
      
      if (!exists("alpha_padj")) alpha_padj <- 0.05
      if (!exists("lfc_thr"))    lfc_thr    <- 1
      
      pal_grupo <- if (exists("pal_grupo")) pal_grupo else c(A="#fdac4f", B="#a2593d", Control="#3d6863")
      cols_dir_c <- c("down"="#aa4d6f", "up"="#a5923d")
      col_neg_c  <- "#aa4d6f";  col_pos_c <- "#a5923d"
      
      keep_paths_c <- NULL
      drop_paths_c <- character(0)
      
      # --- 1) DEGs y FDR por pathway
      resdf_c     <- as.data.frame(res_c)
      deg_mask_c  <- !is.na(resdf_c$padj) & resdf_c$padj < alpha_padj
      deg_genes_c <- rownames(resdf_c)[deg_mask_c]
      
      fdr_path_c <- enriched.FDR_c[["KEGG_Pathway"]]
      stopifnot(!is.null(fdr_path_c))
      
      sig_paths_all_c <- names(fdr_path_c)[fdr_path_c < alpha_padj]
      sig_paths_all_c <- sub("^ko", "map", sig_paths_all_c)
      sig_paths_all_c <- sig_paths_all_c[!duplicated(sig_paths_all_c)]
      
      sig_paths_c <- sig_paths_all_c
      if (!is.null(keep_paths_c)) sig_paths_c <- intersect(sig_paths_c, keep_paths_c)
      if (length(drop_paths_c))   sig_paths_c <- setdiff(sig_paths_c, drop_paths_c)
      
      if (length(sig_paths_c) == 0) {
        message("[B_vs_C] No hay pathways enriquecidos (o fueron filtrados).")
      } else {
        
        # --- 2) ruta -> genes DEG
        path2genes_c <- setNames(vector("list", length(sig_paths_c)), sig_paths_c)
        for (p in sig_paths_c) {
          idx <- grep(p, annotation[,"KEGG_Pathway"], fixed = TRUE)
          g   <- intersect(rownames(annotation)[idx], deg_genes_c)
          path2genes_c[[p]] <- g
        }
        
        # --- 3) resumen + dirección por mediana(LFC)
        lfc_c <- resdf_c$log2FoldChange; names(lfc_c) <- rownames(resdf_c)
        
        path_summary_c <- do.call(rbind, lapply(names(path2genes_c), function(p){
          g <- path2genes_c[[p]]
          if (length(g) == 0) return(NULL)
          lfc_g <- lfc_c[g]
          data.frame(
            pathway    = p,
            FDR        = unname(fdr_path_c[p]),
            n_genes    = length(g),
            n_up       = sum(lfc_g >=  lfc_thr, na.rm=TRUE),
            n_down     = sum(lfc_g <= -lfc_thr, na.rm=TRUE),
            median_lfc = stats::median(lfc_g, na.rm=TRUE),
            dir_lab    = ifelse(stats::median(lfc_g, na.rm=TRUE) > 0, "up", "down"),
            stringsAsFactors = FALSE
          )
        }))
        rownames(path_summary_c) <- NULL
        
        #### Dotplot ----------------------------------------------
        path_summary_c$base_id <- sub("^ko", "map", path_summary_c$pathway)
        pth_c <- path_summary_c[!duplicated(path_summary_c$base_id), ]
        pth_c$mlogFDR <- -log10(pth_c$FDR)
        
        dotp_c <- ggplot(pth_c,
                         aes(x = mlogFDR, y = reorder(base_id, mlogFDR),
                             size = n_genes, color = dir_lab)) +
          geom_point() +
          scale_color_manual(limits = c("down","up"),
                             values = cols_dir_c,
                             name = "Dirección (median LFC)") +
          scale_size_area(max_size = 9, name = "DEGs en ruta") +
          labs(title = "Pathways enriquecidos (B vs C)",
               subtitle = sprintf("FDR < %.2f  •  color = dirección por mediana(log2FC)", alpha_padj),
               x = expression(-log[10]("FDR")), y = NULL) +
          theme_classic(base_size = 12) +
          theme(plot.title = element_text(hjust = .5, face = "bold"),
                plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30"))
        dotp_c
        
        #### Heatmaps por pathway ----------------------------------------
        ann_c <- data.frame(Group = coldata_c$Group, row.names = rownames(coldata_c))
        ann_c$Group <- factor(as.character(ann_c$Group), levels = c(contr_1c, contr_2c))  # C, B
        ann_colors_c <- list(Group = pal_grupo[levels(ann_c$Group)])
        
        nbreaks_c <- 256
        make_hm_c <- function(p){
          g <- path2genes_c[[p]]
          if (length(g) < 2) { message("Ruta con <2 genes DEG: ", p); return(invisible(NULL)) }
          m0 <- assay(vst_c)[g, , drop=FALSE]
          m  <- m0 - rowMeans(m0)
          lim <- max(abs(m))
          brk <- seq(-lim, lim, length.out = nbreaks_c)
          cols <- colorRampPalette(c(col_neg_c, "#f7efe6", col_pos_c))(nbreaks_c - 1)
          
          pheatmap::pheatmap(m, color = cols, breaks = brk,
                             scale="none", border_color=NA,
                             clustering_distance_rows="euclidean",
                             clustering_distance_cols="correlation",
                             clustering_method="ward.D2",
                             show_rownames=FALSE, show_colnames=TRUE,
                             annotation_col=ann_c, annotation_colors=ann_colors_c,
                             main=paste0("Heatmap (DEGs) — ", p, "  [B vs C]"),
                             fontsize_col=10, treeheight_row=30, treeheight_col=30)
        }
        
        sig_to_plot_c <- pth_c$pathway
        invisible(lapply(sig_to_plot_c, make_hm_c))
      }
      
      
      
      
      
      
      
      
      
      
      
      # PLOTING B/A — Heatmap top-N DEGs + Universos + Anotación -----------------
      
      ### Heatmap top-N DEGs ----------------------------
      
      ## --- Parámetros del contraste B/A ---
      if (!exists("alpha_padj")) alpha_padj <- 0.05
      if (!exists("lfc_thr"))    lfc_thr    <- 1
      topN_a <- 300     # <- cambia aquí cuántos DEGs (más variables) quieres en el heatmap
      
      ## --- 0) Objetos base (DEGs por FDR) ---
      resdf_a <- as.data.frame(res_a)
      deg_mask_a <- !is.na(resdf_a$padj) & resdf_a$padj < alpha_padj
      deg_ids_a  <- rownames(resdf_a)[deg_mask_a]
      cat("[B/A] #DEGs (padj <", alpha_padj, "):", length(deg_ids_a), "\n")
      
      ## --- 1) Heatmap global de DEGs más variables (B/A) ---
      # Qué se grafica: SOLO DEGs del contraste B/A; seleccionamos los top-N por varianza tras VST.
      vst_mat_a <- SummarizedExperiment::assay(vst_a)     # genes x muestras
      mat_deg_a <- vst_mat_a[deg_ids_a, , drop = FALSE]   # solo DEGs
      
      # si no hay DEGs, salimos temprano
      if (nrow(mat_deg_a) == 0) {
        stop("[B/A] No hay DEGs con el umbral actual; ajusta alpha_padj o revisa el contraste.")
      }
      
      # top-N por varianza entre DEGs
      library(matrixStats)
      vars_a <- rowVars(mat_deg_a)
      ord_a  <- order(vars_a, decreasing = TRUE)
      sel_a  <- head(rownames(mat_deg_a)[ord_a], min(topN_a, length(ord_a)))
      
      # centrar por gen (patrón, no nivel)
      mat_hm_a0 <- mat_deg_a[sel_a, , drop = FALSE]
      mat_hm_a  <- mat_hm_a0 - rowMeans(mat_hm_a0)
      
      # anotación de columnas (A/B)
      ann_a <- data.frame(Group = coldata_a$Group, row.names = rownames(coldata_a))
      ann_a$Group <- droplevels(ann_a$Group)
      
      # paleta de grupos
      if (!exists("pal_grupo")) pal_grupo <- c(A="#fdac4f", B="#a2593d", Control="#3d6863")
      ann_colors_a <- list(Group = pal_grupo[levels(ann_a$Group)])
      
      # paleta continua simétrica (centrada en 0)
      nbreaks_a <- 256
      lim_a     <- max(abs(mat_hm_a))
      breaks_a  <- seq(-lim_a, lim_a, length.out = nbreaks_a)
      col_neg_a <- "#aa4d6f";  col_pos_a <- "#a5923d"
      hm_cols_a <- colorRampPalette(c(col_neg_a, "#f7efe6", col_pos_a))(nbreaks_a - 1)
      ticks_a   <- pretty(c(-lim_a, lim_a), n = 5)
      
      # plot
      pheatmap::pheatmap(
        mat_hm_a,
        color = hm_cols_a, breaks = breaks_a,
        scale = "none", border_color = NA,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "correlation",
        clustering_method = "ward.D2",
        show_rownames = FALSE, show_colnames = TRUE,
        annotation_col = ann_a, annotation_colors = ann_colors_a,
        legend_breaks = ticks_a, legend_labels = sprintf("%.1f", ticks_a),
        main = sprintf("B vs A — Heatmap top-%d DEGs (centrado por gen)", length(sel_a)),
        fontsize_col = 10, treeheight_row = 30, treeheight_col = 30
      )
      
      
      
      ## Universos + export (B/A) — sufijo a
      ## Requiere: res_a, counts_a, dds_a, coldata_a, annotation,
      ##           contr_1a (A), contr_2a (B)
      
      
      # Umbrales (si no existen, definir)
      if (!exists("alpha_padj")) alpha_padj <- 0.05
      if (!exists("lfc_thr"))    lfc_thr    <- 1
      
      # Data frame base de resultados
      resdf_a <- as.data.frame(res_a)
      resdf_a$gene <- rownames(resdf_a)
      
      # Máscara de significancia y dirección (respecto a contr_2a = B)
      sig_mask_a <- !is.na(resdf_a$padj) & (resdf_a$padj < alpha_padj)
      resdf_a$direction_raw <- ifelse(resdf_a$log2FoldChange >  0, "up", 
                                      ifelse(resdf_a$log2FoldChange <  0, "down", "no_change"))
      resdf_a$direction <- ifelse(sig_mask_a & abs(resdf_a$log2FoldChange) >= lfc_thr,
                                  ifelse(resdf_a$log2FoldChange > 0,  paste0("up in ", contr_2a),
                                         paste0("down in ", contr_2a)),
                                  "ns")
      
      ### Export DEGs ---------------------------
      
      res_sig_a <- subset(resdf_a, sig_mask_a)
      # Columna binaria up/down (solo para sig); si quieres mantener “ns”, deja direction completo
      res_sig_a$up_down <- ifelse(res_sig_a$log2FoldChange > 0, "up", "down")
      
      out_deg_csv_a <- paste0("DEG_gene_results_", contr_2a, "_vs_", contr_1a, ".csv")
      write.csv(res_sig_a[, c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","up_down")],
                file = out_deg_csv_a, row.names = FALSE)
      cat(">> Guardado:", out_deg_csv_a, "\n")
      
      
      ### Universos (B/A) ------------------------------------------------------
      
      # Detectados = genes en la matriz del contraste
      detected_a <- nrow(counts_a)
      
      # Testeados = genes con un resultado estadístico
      tested_ids_a <- rownames(resdf_a)
      tested_a     <- length(tested_ids_a)
      
      # DEGs totales y por sentido (con |log2FC| >= lfc_thr)
      deg_ids_a     <- resdf_a$gene[sig_mask_a & abs(resdf_a$log2FoldChange) >= lfc_thr]
      deg_up_ids_a  <- resdf_a$gene[sig_mask_a & resdf_a$log2FoldChange >=  lfc_thr]
      deg_dn_ids_a  <- resdf_a$gene[sig_mask_a & resdf_a$log2FoldChange <= -lfc_thr]
      
      n_deg_a   <- length(deg_ids_a)
      n_up_a    <- length(deg_up_ids_a)
      n_down_a  <- length(deg_dn_ids_a)
      
      # Anotación KEGG (usamos columna KEGG_Pathway como proxy de “anotado en KEGG”)
      # (ajusta si prefieres KEGG_Module o considerar cualquiera de las dos)
      anno_ids_a <- intersect(tested_ids_a, rownames(annotation))
      kegg_ok_tested <- with(annotation[anno_ids_a, , drop=FALSE],
                             !is.na(KEGG_Pathway) & KEGG_Pathway != "")
      n_kegg_tested  <- sum(kegg_ok_tested)
      
      # Entre los DEGs (up/down)
      anno_up_ids_a   <- intersect(deg_up_ids_a, rownames(annotation))
      anno_down_ids_a <- intersect(deg_dn_ids_a, rownames(annotation))
      
      n_kegg_up_a   <- sum(!is.na(annotation[anno_up_ids_a,  "KEGG_Pathway"])   &
                             annotation[anno_up_ids_a,  "KEGG_Pathway"]   != "")
      n_kegg_down_a <- sum(!is.na(annotation[anno_down_ids_a,"KEGG_Pathway"])   &
                             annotation[anno_down_ids_a,"KEGG_Pathway"]   != "")
      
      # Armar tabla “universo”
      universe_a <- data.frame(
        contrast          = paste0(contr_2a, "_vs_", contr_1a),  # "B_vs_A"
        detected_genes    = detected_a,
        tested_genes      = tested_a,
        DEGs_total        = n_deg_a,
        DEGs_up           = n_up_a,
        DEGs_down         = n_down_a,
        KEGG_annot_tested = n_kegg_tested,
        KEGG_annot_up     = n_kegg_up_a,
        KEGG_annot_down   = n_kegg_down_a,
        alpha_padj        = alpha_padj,
        lfc_threshold     = lfc_thr,
        stringsAsFactors  = FALSE
      )
      
      out_universe_csv_a <- paste0("Universe_", contr_2a, "_vs_", contr_1a, ".csv")
      write.csv(universe_a, file = out_universe_csv_a, row.names = FALSE)
      cat(">> Guardado:", out_universe_csv_a, "\n")
      
      # (Opcional) imprime resumen en consola
      print(universe_a)
      
      
      
      
      # PLOTING A/C — Heatmap top-N DEGs + Universos + Anotación --------------------
      
      ## Parámetros
      if (!exists("alpha_padj")) alpha_padj <- 0.05
      if (!exists("lfc_thr"))    lfc_thr    <- 1
      topN_b <- 300
      
      ## 0) Objetos base
      resdf_b   <- as.data.frame(res_b)
      deg_mask_b<- !is.na(resdf_b$padj) & resdf_b$padj < alpha_padj
      deg_ids_b <- rownames(resdf_b)[deg_mask_b]
      cat("[A/C] #DEGs (padj <", alpha_padj, "):", length(deg_ids_b), "\n")
      
      ### Heatmap top-N (solo DEGs) ------------------------
      vst_mat_b <- SummarizedExperiment::assay(vst_b)
      mat_deg_b <- vst_mat_b[deg_ids_b, , drop = FALSE]
      if (nrow(mat_deg_b) == 0) stop("[A/C] No hay DEGs con el umbral actual.")
      
      library(matrixStats)
      vars_b <- rowVars(mat_deg_b)
      ord_b  <- order(vars_b, decreasing = TRUE)
      sel_b  <- head(rownames(mat_deg_b)[ord_b], min(topN_b, length(ord_b)))
      
      mat_hm_b0 <- mat_deg_b[sel_b, , drop = FALSE]
      mat_hm_b  <- mat_hm_b0 - rowMeans(mat_hm_b0)
      
      ann_b <- data.frame(Group = coldata_b$Group, row.names = rownames(coldata_b))
      ann_b$Group <- droplevels(ann_b$Group)
      
      if (!exists("pal_grupo")) pal_grupo <- c(A="#fdac4f", B="#a2593d", Control="#3d6863")
      ann_colors_b <- list(Group = pal_grupo[levels(ann_b$Group)])
      
      nbreaks_b <- 256
      lim_b     <- max(abs(mat_hm_b)) + 1e-8
      breaks_b  <- seq(-lim_b, lim_b, length.out = nbreaks_b)
      col_neg_b <- "#aa4d6f";  col_pos_b <- "#a5923d"
      hm_cols_b <- colorRampPalette(c(col_neg_b, "#f7efe6", col_pos_b))(nbreaks_b - 1)
      ticks_b   <- pretty(c(-lim_b, lim_b), n = 5)
      
      pheatmap::pheatmap(
        mat_hm_b,
        color = hm_cols_b, breaks = breaks_b,
        scale = "none", border_color = NA,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "correlation",
        clustering_method = "ward.D2",
        show_rownames = FALSE, show_colnames = TRUE,
        annotation_col = ann_b, annotation_colors = ann_colors_b,
        legend_breaks = ticks_b, legend_labels = sprintf("%.1f", ticks_b),
        main = sprintf("A vs C — Heatmap top-%d DEGs (centrado por gen)", length(sel_b)),
        fontsize_col = 10, treeheight_row = 30, treeheight_col = 30
      )
      
      ## Universos + export (A/C) — sufijo b
      if (!exists("annotation")) annotation <- read.delim("fullAnnotation.tsv.txt", stringsAsFactors = FALSE, row.names = 1)
      
      resdf_b$gene <- rownames(resdf_b)
      sig_mask_b <- !is.na(resdf_b$padj) & (resdf_b$padj < alpha_padj)
      resdf_b$direction_raw <- ifelse(resdf_b$log2FoldChange > 0, "up",
                                      ifelse(resdf_b$log2FoldChange < 0, "down", "no_change"))
      resdf_b$direction <- ifelse(sig_mask_b & abs(resdf_b$log2FoldChange) >= lfc_thr,
                                  ifelse(resdf_b$log2FoldChange > 0,  paste0("up in ", contr_2b),
                                         paste0("down in ", contr_2b)),
                                  "ns")
      
      res_sig_b <- subset(resdf_b, sig_mask_b)
      res_sig_b$up_down <- ifelse(res_sig_b$log2FoldChange > 0, "up", "down")
      
      ### Export DEGs ---------------------------
      out_deg_csv_b <- paste0("DEG_gene_results_", contr_2b, "_vs_", contr_1b, ".csv")
      write.csv(res_sig_b[, c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","up_down")],
                file = out_deg_csv_b, row.names = FALSE)
      cat(">> Guardado:", out_deg_csv_b, "\n")
      
      detected_b <- nrow(counts_b)
      tested_ids_b <- rownames(resdf_b)
      tested_b     <- length(tested_ids_b)
      
      deg_ids_b     <- resdf_b$gene[sig_mask_b & abs(resdf_b$log2FoldChange) >= lfc_thr]
      deg_up_ids_b  <- resdf_b$gene[sig_mask_b & resdf_b$log2FoldChange >=  lfc_thr]
      deg_dn_ids_b  <- resdf_b$gene[sig_mask_b & resdf_b$log2FoldChange <= -lfc_thr]
      
      n_deg_b  <- length(deg_ids_b)
      n_up_b   <- length(deg_up_ids_b)
      n_down_b <- length(deg_dn_ids_b)
      
      anno_ids_b <- intersect(tested_ids_b, rownames(annotation))
      kegg_ok_tested_b <- with(annotation[anno_ids_b, , drop=FALSE],
                               !is.na(KEGG_Pathway) & KEGG_Pathway != "")
      n_kegg_tested_b  <- sum(kegg_ok_tested_b)
      
      anno_up_ids_b   <- intersect(deg_up_ids_b, rownames(annotation))
      anno_down_ids_b <- intersect(deg_dn_ids_b, rownames(annotation))
      n_kegg_up_b   <- sum(!is.na(annotation[anno_up_ids_b,  "KEGG_Pathway"]) &
                             annotation[anno_up_ids_b,  "KEGG_Pathway"] != "")
      n_kegg_down_b <- sum(!is.na(annotation[anno_down_ids_b,"KEGG_Pathway"]) &
                             annotation[anno_down_ids_b,"KEGG_Pathway"] != "")
      
      ### Universos (C/A) ------------------------
      universe_b <- data.frame(
        contrast          = paste0(contr_2b, "_vs_", contr_1b),  # "A_vs_C"
        detected_genes    = detected_b,
        tested_genes      = tested_b,
        DEGs_total        = n_deg_b,
        DEGs_up           = n_up_b,
        DEGs_down         = n_down_b,
        KEGG_annot_tested = n_kegg_tested_b,
        KEGG_annot_up     = n_kegg_up_b,
        KEGG_annot_down   = n_kegg_down_b,
        alpha_padj        = alpha_padj,
        lfc_threshold     = lfc_thr,
        stringsAsFactors  = FALSE
      )
      
      out_universe_csv_b <- paste0("Universe_", contr_2b, "_vs_", contr_1b, ".csv")
      write.csv(universe_b, file = out_universe_csv_b, row.names = FALSE)
      cat(">> Guardado:", out_universe_csv_b, "\n")
      print(universe_b)
      
      
      
      
      # PLOTING B/C — Heatmap top-N DEGs + Universos + Anotación -----------------------
      
      ## Parámetros
      if (!exists("alpha_padj")) alpha_padj <- 0.05
      if (!exists("lfc_thr"))    lfc_thr    <- 1
      topN_c <- 300
      
      ## 0) Objetos base
      resdf_c   <- as.data.frame(res_c)
      deg_mask_c<- !is.na(resdf_c$padj) & resdf_c$padj < alpha_padj
      deg_ids_c <- rownames(resdf_c)[deg_mask_c]
      cat("[B/C] #DEGs (padj <", alpha_padj, "):", length(deg_ids_c), "\n")
      
      ## 1) Heatmap top-N (solo DEGs)
      vst_mat_c <- SummarizedExperiment::assay(vst_c)
      mat_deg_c <- vst_mat_c[deg_ids_c, , drop = FALSE]
      if (nrow(mat_deg_c) == 0) stop("[B/C] No hay DEGs con el umbral actual.")
      
      library(matrixStats)
      vars_c <- rowVars(mat_deg_c)
      ord_c  <- order(vars_c, decreasing = TRUE)
      sel_c  <- head(rownames(mat_deg_c)[ord_c], min(topN_c, length(ord_c)))
      
      mat_hm_c0 <- mat_deg_c[sel_c, , drop = FALSE]
      mat_hm_c  <- mat_hm_c0 - rowMeans(mat_hm_c0)
      
      ann_c <- data.frame(Group = coldata_c$Group, row.names = rownames(coldata_c))
      ann_c$Group <- droplevels(ann_c$Group)
      
      if (!exists("pal_grupo")) pal_grupo <- c(A="#fdac4f", B="#a2593d", Control="#3d6863")
      ann_colors_c <- list(Group = pal_grupo[levels(ann_c$Group)])
      
      nbreaks_c <- 256
      lim_c     <- max(abs(mat_hm_c)) + 1e-8
      breaks_c  <- seq(-lim_c, lim_c, length.out = nbreaks_c)
      col_neg_c <- "#aa4d6f";  col_pos_c <- "#a5923d"
      hm_cols_c <- colorRampPalette(c(col_neg_c, "#f7efe6", col_pos_c))(nbreaks_c - 1)
      ticks_c   <- pretty(c(-lim_c, lim_c), n = 5)
      
      ### Heatmap top-N (solo DEGs) -------------------------------------------
      pheatmap::pheatmap(
        mat_hm_c,
        color = hm_cols_c, breaks = breaks_c,
        scale = "none", border_color = NA,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "correlation",
        clustering_method = "ward.D2",
        show_rownames = FALSE, show_colnames = TRUE,
        annotation_col = ann_c, annotation_colors = ann_colors_c,
        legend_breaks = ticks_c, legend_labels = sprintf("%.1f", ticks_c),
        main = sprintf("B vs C — Heatmap top-%d DEGs (centrado por gen)", length(sel_c)),
        fontsize_col = 10, treeheight_row = 30, treeheight_col = 30
      )
      
      ## Universos + export (B/C) — sufijo c
      if (!exists("annotation")) annotation <- read.delim("fullAnnotation.tsv.txt", stringsAsFactors = FALSE, row.names = 1)
      
      resdf_c$gene <- rownames(resdf_c)
      sig_mask_c <- !is.na(resdf_c$padj) & (resdf_c$padj < alpha_padj)
      resdf_c$direction_raw <- ifelse(resdf_c$log2FoldChange > 0, "up",
                                      ifelse(resdf_c$log2FoldChange < 0, "down", "no_change"))
      resdf_c$direction <- ifelse(sig_mask_c & abs(resdf_c$log2FoldChange) >= lfc_thr,
                                  ifelse(resdf_c$log2FoldChange > 0,  paste0("up in ", contr_2c),
                                         paste0("down in ", contr_2c)),
                                  "ns")
      ### Export DEGs ------------------------------
      res_sig_c <- subset(resdf_c, sig_mask_c)
      res_sig_c$up_down <- ifelse(res_sig_c$log2FoldChange > 0, "up", "down")
      
      out_deg_csv_c <- paste0("DEG_gene_results_", contr_2c, "_vs_", contr_1c, ".csv")
      write.csv(res_sig_c[, c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","up_down")],
                file = out_deg_csv_c, row.names = FALSE)
      cat(">> Guardado:", out_deg_csv_c, "\n")
      
      detected_c <- nrow(counts_c)
      tested_ids_c <- rownames(resdf_c)
      tested_c     <- length(tested_ids_c)
      
      deg_ids_c     <- resdf_c$gene[sig_mask_c & abs(resdf_c$log2FoldChange) >= lfc_thr]
      deg_up_ids_c  <- resdf_c$gene[sig_mask_c & resdf_c$log2FoldChange >=  lfc_thr]
      deg_dn_ids_c  <- resdf_c$gene[sig_mask_c & resdf_c$log2FoldChange <= -lfc_thr]
      
      n_deg_c  <- length(deg_ids_c)
      n_up_c   <- length(deg_up_ids_c)
      n_down_c <- length(deg_dn_ids_c)
      
      anno_ids_c <- intersect(tested_ids_c, rownames(annotation))
      kegg_ok_tested_c <- with(annotation[anno_ids_c, , drop=FALSE],
                               !is.na(KEGG_Pathway) & KEGG_Pathway != "")
      n_kegg_tested_c  <- sum(kegg_ok_tested_c)
      
      anno_up_ids_c   <- intersect(deg_up_ids_c, rownames(annotation))
      anno_down_ids_c <- intersect(deg_dn_ids_c, rownames(annotation))
      n_kegg_up_c   <- sum(!is.na(annotation[anno_up_ids_c,  "KEGG_Pathway"]) &
                             annotation[anno_up_ids_c,  "KEGG_Pathway"] != "")
      n_kegg_down_c <- sum(!is.na(annotation[anno_down_ids_c,"KEGG_Pathway"]) &
                             annotation[anno_down_ids_c,"KEGG_Pathway"] != "")
      ### Universos (C/A) ----------------------------------------
      universe_c <- data.frame(
        contrast          = paste0(contr_2c, "_vs_", contr_1c),  # "B_vs_C"
        detected_genes    = detected_c,
        tested_genes      = tested_c,
        DEGs_total        = n_deg_c,
        DEGs_up           = n_up_c,
        DEGs_down         = n_down_c,
        KEGG_annot_tested = n_kegg_tested_c,
        KEGG_annot_up     = n_kegg_up_c,
        KEGG_annot_down   = n_kegg_down_c,
        alpha_padj        = alpha_padj,
        lfc_threshold     = lfc_thr,
        stringsAsFactors  = FALSE
      )
      
      out_universe_csv_c <- paste0("Universe_", contr_2c, "_vs_", contr_1c, ".csv")
      write.csv(universe_c, file = out_universe_csv_c, row.names = FALSE)
      cat(">> Guardado:", out_universe_csv_c, "\n")
      print(universe_c)
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      # VENN (tres contrastes) – UP y DOWN -----------------------------------
      
      
      # --- Umbrales por si no existen 
      if (!exists("alpha_padj")) alpha_padj <- 0.05
      if (!exists("lfc_thr"))    lfc_thr    <- 1
      
      # helper: IDs UP/DOWN a partir de un DESeqResults 
      get_up_down_ids <- function(res, alpha = alpha_padj, lfc = lfc_thr) {
        df <- as.data.frame(res)
        df$gene <- rownames(df)
        df <- df[!is.na(df$padj), , drop = FALSE]
        up_ids   <- df$gene[df$padj < alpha & df$log2FoldChange >=  lfc]
        down_ids <- df$gene[df$padj < alpha & df$log2FoldChange <= -lfc]
        list(up = unique(up_ids), down = unique(down_ids))
      }
      
      # helper: tabla de cuentas por región (A solo, B solo, C solo, etc.)
      venn_region_counts <- function(sets_named) {
        stopifnot(length(sets_named) == 3)
        A <- sets_named[[1]]; B <- sets_named[[2]]; C <- sets_named[[3]]
        Aonly <- setdiff(A, union(B, C))
        Bonly <- setdiff(B, union(A, C))
        Conly <- setdiff(C, union(A, B))
        AB    <- setdiff(intersect(A, B), C)
        AC    <- setdiff(intersect(A, C), B)
        BC    <- setdiff(intersect(B, C), A)
        ABC   <- Reduce(intersect, list(A, B, C))
        data.frame(
          region = c(names(sets_named)[1], names(sets_named)[2], names(sets_named)[3],
                     paste(names(sets_named)[1:2], collapse = "∩"),
                     paste(names(sets_named)[c(1,3)], collapse = "∩"),
                     paste(names(sets_named)[2:3], collapse = "∩"),
                     "A∩B∩C"),
          n = c(length(Aonly), length(Bonly), length(Conly),
                length(AB), length(AC), length(BC), length(ABC)),
          stringsAsFactors = FALSE
        )
      }
      
      ### Plot ---------
      #  plot + guardar (Opción A: outfile requerido) 
      plot_venn_3 <- function(sets, title, subtitle,
                              lab_a, lab_b, lab_c,
                              fill_cols, outfile,
                              width = 7, height = 5, dpi = 300) {
        
        stopifnot(length(sets) == 3)
        names(sets) <- c(lab_a, lab_b, lab_c)
        
        p <- ggvenn(sets,
                    fill_color = fill_cols,
                    show_percentage = FALSE, stroke_size = 0.6) +
          labs(title = title, subtitle = subtitle) +
          theme_classic(base_size = 12) +
          theme(
            plot.title    = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 1, size = 9, colour = "grey30")
          )
        
        ggsave(outfile, p, width = width, height = height, dpi = dpi)
        p
      }
      
      # etiquetas legibles por contraste 
      lab_a <- paste0(contr_2a, "_vs_", contr_1a)  # B_vs_A
      lab_b <- paste0(contr_2b, "_vs_", contr_1b)  # A_vs_C
      lab_c <- paste0(contr_2c, "_vs_", contr_1c)  # B_vs_C
      
      #  paleta para los tres juegos (contrastes) 
      # orden: a, b, c
      venn_cols <- c("#a2593d", "#fdac4f", "#3d6863")  # B_vs_A, A_vs_C, B_vs_C
      
      # obtener sets por contraste 
      ud_a <- get_up_down_ids(res_a)  # B/A
      ud_b <- get_up_down_ids(res_b)  # A/C
      ud_c <- get_up_down_ids(res_c)  # B/C
      
      sets_up <- list(ud_a$up,   ud_b$up,   ud_c$up)
      sets_down <- list(ud_a$down, ud_b$down, ud_c$down)
      
      ### exportar tablas de regiones ----------
      cnt_up   <- venn_region_counts(setNames(sets_up,   c(lab_a, lab_b, lab_c)))
      cnt_down <- venn_region_counts(setNames(sets_down, c(lab_a, lab_b, lab_c)))
      write.csv(cnt_up,   file = "Venn_regions_counts_UP.csv",   row.names = FALSE)
      write.csv(cnt_down, file = "Venn_regions_counts_DOWN.csv", row.names = FALSE)
      cat(">> Guardados: Venn_regions_counts_UP.csv, Venn_regions_counts_DOWN.csv\n")
      
      ### plot printing ----------
      subtxt <- sprintf("padj < %.2f  •  (log2FC ≥ %g)", alpha_padj, lfc_thr)
      
      venn_up_plot <- plot_venn_3(
        sets      = setNames(sets_up, c(lab_a, lab_b, lab_c)),
        title     = "Venn — UP (B/A, A/C, B/C)",
        subtitle  = subtxt,
        lab_a     = lab_a, lab_b = lab_b, lab_c = lab_c,
        fill_cols = venn_cols,
        outfile   = "Venn_UP.png"
      )
      print(venn_up_plot)
      cat(">> Guardado: Venn_UP.png\n")
      
      venn_down_plot <- plot_venn_3(
        sets      = setNames(sets_down, c(lab_a, lab_b, lab_c)),
        title     = "Venn — DOWN (B/A, A/C, B/C)",
        subtitle  = subtxt,
        lab_a     = lab_a, lab_b = lab_b, lab_c = lab_c,
        fill_cols = venn_cols,
        outfile   = "Venn_DOWN.png"
      )
      print(venn_down_plot)
      cat(">> Guardado: Venn_DOWN.png\n")
      
      
      
      
      
      
      
      
      
      
      
      
      
      # DEGs avec Pfams & log2FC ---------
      
      ##  DEG + Anotación (Pfams/Description/GOs/KEGG) por contraste 
      
      # 0) Cargar anotación si hace falta
      if (!exists("annotation")) {
        annotation <- read.delim("fullAnnotation.tsv.txt",
                                 stringsAsFactors = FALSE,
                                 row.names = 1, check.names = FALSE)
      }
      
      # 1) Helper para unir y exportar (mantiene el orden de res_sig_*)
      annotate_and_export <- function(res_sig, contr_2, contr_1, suffix) {
        stopifnot(is.data.frame(res_sig), "gene" %in% colnames(res_sig))
        # columnas útiles si existen en tu anotación
        ann_cols <- intersect(
          c("Preferred_name","Description","PFAMs","PFAM","GOs","GO",
            "KEGG_ko","KEGG_Pathway","KEGG_Module"),
          colnames(annotation)
        )
        idx <- match(res_sig$gene, rownames(annotation))      # alinea por ID de gen
        ann_part <- annotation[idx, ann_cols, drop = FALSE]   # puede traer NAs si falta
        out <- cbind(res_sig, ann_part)
        
        
        # crea un objeto con sufijo (a/b/c) en el GlobalEnv y exporta CSV
        assign(paste0("DEG_annotated_", suffix), out, envir = .GlobalEnv)
        out_file <- paste0("DEG_annotated_", contr_2, "_vs_", contr_1, ".csv")
        write.csv(out, out_file, row.names = FALSE)
        cat(">> Guardado:", out_file, "| rows:", nrow(out), "\n")
        invisible(out)
      }
      
      
      # 2) Ejecutar por contraste (requiere res_sig_a/b/c y contr_1*/contr_2* ya definidos)
      
      # a) B vs A
      DEG_annotated_a <- annotate_and_export(res_sig_a, contr_2a, contr_1a, "a")
      
      # b) A vs C
      DEG_annotated_b <- annotate_and_export(res_sig_b, contr_2b, contr_1b, "b")
      
      # c) B vs C
      DEG_annotated_c <- annotate_and_export(res_sig_c, contr_2c, contr_1c, "c")
      
      
      
      
      
      # OTROS COMENTARIOS PARA DOCUMENTOS EN EXCEL -------------------------
      
      ### BLUE-RED KEGG COLOR----------------------
      #=SI(L2>0, "green", "purple")
