# Challenge
*Coffea arabica* challenge for RUST stress

**_Coffea arabica_ challenge adressed by RNAseq under different types of mutagenesis at different levels of RUST infection**

## CODE

Using DESeq2 normalization

Unsupervised plots: (PCA, clustering/heatmaps) 

Don’t use the DESeq2 design. They depend only on the expression matrix we feed

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

### Heatmao Gene x Sample

- VST (Variance Stabilizing Transform) → removes dependence of variance on mean counts, making genes comparable across samples

- Row-z (row scaling) → each gene’s expression values are z-scored (mean = 0, sd = 1 across samples) so patterns are relative per gene, not absolute counts.

- Correlation distance (1 – correlation) → similarity between samples is measured by **how similar their expression profiles are** (ignores absolute magnitude, focuses on shape).

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

This allows us to prepare data so that genes are comparable across samples and patterns (clusters, SOM, heatmaps) highlight relative up/down regulation, not absolute counts

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


### Heatmaps 


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
