# Challenge
*Coffea arabica* challenge for RUST stress

**_Coffea arabica_ challenge adressed by RNAseq under different types of mutagenesis at different levels of RUST infection**

## CODE

Using DESeq2 normalization

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

vst → applies DESeq2’s variance-stabilizing transformation

"Raw counts" follow a negative binomial distribution with variance dependent on mean, what is thy that makes PCA and heatmaps messy

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

### PCA plots

Principal Component Analysis (PCA) is a dimensionality reduction technique

In RNA-seq:

→ Each sample (T1…H16) has thousands of gene expression values (a high-dimensional space)

→ PCA finds new axes (principal components) that capture the largest sources of variation in the data

- PC1 = axis that explains the largest % of total variance.

- PC2 = second largest, orthogonal to PC1.

→ Plotting samples on PC1–PC2 shows which samples are more similar or different

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
      
      ## 3) PLOT
      ggplot(df_pca, aes(x, y, color = time, label = name)) +
        geom_point(size = 2.6) +
        ggrepel::geom_text_repel(show.legend = FALSE, max.overlaps = Inf, size = 3) +     # Size of each dot
        scale_color_manual(
          values = pal_time,                # <- uses your exact palette
          limits = c("T0","H24"),           # <- fixes legend order & ensures mapping
          name   = "Time"                   # <- legend title to match base plot
        ) +
        labs(
          title = "PCA – Samples (VST z-scored)",
          x = sprintf("PC%d (%.1f%%)", var_1, var_expl[var_1]),
          y = sprintf("PC%d (%.1f%%)", var_2, var_expl[var_2])
        ) +
        theme_classic(base_size = 12) +     # Relación gráfico vs títulos 
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5),      # Moves tittle
          legend.position = "right"
        )
      
      ## 4) save PCA matrices
      # write.csv(pca$x, file = "PCA_samples_scores.csv")      # samples in PC space
      # write.csv(pca$rotation, file = "PCA_gene_loadings.csv") # genes loadings

<img width="600" height="350" alt="image" src="https://github.com/user-attachments/assets/3f3b560a-073e-4185-897b-3920e5ed2f79" />







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
