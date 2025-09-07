# Challenge
*Coffea arabica* challenge for RUST stress

**_Coffea arabica_ challenge adressed by RNAseq under different types of mutagenesis at different levels of RUST infection**

## CODE

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
      counts <- counts[rowSums(counts) > 0, , drop = FALSE]
      cat("After zero-row filter:", nrow(counts), "genes and", ncol(counts), "samples\n")
            # cat() → prints a quick summary

### Grouping into T0 & H24

      ## Build colData for grouping (T0 vs H24) 
      samples <- colnames(counts)
      time <- ifelse(grepl("^T", samples), "T0", "H24")
            # grepl("^T", samples) → TRUE if the sample name starts with "T"
            # it TRUE the it assigns T0, else H24
      time <- factor(time, levels = c("T0", "H24"))
            # factor() → categorical variable
      coldata <- data.frame(time = time, row.names = samples)
            # data.frame() → Build up the metadata "table"

### 

      ## ========= STEP 2: FPM screen (>1 in ≥ half of reps per group) =========
      # DESeq2 needs integer counts; round if your file has decimals
      dds_raw <- DESeqDataSetFromMatrix(countData = as.matrix(round(counts)),
                                        colData = coldata,
                                        design = ~ time)
      
      # Raw FPM/CPM-like (before size-factor normalization) for *presence* filtering
      cpm_raw <- fpm(dds_raw, robust = FALSE)   # keep FALSE to match our pipeline
      
      # Keep a gene if FPM>1 in at least half of the replicates in ANY group (T0 or H24)
      grp <- coldata$time
      mat_eval <- cpm_raw > 1
      keep <- rep(FALSE, nrow(mat_eval))
      cutoff <- ceiling(table(grp) / 2)  # half of replicates per group
      
      for (g in levels(grp)) {
        idx <- which(grp == g)
        keep <- keep | (rowSums(mat_eval[, idx, drop = FALSE]) >= cutoff[g])
      }
      
      cat("Genes kept by FPM>1 prevalence rule:", sum(keep), "of", length(keep),
          sprintf("(%.1f%%)\n", 100 * sum(keep) / length(keep)))
      
      # Subset your objects for downstream steps
      counts_filt <- counts[keep, , drop = FALSE]
      dds_1 <- dds_raw[keep, ]  # we’ll use this for size-factor normalization next
      
        # Quick sanity check
        stopifnot(identical(colnames(counts_filt), rownames(coldata)))


















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
