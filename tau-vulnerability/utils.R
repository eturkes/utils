#    This file is part of tau-vulnerability.
#    Copyright (C) 2019-2020  Emir Turkes, UK DRI at UCL, Columbia University Medical Center
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Emir Turkes can be contacted at emir.turkes@eturkes.com

# This file holds common functions and methods.

#' ggplot2 function providing custom aesthetics and automatic placement of categorical labels.
#' For continuous data, a colorbar is implemented.
#'
#' @param data SingleCellExperiment or Seurat object.
#' @param x,y Dimensionality reduction coordinates.
#' @param color Column metadata to color points by.
#' @param type \code{"cat"} is categorical, \code{"cont"} is continuous, \code{"NULL"} is generic.
#' @examples
#' red_dim_plot(data = sce, x = "tsne1", y = "tsne2", color = "cluster", type = "cat")
#' red_dim_plot(data = seurat, x = "umap1", y = "umap2", color = "nUMI", type = "cont")
#'
red_dim_plot <- function(data, x, y, color, type = NULL) {

  if ((class(data))[1] == "SingleCellExperiment") {
    gg_df <- data.frame(colData(data)[ , c(x, y, color)])
  } else if ((class(data))[1] == "Seurat") {
    gg_df <- data.frame(data[[x]], data[[y]], data[[color]])
  }
  rownames(gg_df) <- NULL
  gg_df[[color]] <- factor(gg_df[[color]])

  gg <- ggplot(gg_df, aes_string(x, y, col = color)) +
    geom_point(alpha = 0.35, stroke = 0.05, shape = 21, aes_string(fill = color)) +
    theme_classic() +
    theme(
      legend.position = "right", plot.title = element_text(hjust = 0.5),
      legend.title = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1)))

  if (is.null(type)) {
    return(gg)

  } else if (type == "cat") {
    label_df <- gg_df %>% group_by_at(color) %>% summarise_at(vars(x:y), median)
    label_df <- cbind(label_df[[1]], label_df)
    names(label_df) <- c("label", color, x, y)
    gg <- gg + geom_label_repel(data = label_df, aes(label = label), show.legend = FALSE)

  } else if (type == "cont") {
    if ((class(data))[1] == "SingleCellExperiment") {
      gg_df <- data.frame(colData(data)[ , c(x, y, color)])
    } else if ((class(data))[1] == "Seurat") {
      gg_df <- data.frame(data[[x]], data[[y]], data[[color]])
    }
    rownames(gg_df) <- NULL

    gg <- ggplot(gg_df, aes_string(x, y)) +
      geom_point(alpha = 0.35, stroke = 0.05, aes_string(color = color)) +
      theme_classic() +
      theme(
        legend.position = "right", plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()
      ) +
      scale_color_viridis()
  }
  gg
}

#' Adds download buttons and horizontal scrolling to \code{"DT::datatable"}.
#'
#' @param dt A data.table object.
#' @examples
#' datatable_download(dt = data_table)
#'
datatable_download <- function(dt) {

  datatable(
    dt,
    list(
      scrollX = TRUE, dom = "Blfrtip",
      buttons = list(
        "copy", "print",
        list(extend = "collection", buttons = c("csv", "excel", "pdf"), text = "Download")
      )
    ),
    extensions = "Buttons"
  )
}

#' Adds download buttons, horizontal scrolling, exponential values to \code{"DT::datatable"}.
#'
#' @param dt A data.table object.
#' @examples
#' datatable_download_exp(dt = data_table)
#'
datatable_download_exp <- function(dt) {

  datatable(
    dt,
    list(
      scrollX = TRUE,
      dom = "Blfrtip",
      buttons = list(
        "copy", "print",
        list(extend = "collection", buttons = c("csv", "excel", "pdf"), text = "Download")
      ),
      rowCallback = JS(
        "function(row, data) {",
        "for (i = 1; i < data.length; i++) {",
        "if (data[i]>=1000 | data[i]<1000) {",
        "$('td:eq('+i+')', row).html(data[i].toExponential(2));}}}"
      )
    ),
    extensions = "Buttons"
  )
}

#' Convert human to mouse gene names.
#' Adapted from:
#' https://rjbioinformatics.com/2016/10/14/converting-mouse-to-human-gene-names-with-biomart-package/
#'
#' @param genes A vector of human genes.
#' @examples
#' human_to_mouse_genes(genes = gene_list)
#'
human_to_mouse_genes <- function(genes) {

  human <- useMart("ensembl", "hsapiens_gene_ensembl")
  mouse <- useMart("ensembl", "mmusculus_gene_ensembl")

  new_genes <- getLDS(
    attributes = "external_gene_name", filters = "external_gene_name", values = genes,
    mart = human, attributesL = "external_gene_name", martL = mouse
  )
  new_genes <- unique(new_genes[ , 2])
}

#' Pipeline for normalization, dimensionality reduction, and clustering of post-QC scRNA-seq data.
#'
#' @param seurat Post-QC Seurat object.
#' @param cache_dir Directory to save post-processed Seurat object.
#' @param sub_name Subset level for naming of cache object.
#' @param protocol Vector with the following elements in this order: \code{"human"} or
#' \code{"mouse"}. \code{"droplet"} or \code{"smart-seq"}. \code{"single-cell"} or
#' \code{"single-nuc"}. \code{"umis"} or \code{"reads"}.
#' @param vars_to_regress Vector of nuisance variables for sctransform to regress out.
#' @param parallel_override See function \code{"parallel_plan"}.
#' @param cc Logical, whether to perform cell-cycle scoring.
#' @examples
#' cluster_pipeline(
#'   seurat = seurat, cache_dir = cache_dir, sub_name = "neuronal", protocol = protocol,
#'   vars_to_regress = "mito_percent", parallel_override = NULL, cc = FALSE)
#' )
cluster_pipeline <- function(
  seurat, cache_dir, sub_name, protocol, vars_to_regress, parallel_override, cc = TRUE
) {

  rds <- file.path(cache_dir, paste0(sub_name, "_seurat.rds"))
  if (file.exists(rds)) {
    seurat <- readRDS(rds)
    return(seurat)
  } else {

    if (protocol[4] == "umis") {
      # Run sctransform.
      # Note that this function produces many iterations of the following benign warning:
      # Warning in theta.ml(y = y, mu = fit$fitted): iteration limit reached
      # ---------------------------------------------------------------------------------
      parallel_plan(seurat, parallel_override)
      seurat <- suppressWarnings(
        SCTransform(seurat, vars.to.regress = vars_to_regress, verbose = FALSE)
      )
      # ---------------------------------------------------------------------------------

      # Perform PCA.
      # ------------
      seurat <- RunPCA(seurat, verbose = FALSE)
      add_df <- data.frame(Embeddings(seurat)[ , 1:2])
      names(add_df) <- paste0("pca", seq(ncol(add_df)))
      seurat$pca1 <- add_df$pca1
      seurat$pca2 <- add_df$pca2
      reduction <- "pca"
      dims <- 1:30 # Dimensions for downstream computations.
      # ------------

    } else if (protocol[4] == "reads") {
      sce <- as.SingleCellExperiment(seurat) # ZINB-WaVE can directly take an SCE but not Seurat.
      logcounts(sce) <- NULL

      # Use top 1,000 variable features for downstream computations.
      # ------------------------------------------------------------
      seurat <- NormalizeData(seurat, verbose = FALSE)
      seurat <- FindVariableFeatures(seurat, nfeatures = 1000, verbose = FALSE)
      sce <- sce[rownames(sce) %in% VariableFeatures(seurat), ]
      # ------------------------------------------------------------

      # Get ENSEMBL annotations for more accurate retrieval of gene length and GC content.
      # ----------------------------------------------------------------------------------
      if (protocol[1] == "human") {
        dataset <- "hsapiens_gene_ensembl"
      } else if (protocol[1] == "mouse") {
        dataset <- "mmusculus_gene_ensembl"
      }
      mart <- useEnsembl("ensembl", dataset)
      attributes <- c("external_gene_name", "ensembl_gene_id")
      gene_anno <- getBM(attributes, "external_gene_name", rownames(sce), mart)
      # ----------------------------------------------------------------------------------

      # For gene symbols with multiple ENSEMBL IDs, duplicate the gene symbol to have an identical
      # row for each ENSEMBL ID.
      # ------------------------------------------------------------------------------------------
      dup <- gene_anno[duplicated(gene_anno$external_gene_name), ]
      for (i in 1:dim(dup)[1]) {
        for (j in 1:dim(gene_anno)[1]) {
          if (dup$ensembl_gene_id[i] == gene_anno$ensembl_gene_id[j]) {
            gene_anno$external_gene_name[j] <- paste0(gene_anno$external_gene_name[j], "-dup")
          }
        }
      }
      sce <- sce[rownames(sce) %in% gene_anno$external_gene_name, ]
      new_mat <- counts(sce)
      for (i in 1:dim(dup)[1]) {
        for (j in 1:dim(sce)[1]) {
          if (dup$external_gene_name[i] == rownames(sce)[j]) {
            new_row <- t(counts(sce)[j, ])
            rownames(new_row) <- paste0(rownames(counts(sce))[j], "-dup")
            new_mat <- rbind(new_mat, new_row)
          }
        }
      }
      gene_anno <- gene_anno[gene_anno$external_gene_name %in% rownames(new_mat), ]
      gene_anno <- gene_anno[order(match(gene_anno$external_gene_name, rownames(new_mat))), ]
      rownames(new_mat) <- gene_anno$ensembl_gene_id
      sce <- SingleCellExperiment(list(counts = new_mat), colData = colData(sce))
      # ------------------------------------------------------------------------------------------

      # Get gene lengths and GC content.
      # --------------------------------
      row_data <- data.frame(getGeneLengthAndGCContent(rownames(sce), dataset, "biomart"))
      rowData(sce) <- data.frame(gc_content = row_data$gc, length = row_data$length)
      # --------------------------------

      # Run ZINB-WaVE.
      # TODO: Fix hardcoding of `vars_to_regress`.
      # ------------------------------------------
      counts(sce) <- as.matrix(counts(sce)) # ZINB-WaVE is incompatible with dgCMatrix.
      sce <- zinbwave(
        sce, paste0("~ ", vars_to_regress[1], " + ", vars_to_regress[2]), "~ gc_content + length",
        10, epsilon = 1e12, BPPARAM = MulticoreParam()
      )

      zinb <- data.frame(reducedDim(sce, "zinbwave"))
      seurat$zinb1 <- zinb$W1
      seurat$zinb2 <- zinb$W2
      zinb <- as.matrix(zinb)
      seurat[["zinb"]] <- CreateDimReducObject(zinb, key = "W", assay = DefaultAssay(seurat))
      reduction <- "zinb"
      dims <- 1:10 # Dimensions for downstream computations.
      # ------------------------------------------
    }

    # Perform cell cycle scoring.
    # ---------------------------
    if (cc == TRUE) {
      if (protocol[1] == "human") {
        s_genes <- cc.genes.updated.2019$s.genes
        g2m_genes <- cc.genes.updated.2019$g2m.genes
      } else if (protocol[1] == "mouse") {
        s_genes <- human_to_mouse_genes(cc.genes.updated.2019$s.genes)
        g2m_genes <- human_to_mouse_genes(cc.genes.updated.2019$g2m.genes)
      }
      seurat <- CellCycleScoring(seurat, s_genes, g2m_genes)
      seurat$cc_diff <- seurat$S.Score - seurat$G2M.Score # Combined proliferating cell signal.
    }
    # ---------------------------

    # Perform UMAP reduction.
    # -----------------------
    seurat <- RunUMAP(seurat, dims, reduction, min.dist = 0.75, verbose = FALSE)
    add_df <- data.frame(Embeddings(seurat, "umap"))
    names(add_df) <- paste0("umap", seq(ncol(add_df)))
    seurat$umap1 <- add_df$umap1
    seurat$umap2 <- add_df$umap2
    # -----------------------

    # Perform Louvain clustering.
    # ---------------------------
    resolution <- (dim(seurat)[2] / 3000) * 0.8 # Default is optimal for 3K cells so we scale it.
    seurat <- FindNeighbors(seurat, reduction, dims, verbose = FALSE)
    seurat <- FindClusters(seurat, resolution = resolution, verbose = FALSE)
    # ---------------------------

    saveRDS(seurat, rds)
  }
  seurat
}

#' Set the \code{"plan"} for \code{"future"} based on free memory and object size with the option to
#' override.
#'
#' @param object Object to check if \code{"future.globals.maxSize"} large enough to parallelize.
#' @param parallel_override \code{"NULL"} to calculate plan decision, \code{0} for sequential, a
#' non-zero integer for multiprocess and to set \code{"future.globals.maxSize"}.
#' @examples
#' parallel_plan(object = seurat, parallel_override = 5368709120)
#'
parallel_plan <- function(object, parallel_override = NULL) {

  if (is.null(parallel_override)) {
    # Get free memory.
    # ----------------
    gc()
    mem <- as.numeric(unlist(strsplit(system("free -b", TRUE)[2], " "))[7])
    # ----------------

    # Distribute free memory (minus 10 GiB) across available cores.
    # -------------------------------------------------------------
    mem <- mem - 10 * 1024 ^ 3
    mem <- mem / detectCores()
    # -------------------------------------------------------------

    # Enable parallelization only if `object` can fit in `future.globals.maxSize` (plus 1 Gib).
    # -----------------------------------------------------------------------------------------
    if (mem > object.size(object) + 1 * 1024 ^ 3) {
      plan("multiprocess")
      options(future.globals.maxSize = mem)
    } else {
      plan("sequential")
    }
    # -----------------------------------------------------------------------------------------

  } else if (parallel_override == 0) {
    plan("sequential")

  } else {
    plan("multiprocess")
    options(future.globals.maxSize = parallel_override)
  }
}


#' Plot genes of a Seurat object in various ways.
#'
#' @param subset_list A list containing Seurat objects.
#' @param subset_names A character vector of names for each Seurat object.
#' @param genes A character vector of genes to plot.
#' @examples
#' gene_plot(subset_list = sublist, subset_names = c("seurat1, "seurat2"), genes = c("TBR1, "GAD1))
#'
gene_plot <- function(subset_list, subset_names, genes) {

  for (i in 1:length(subset_list)) {

    if (length(which(genes %in% rownames(subset_list[[i]]))) > 4) {
      ncol <- 3
    } else {
      ncol <- 2
    }

    print(subset_names[i])
    if (length(genes) == 2 && all(genes %in% rownames(subset_list[[i]]))) {
      plots <- FeaturePlot(
        subset_list[[i]], genes, order = TRUE, blend = TRUE,
        combine = FALSE, blend.threshold = 0, max.cutoff = "q10"
      )
      print(CombinePlots(plots[3:4], legend = "none") + NoLegend())
    } else if (any(genes %in% rownames(subset_list[[i]]))) {
      print(
        FeaturePlot(
          subset_list[[i]], genes, order = TRUE, cols = c("lightgrey", "red"), ncol = ncol
        )
      )
    }

    if (any(genes %in% rownames(subset_list[[i]]))) {
      print(DotPlot(subset_list[[i]], features = genes, cols = c("blue", "red")) + RotatedAxis())
      print(
        DoHeatmap(subset_list[[i]], genes, slot = "data", size = 3) +
          NoLegend() +
          scale_fill_gradient2(
            low = rev(c('#D1E5F0','#67A9CF','#2166AC')), mid = "white",
            high = rev(c('#B2182B','#EF8A62','#FDDBC7')), midpoint = 0,
            guide = "colourbar", aesthetics = "fill"
          )
      )
      print(VlnPlot(subset_list[[i]], genes, ncol = ncol - 1))
    }
  }
}

############################################################
#
# author: Ludwig Geistlinger
# date: 2015-03-10 13:32:37
#
# descr: get.gene.length.and.gc.content
# update: 2015-06-14 exonic sequences
############################################################

#
# @input:
#   - id: one or more gene IDs (ensembl or entrez)
#   - org: organism three letter code, e.g. 'hsa' for 'Homo sapiens'
#   - mode: 1. biomart (supports all ensembl organisms, but might be a time-consuming)
#           2. org.db (based on BioC annotation, which is much faster but only
#                    for organisms with a respective TxDb, BSgenome, and OrgDb package)
getGeneLengthAndGCContent <- function(id, org, mode=c("biomart", "org.db"))
{
  id.type <- .autoDetectGeneIdType(id[1])
  if(is.na(id.type))
    stop("Only ENTREZ or ENSEMBL gene IDs are supported.")

  mode <- match.arg(mode)
  inp.id <- id

  # (a) based on BioC annotation utilities:
  #       (0) OrgDb: map between identifiers (if necessary)
  #       (1) TxDB: get genomic coordinates of genes
  #       (2) BSgenome: get sequences of genomic coordinates
  #
  if(mode=="org.db")
  {
    # check for TxDb package
    txdb.pkg <- .org2pkg(org, type="TxDb")
    .isAvailable(txdb.pkg, type="TxDb")

    # check for BSgenome package
    bsgen.pkg <- .org2pkg(org, type="BSgenome")
    .isAvailable(bsgen.pkg, type="BSgenome")

    txdb.spl <- unlist(strsplit(txdb.pkg, "\\."))
    txdb.id.type <- txdb.spl[length(txdb.spl)]
    if(txdb.id.type == "ensGene") {
      txdb.id.type <- "ensembl"
    } else if(txdb.id.type == "knownGene") {
      txdb.id.type <- "entrez"
    } else if(txdb.id.type == "sgdGene") {
      txdb.id.type <- "sgd"
    } else {
      stop(paste("TxDb does not use ENSEMBL or ENTREZ gene IDs"))
    }

    # (0) map ensembl <-> entrez,
    # if given id.type is entrez, but Txdb uses ensembl (or vice versa)
    if(id.type != txdb.id.type)
    {
      orgdb.pkg <- .org2pkg(org)
      .isAvailable(orgdb.pkg)
      orgdb.pkg <- get(orgdb.pkg)
      id.map <- mapIds(orgdb.pkg, keys=id,
                       column=ifelse(id.type == "entrez", "ENSEMBL", "ENTREZID"),
                       keytype=ifelse(id.type == "entrez", "ENTREZID", "ENSEMBL"))
      id <- id.map[!is.na(id.map)]
    }

    # (1) get genomic coordinates
    txdb.pkg <- get(txdb.pkg)
    coords <- exonsBy(txdb.pkg, by="gene")
    id <- id[id %in% names(coords)]
    coords <- reduce(coords[id])
    len <- sum(width(coords))

    # (2) get sequences
    bsgen.pkg <- get(bsgen.pkg)
    seqs <- getSeq(bsgen.pkg, coords)
    af <- alphabetFrequency(unlist(seqs, use.names=FALSE), baseOnly=TRUE, as.prob=TRUE)
    gc.cont <- mean(relist(rowSums(af[,c("C", "G")]), seqs))
  }
  # (b) based on BioMart
  #
  #
  else
  {
    id.type <- paste0(id.type, ifelse(id.type=="entrez", "gene", "_gene_id"))

    # setting mart
    message("Connecting to BioMart ...")
    ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
    ds <- listDatasets(ensembl)[,"dataset"]
    ds <- grep(paste0("^", org), ds, value=TRUE)
    if(length(ds) == 0)
      stop(paste("Mart not found for:", org))
    else if(length(ds) > 1)
    {
      message("Found several marts")
      sapply(ds, function(d)
        message(paste(which(ds==d), d, sep=": ")))
      n <- readline(paste0("Choose mart (1-", length(ds),") : "))
      ds <- ds[as.integer(n)]
    }

    ensembl <- useDataset(ds, mart=ensembl)

    message( paste0( "Downloading sequence",
                     ifelse(length(id) > 1, "s", ""), " ..."))
    if(length(id) > 100) message("This may take a few minutes ...")

    # download sequence
    # (1) get exon coordinates
    attrs <- c(id.type, "ensembl_exon_id",
               "chromosome_name", "exon_chrom_start", "exon_chrom_end")
    coords <- getBM(filters=id.type, attributes=attrs, values=id, mart=ensembl)
    id <- unique(coords[,id.type])
    coords <- GRangesList(sapply(id,
                                 function(i)
                                 {
                                   i.coords <- coords[coords[,1]== i, 3:5]
                                   g <- GRanges(i.coords[,1], IRanges(i.coords[,2],i.coords[,3]))
                                   return(g)
                                 }), compress=FALSE)
    coords <- reduce(coords)
    len <- sum(width(coords))

    # (2) get genes and sequences
    sel <- c(id.type, "start_position", "end_position")
    gene.pos <- getBM(attributes = sel, filters=id.type, values=id,
                      mart=ensembl)
    gene.seqs <- getSequence(id=id,
                             type=id.type, seqType="gene_exon_intron", mart=ensembl)

    # (3) get exonic sequences and correspondig GC content
    gc.cont <- sapply(id,
                      function(i)
                      {
                        # exon coordinates, gene position & sequence for current id i
                        ecoords <- coords[[i]]
                        gpos <- gene.pos[gene.pos[,id.type] == i,
                                         c("start_position", "end_position")]
                        gseq <- DNAString(
                          gene.seqs[gene.seqs[,id.type] == i, "gene_exon_intron"])

                        # exon coordinates relative to gene position
                        start <- start(ranges(ecoords)) - gpos[1,1] + 1
                        end <- end(ranges(ecoords)) - gpos[1,1] + 1
                        eseq <- gseq[IRanges(start, end)]
                        gc.cont <- sum(alphabetFrequency(eseq, as.prob=TRUE)[c("C","G")])
                        return(gc.cont)
                      }
    )
  }

  res <- cbind(len, gc.cont)
  colnames(res) <- c("length", "gc")
  rownames(res) <- id

  # (4) order according to input ids
  if(mode == "org.db")
    if(id.type != txdb.id.type)
      rownames(res) <- names(id)

  not.found <- !(inp.id %in% rownames(res))
  na.col <- rep(NA, sum(not.found))
  rn <- c(rownames(res), inp.id[not.found])
  res <- rbind(res, cbind(na.col, na.col))
  rownames(res) <- rn
  res <- res[inp.id,]
  return(res)
}


.isAvailable <- function(pkg, type="annotation")
{
  if(!(pkg %in% .packages(all.available=TRUE)))
  {
    message(paste0("Corresponding ", type,  " package not found: ",
                   pkg, "\nMake sure that you have it installed."))
    choice <- readline("Install it now? (y/n): ")
    if(choice == "y")
    {
      if (!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg)
    }
    else stop(paste("Package", pkg, "is not available"))
  }
  require(pkg, character.only = TRUE)
}

.autoDetectGeneIdType <- function(id)
{
  type <- NA
  if(grepl("^[Ee][Nn][Ss][A-Za-z]{0,3}[Gg][0-9]+", id)) type <- "ensembl"
  else if(grepl("^[0-9]+$", id)) type <- "entrez"
  else if(grepl("^[Yy][A-Za-z]{2}[0-9]{3}[A-Za-z]", id)) type <- "sgd"
  else if(grepl("^[Aa][Tt][0-9][A-Za-z][0-9]{5}", id)) type <- "tair"
  return(type)
}

.getOrgIdType <- function(org)
{
  it <- "eg"
  if(org == "At") it <- "tair"
  else if(org == "Pf") it <- "plasmo"
  else if(org == "Sc") it <- "sgd"
  return(it)
}

.supportedOrganisms <- function() sub(".db0$", "", .availableOrgPkgs())

.availableOrgPkgs <- function(type=c("OrgDb", "TxDb", "BSgenome"), local=TRUE)
{
  if(local) pkgs <- .packages(all.available=TRUE)
  else pkgs <- available.packages(paste0("http://bioconductor.org/",
                                         "packages/release/data/annotation/src/contrib"))[, "Package"]

  type <- match.arg(type)
  org.string <- "^org.[A-z][a-z]+.[a-z]+.db$"
  if(type == "TxDb")
    org.string <- "^TxDb.[A-Z][a-z]+.UCSC.[a-z]{2}[A-Za-z]*[0-9]{1,3}.[a-z]{3,5}Gene$"
  else if(type == "BSgenome")
    org.string <- "^BSgenome.[A-Z][a-z]+.UCSC.[a-z]{2}[A-Za-z]*[0-9]{1,3}$"
  org.pkgs <- grep(org.string, pkgs, value=TRUE)
  names(org.pkgs) <- NULL
  return(org.pkgs)
}

.org2pkg <- function(org, type=c("OrgDb", "TxDb", "BSgenome"))
{
  type <- match.arg(type)

  SPECIES <- rbind(
    c("anopheles", "Anopheles gambiae", "Ag", "aga", "anoGam", "7165"),
    c("arabidopsis", "Arabidopsis thaliana", "At", "ath", NA, "3702"),
    c("bovine", "Bos taurus", "Bt", "bta", "bosTau", "9913"),
    c("canine", "Canis familiaris", "Cf", "cfa", "canFam", "9615"),
    c("chicken", "Gallus gallus", "Gg", "gga", "galGal", "9031"),
    c("chimp", "Pan troglodytes", "Pt", "ptr", "PanTro", "9598"),
    c("ecoliK12", "Escherichia coli K12", "EcK12", "eco", NA, "562,83333,511145"),
    c("ecoliSakai", "Escherichia coli Sakai", "EcSakai", "ecs", NA, "83334"),
    c("fly", "Drosophila melanogaster", "Dm", "dme", "dm", "7227"),
    c("human", "Homo sapiens", "Hs", "hsa", "hg", "9606"),
    c("malaria", "Plasmodium falciparum", "Pf", "pfa", NA, "5833"),
    c("mouse", "Mus musculus", "Mm", "mmu", "mm", "10090"),
    c("pig", "Sus scrofa", "Ss", "ssc", "susScr", "9823"),
    c("rat", "Rattus norvegicus", "Rn", "rno", "rn", "10116"),
    c("rhesus", "Macaca mulatta", "Mmu", "mcc", "rheMac", "9544"),
    c("worm", "Caenorhabditis elegans", "Ce", "cel", "ce", "6239"),
    c("xenopus", "Xenopus laevis", "Xl", "xla", "NA", "8355"),
    c("yeast", "Saccharomyces cerevisiae", "Sc", "sce", "sacCer", "4932,559292"),
    c("zebrafish", "Danio rerio", "Dr", "dre", "danRer", "7955")
  )
  colnames(SPECIES) <- c("common", "tax", "bioc", "kegg", "ucsc", "ncbi")


  # org specification via
  # (a) 3-letter code, e.g. 'hsa'
  # (b) genome assembly, e.g. 'hg38'
  is.genome <- sub("[0-9]+$", "", org) %in% SPECIES[,"ucsc"]
  if(is.genome)
  {
    ucsc.id <- org
    i <- grep(sub("[0-9]+$", "", org), SPECIES[,"ucsc"])
    bioc.id <- SPECIES[i, "bioc"]
  }
  else
  {
    ind <- apply(SPECIES, 1, function(r) org %in% r)
    if(any(ind)) i <- which(ind)[1]
    else stop(paste0("unrecognized organism ID \'", org, "\'"))
    bioc.id <- SPECIES[i, "bioc"]
    ucsc.id <- SPECIES[i, "ucsc"]
  }

  # TxDB, BSgenome, or OrgDB package?
  if(type %in% c("TxDb", "BSgenome"))
  {
    pkg.string <- paste0("^", type, ".", bioc.id, "[a-z]+.UCSC.", ucsc.id)
    pkg <- grep(pkg.string, .availableOrgPkgs(type), value=TRUE)
    if(length(pkg) == 0)
      pkg <- grep(pkg.string, .availableOrgPkgs(type, local=FALSE), value=TRUE)
    if(length(pkg) == 0)
      stop(paste("No corresponding", type, "package for", org))
    else if(length(pkg) > 1)
    {
      message("Found several genome assemblies")
      sapply(pkg, function(p)
        message(paste(which(pkg==p), p, sep=": ")))
      n <- readline(paste0("Choose assembly (1-", length(pkg),") : "))
      pkg <- pkg[as.integer(n)]

      #message("Found several genome assemblies")
      #message(paste("Using latest:", pkg))
      #ver <- sapply(pkg,
      #    function(p)
      #    {
      #        spl <- unlist(strsplit(p, "\\."))
      #        ind <- length(spl)
      #        if(type == "TxDb") ind <- ind - 1
      #        ass <- spl[ind]
      #        ver <- sub("^[a-zA-Z]+", "", ass)
      #        return(as.integer(ver))
      #    })
      #pkg <- pkg[which.max(ver)]
    }
  }
  else
  {
    id.type <- .getOrgIdType(bioc.id)
    pkg <- paste("org", bioc.id, id.type, "db", sep=".")
  }
  return(pkg)
}

word_cloud = function(x, width = NULL){
  t.rW = c("cell", "process", "negative", "positive",
           "activity", "protein", "involved",
           "component", "level", "event", "organismal",
           "cellular", "pathway", "mediated", "dependent",
           "group", "target", "biocarta", "kegg",
           "reactome", "system", "nervous", "cells",
           "time", "structure", "whose", "progression",
           "formation", "divided", "can", "specific",
           "outcome", "two", "form", "one",
           "size", "forms", "becomes", "become",
           "generated", "frequency", "rate", "extent",
           "will", "organism", "inner", "wall",
           "walls", "mammals", "organized", "anatomical",
           "concentration", "directed", "towards", "higher",
           "gradient", "results", "change", "state",
           "terms", "etc", "result", "stops",
           "prevents", "reduces", "activates", "increases",
           "activation", "within", "series", "molecular",
           "modulates", "introducing", "thin", "thick",
           "past", "located", "generation", "molecule",
           "reactions", "pathways", "resulting", "compounds",
           "population", "body", "regulation", "early",
           "commonpartner", "region", "regulated", "means",
           "agent", "structures", "involving", "family",
           "decreases", "together", "set", "daughter",
           "subsequent", "outer", "via", "molecules",
           "comprising", "diameter", "occurring", "removes",
           "another", "adjacent", "presence", "events",
           "specialized", "features", "acquires", "association",
           "type", "include", "consist", "linked",
           "containing", "generate", "member", "widely",
           "usually", "position", "often", "tissue",
           "responsible", "diverse", "range", "leading",
           "contributes", "upper", "part", "slightly",
           "expansion", "internal", "steady", "external",
           "middle", "outflow", "site", "orderly",
           "switching", "also", "known", "effects",
           "defined", "larger", "organisms", "along",
           "powered", "action", "effects", "across",
           "asymmetry", "host", "animal", "pertains",
           "creation", "formed", "tissues", "found",
           "certain", "localized", "distinct", "increase",
           "appearance", "due", "following", "levels",
           "consists", "filters", "end", "intial",
           "either", "double", "introduced", "initiation",
           "addition", "contain", "sorting", "wide",
           "ability", "similar", "act", "primarily",
           "passed", "separation", "received", "steps",
           "initial", "production", "somatic", "shaping",
           "organ", "portion", "systems", "components",
           "key", "mature", "transition", "selection",
           "promotor", "factor", "processing", "cleavage",
           "response", "stimulus", "secretion", "produced",
           "physiologic", "synthesized", "bind", "receptor",
           "trigger", "signaling", "signal", "release",
           "downstream", "communication", "using", "contains",
           "secretion", "transport", "transporter", "movement",
           "response", "signaling", "adhesion", "maintenance",
           "force", "multiplication", "common", "flow",
           "stream", "specialization", "charged", "multicellular",
           "fluid", "proteins", "stem", "converted",
           "transported", "arm", "superfamily", "air",
           "small", "maintained", "parts", "groups",
           "xray", "element", "beam", "strut",
           "rod", "preb", "adaptation", "charges",
           "required", "receive", "characterize", "function",
           "central", "sprouting", "long", "carries",
           "outgoing", "left", "longterm", "neuronal",
           "relatively", "unspecialized", "multiple", "high",
           "number", "tolerance", "induction", "compound",
           "neuron", "work", "perform", "functions",
           "electrical", "biological", "condition", "composed",
           "gives", "neural", "begins", "ends",
           "carry", "innermost", "densely", "packed",
           "mostly", "border", "send", "parallel",
           "brush", "provide", "supplies", "twothirds",
           "principal", "main", "supply", "connecting",
           "commonly", "observed", "visibly", "may",
           "exist", "loosely", "associated", "clusters",
           "combination", "attractive", "fully", "functional",
           "strongly", "several", "enclosed", "individual",
           "give", "rise", "attains", "length",
           "joining", "class", "consisting", "occurs",
           "closely", "related", "may", "many",
           "positions", "classes", "including", "occurs",
           "outside", "location", "connected", "various",
           "decrease", "causes", "makes", "consequence",
           "proceeds", "begins", "circumstances", "require",
           "actual", "numbers", "difficult", "numerous",
           "species", "acts", "upon", "lower",
           "carried", "regular", "longitudinal", "array",
           "assisting", "correct", "pressure", "typically",
           "presentation", "expresses", "passes", "towards",
           "existing", "arising", "present", "side",
           "several", "brown", "polarity", "andor",
           "marginal", "tract", "zone", "indicating",
           "cord", "units", "organic", "chain",
           "distributed", "greater", "without", "concomitant",
           "substance", "dense", "core", "late",
           "substances", "establishment", "initiated", "surface",
           "combining", "involves", "solutes", "role",
           "interaction", "coupled", "classical", "entry",
           "clustering", "domain", "insertion", "disassembly",
           "determination", "sliding", "ion", "potential",
           "acid", "guidance", "assembly", "complex",
           "import", "gland", "fusion", "sequestered",
           "cytosol", "propagation", "channels", "membrane",
           "chemical", "acids", "longchain", "metabolic",
           "unsaturated", "sugar", "carbon", "atoms",
           "base", "signals", "binding", "bonds",
           "derived", "sequestering", "separated", "adaptive",
           "commitment", "signals", "potentials", "extension",
           "ending", "replication", "new", "storage",
           "strands", "sequence", "extrusion", "arrangement",
           "eukaryotic", "organelle", "bounded", "anchored",
           "basally", "elongation", "channel", "smooth",
           "pore", "ions", "systemic", "plasma",
           "anion", "much", "organs", "endogenous",
           "connective", "tuft", "surrounded", "capsule",
           "vertebrate", "intracellular", "precursor", "extracellular",
           "secreted", "cellcell", "proximal", "acute",
           "hydrogen", "alpha", "produces", "transducers",
           "activators", "convey", "cascade", "toward",
           "nuclear", "translocation", "decline", "transduction",
           "accomplished", "retraction", "noncanonical", "bonding",
           "transforming", "right", "ventral", "phase",
           "repulsive", "cues", "canal", "sacs",
           "plant", "natural", "pesticide", "insects",
           "derivative", "linkage", "saturated", "lengths",
           "step", "cycle", "removed", "three",
           "remain", "respectively", "reaction", "ring",
           "transmission", "increased", "white", "achieved",
           "decreased", "attachment", "anterior", "secretory",
           "temporary", "negatively", "nucleus", "droplets",
           "residues", "linkages", "bonding", "incorporation",
           "context", "transfer", "multisubunit", "subunits",
           "processes", "identical", "sister", "environment",
           "density", "detection", "stimuli", "exposure",
           "soluble", "solvents", "powers", "defense",
           "immature", "break", "coat", "exogenous",
           "origin", "tertiary", "cycles", "novo",
           "always", "dorsal", "processes", "nuclear",
           "transmission", "loss", "corpus", "removing",
           "input", "acquire", "modulation", "generally",
           "hemispheres", "covalently", "bonded", "bonding",
           "complete", "apical", "protrusion", "subcellular",
           "compartment", "striated", "examples", "driven",
           "doublet", "tube", "planar", "closure",
           "coordinated", "plane", "partly", "axis",
           "activated", "removal", "includes", "identity",
           "sheath", "stage", "occur", "chains", "cushion",
           "glands", "focal", "stages", "evoked",
           "foam", "confining", "lateral", "cavity",
           "third", "interactions", "node", "mediates",
           "surroundings", "cyclic", "intermediate", "mass",
           "selfpropelled", "bound", "excluding", "point",
           "modification", "particle", "reactive", "active",
           "segregation", "substrate", "gamma", "execution",
           "facial", "branches", "fibers", "metal",
           "residue", "radial", "satellite", "restore",
           "interacting", "selectively", "noncovalently", "biosynthetic",
           "inheritance", "links", "light", "link",
           "junction", "enables", "innate", "biosynthesis",
           "tubule", "bulb", "cohesion", "replicated",
           "lowdensity", "directly", "solute", "lens",
           "human", "roof", "poles", "structural",
           "gene", "specification", "allows", "direct",
           "exchange", "smoothened", "secondary", "controlled",
           "intrinsic", "relaxation", "strength", "animals",
           "bodies", "information", "bond", "products",
           "embedded", "direction", "characteristics", "rhythm",
           "changes", "major", "requires", "residing",
           "physiological", "capable", "mediator")
  txt = unlist(strsplit(x, " "))
  txt = Corpus(VectorSource(txt))
  txt = tm_map(txt, PlainTextDocument)
  txt = tm_map(txt, removePunctuation)
  txt = tm_map(txt, removeNumbers)
  txt = tm_map(txt, content_transformer(tolower))
  txt = tm_map(txt, removeWords, c(t.rW, stopwords("english")))
  # corpus = txt
  # txt = tm_map(txt, stemDocument)
  # txt = tm_map(txt, stemCompletion, corpus)
  tdm = TermDocumentMatrix(txt)
  m = as.matrix(tdm)
  word_freqs = sort(rowSums(m), decreasing=TRUE)
  word_freqs = word_freqs[word_freqs>1]
  if (is.null(width)) {
    word_freqs = paste(names(word_freqs), collapse=" ")
  } else {
    word_freqs = paste(names(word_freqs)[1:width], collapse="\n")
  }
  gsub("[[:space:]]?NA[[:space:]]?", "", word_freqs)
}

computeGeneSetsOverlapMax <- function(gSets, uniqGenes, min.sz=1, max.sz=Inf) {
  ## gSetsMembershipMatrix should be a (genes x gene-sets) incidence matrix

  gSetsMembershipMatrix <- incidence(gSets)
  gSetsMembershipMatrix <- t(gSetsMembershipMatrix[, colnames(gSetsMembershipMatrix) %in% uniqGenes])

  lenGsets <- colSums(gSetsMembershipMatrix)

  szFilterMask <- lenGsets >= max(1, min.sz) & lenGsets <= max.sz
  if (!any(szFilterMask))
    stop("No gene set meets the minimum and maximum size filter\n")

  gSetsMembershipMatrix <- gSetsMembershipMatrix[, szFilterMask]
  lenGsets <- lenGsets[szFilterMask]

  totalGsets <- ncol(gSetsMembershipMatrix)

  M <- t(gSetsMembershipMatrix) %*% gSetsMembershipMatrix

  M1 <- matrix(lenGsets, nrow=totalGsets, ncol=totalGsets,
               dimnames=list(colnames(gSetsMembershipMatrix), colnames(gSetsMembershipMatrix)))
  M2 <- t(M1)
  M.max <- matrix(0, nrow=totalGsets, ncol=totalGsets)
  M.max[M1 > M2] <- M1[M1 > M2]
  M.max[M2 >= M1] <- M2[M2 >= M1]
  overlapMatrix <- M / M.max

  return (overlapMatrix)
}
