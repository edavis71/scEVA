# dependencies
suppressMessages(library(SingleCellExperiment))

# pathway class used for input
setClass('Pathway', slot=c(
    genes='character',
    control_variance='numeric',
    treatment_variance='numeric',
    gene_mean='numeric',
    gene_variance='numeric'
))

setMethod('initialize', 'Pathway',
    function(.Object, ...)
    {
        if (is.null(list(...)$control_variance))
            .Object@control_variance <- 0
        if (is.null(list(...)$treatment_variance))
            .Object@treatment_variance <- 0
        if (is.null(list(...)$gene_mean))
            .Object@gene_mean <- 1
        if (is.null(list(...)$gene_variance))
            .Object@gene_variance <- 3

        .Object <- callNextMethod(.Object, ...)
        return(.Object)
    }
)

Pathway <- function(genes, control_variance=0, treatment_variance=0,
gene_mean=1, gene_variance=3)
{
    return(new('Pathway', genes=genes, control_variance=control_variance,
        treatment_variance=treatment_variance, gene_mean=gene_mean,
        gene_variance=gene_variance))
}

# randomly generate multiplication factors from a log-normal distribution
getLNormFactors <- function(n.facs, sel.prob, neg.prob, fac.loc, fac.scale)
{
    is.selected <- as.logical(rbinom(n.facs, 1, sel.prob))
    n.selected <- sum(is.selected)
    dir.selected <- (-1) ^ rbinom(n.selected, 1, neg.prob)
    facs.selected <- rlnorm(n.selected, fac.loc, fac.scale)
    # Reverse directions for factors that are less than one
    dir.selected[facs.selected < 1] <- -1 * dir.selected[facs.selected < 1]
    factors <- rep(1, n.facs)
    factors[is.selected] <- facs.selected ^ dir.selected
    return(factors)
}

# Implementation of the logistic function
logistic <- function(x, x0, k)
{
    1 / (1 + exp(-k * (x - x0)))
}

# based on splatter
differentialVarianceSimulation <- function(pathways, nCellsPerGroup, lib.loc=11,
lib.scale=0.2, out.prob=0.05, out.facLoc=4, out.facScale=0.5, bcv.common=0.1,
bcv.df=60, dropout.mid=0, dropout.shape=-1)
{
    # check inputs
    if (!all(sapply(pathways, function(p) is(p, 'Pathway'))))
        stop("pathways must be a list of Pathway objects")

    # Set up name vectors
    control.names <- paste0("control.", seq_len(nCellsPerGroup))
    treatment.names <- paste0("treatment.", seq_len(nCellsPerGroup))
    cell.names <- c(control.names, treatment.names)
    gene.names <- do.call(c, lapply(pathways, function(p) p@genes))
    group.names <- c(rep("control", nCellsPerGroup),
        rep("treatment", nCellsPerGroup))
    
    # Create SingleCellExperiment to store simulation
    cells <- data.frame(Cell=cell.names)
    features <- data.frame(Gene=gene.names)
    rownames(cells) <- cell.names
    rownames(features) <- gene.names
    sim <- SingleCellExperiment(rowData=features, colData=cells)
    colData(sim)$Group <- group.names

    # simulate expected library sizes (no normalization)
    nCells <- length(cell.names)
    exp.lib.sizes <- rlnorm(nCells, lib.loc, lib.scale)
    colData(sim)$ExpLibSize <- exp.lib.sizes

    # simulate gene means (with outliers)
    nGenes <- length(gene.names)
    pathway.means.gene <- lapply(pathways, function(p)
        {
            mean <- p@gene_mean
            var <- p@gene_variance
            return(rgamma(length(p@genes), shape=mean^2/var, scale=var/mean))
        })
    base.means.gene <- do.call(c, pathway.means.gene)
    outlier.facs <- getLNormFactors(nGenes, out.prob, 0, out.facLoc, out.facScale)
    median.means.gene <- median(base.means.gene)
    outlier.means <- median.means.gene * outlier.facs
    is.outlier <- outlier.facs != 1
    means.gene <- base.means.gene
    means.gene[is.outlier] <- outlier.means[is.outlier]
    rowData(sim)$BaseGeneMean <- base.means.gene
    rowData(sim)$OutlierFactor <- outlier.facs

    # simulate cell means
    cell.means.cell <- matrix(1, ncol=nCells, nrow=nGenes) * means.gene
    rownames(cell.means.cell) <- gene.names
    for (pwy in pathways)
    {
        control_ndx <- 1:nCellsPerGroup
        treatment_ndx <- (nCellsPerGroup + 1):ncol(cell.means.cell)
        size <- length(pwy@genes) * nCellsPerGroup
        control_var <- matrix(rnorm(size, 1, pwy@control_variance), ncol=nCellsPerGroup)
        control_var <- pmax(control_var, 0)
        treatment_var <- pmax(control_var, 0)
        treatment_var <- matrix(rnorm(size, 1, pwy@treatment_variance), ncol=nCellsPerGroup)
        cell.means.cell[pwy@genes,control_ndx] <-
            cell.means.cell[pwy@genes,control_ndx] * control_var
        cell.means.cell[pwy@genes,treatment_ndx] <-
            cell.means.cell[pwy@genes,treatment_ndx] * treatment_var
    }
    cell.means.cell <- pmax(cell.means.cell, 0.1)
    cell.props.gene <- t(t(cell.means.cell) / colSums(cell.means.cell))
    base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)
    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- gene.names
    assays(sim)$BaseCellMeans <- base.means.cell

    #### simulate BCV means
    bcv <- (bcv.common + (1 / sqrt(base.means.cell))) * 
        sqrt(bcv.df / rchisq(nGenes, df = bcv.df))
    means.cell <- matrix(rgamma(nGenes * nCells, shape = 1 / (bcv ^ 2),
        scale = base.means.cell * (bcv ^ 2)), nrow = nGenes, ncol = nCells)
    colnames(means.cell) <- cell.names
    rownames(means.cell) <- gene.names
    assays(sim)$BCV <- bcv
    assays(sim)$CellMeans <- means.cell

    #### simulate true counts
    true.counts <- matrix(rpois(nGenes * nCells, lambda = means.cell),
        nrow = nGenes, ncol = nCells)
    colnames(true.counts) <- cell.names
    rownames(true.counts) <- gene.names
    assays(sim)$TrueCounts <- true.counts

    #### simulate dropout
    dropout.mid <- rep(dropout.mid, nCells)
    dropout.shape <- rep(dropout.shape, nCells)
    drop.prob <- sapply(seq_len(nCells), function(idx)
        {
            eta <- log(means.cell[,idx])
            return(logistic(eta, x0 = dropout.mid[idx], k = dropout.shape[idx]))
        })
    keep <- matrix(rbinom(nCells * nGenes, 1, 1 - drop.prob), nrow = nGenes,
        ncol = nCells)
    counts <- true.counts * keep
    colnames(drop.prob) <- cell.names
    rownames(drop.prob) <- gene.names
    colnames(keep) <- cell.names
    rownames(keep) <- gene.names
    assays(sim)$DropProb <- drop.prob
    assays(sim)$Dropout <- !keep
    BiocGenerics::counts(sim) <- counts

    return(sim)
}
