suppressPackageStartupMessages({
	library(argparse)
	library(Rmpfr)
})
set.seed(0)

# Parse Arguments
parser <- ArgumentParser()
parser$add_argument("--chisq_per_nampc_file",type="character")
parser$add_argument("--ks_file",type="character")
parser$add_argument("--outfile",type="character")
args <- parser$parse_args()

# Vectorized: Calculate p-values for all k values at once
for(k in ks){
    chi_stats <- rowSums(all_res[, 1:k, drop=FALSE])
    all_res[[paste0("k", k, "_P")]] <- pchisq(chi_stats, df=k, lower.tail=TRUE)
}

# Vectorized: Find minimum p-value across all k (vectorized row operation)
p_cols <- paste0("k", ks, "_P")
all_res[["P"]] <- apply(as.matrix(all_res[, p_cols]), 1, min)

# Find small p-values
small_vals <- which(all_res[["P"]] < 1e-15)

# Vectorized: Multiple testing correction
all_res[["P"]] <- 1 - (1 - all_res[["P"]])^length(ks)
idx_result <- apply(as.matrix(all_res[, p_cols]), 1, which.min)

# Vectorized: Which k gave minimum p-value
all_res[["k"]] <- ks[as.integer(idx_result)]

# Pre-compute sum of T-squared for each k
for(k in ks){
    all_res[[paste0("k", k, "_sumTsq")]] <- rowSums(all_res[, 1:k, drop=FALSE])
}

# Vectorized p_mht_corrected function
p_mht_corrected <- function(sumTsq_vec, ks, prec_bits = 10000){
    # Vectorized: Calculate all p-values at once
    uncorr_ps <- exp(mpfr(
        pchisq(sumTsq_vec, df=ks, lower.tail=TRUE, log.p=TRUE),
        precBits=prec_bits
    ))
    sel_k <- as.character(ks[which.min(uncorr_ps)])
    
    # Multiple testing correction
    corr_p <- formatDec(
        mpfr(1, prec_bits) - (mpfr(1, prec_bits) - mpfr(min(uncorr_ps), prec_bits))^length(ks),
        precBits=64
    )
    return(c(corr_p, sel_k))
}

# Apply high-precision correction only to small p-values
if(length(small_vals) > 0){
    sumTsq_cols <- paste0("k", ks, "_sumTsq")
    all_res[small_vals, c("P", "k")] <- t(apply(
        all_res[small_vals, sumTsq_cols, drop=FALSE],
        1,
        p_mht_corrected,
        ks
    ))
}

all_res[small_vals,c("P", "k")] = t(apply(as.matrix(all_res[small_vals,paste0("k",ks,"_sumTsq"), drop=FALSE], ncol=length(ks)), 1, p_mht_corrected, ks))

write.table(all_res[,c("P", "k")], args$outfile, quote=FALSE, row.names=FALSE, sep = "\t")
