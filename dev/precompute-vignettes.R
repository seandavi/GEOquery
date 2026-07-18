#!/usr/bin/env Rscript
# Precompile GEOquery vignettes.
#
# The vignettes hit the live NCBI GEO network, which must NOT happen during
# `R CMD build` / `R CMD check` on CRAN/Bioconductor. So we author each vignette
# in a `vignettes/<name>.qmd.orig` source (live `eval: true` code) and precompile
# it here into a static `vignettes/<name>.qmd` with the real output baked in.
# The shipped `.qmd` executes nothing at build time; the `quarto::html` engine
# just renders markdown -> HTML. This is the knitr `*.Rmd.orig` pattern adapted
# for Quarto vignettes.
#
# Usage (from the package root, with network + Suggests deps installed):
#
#     Rscript dev/precompute-vignettes.R            # all vignettes
#     Rscript dev/precompute-vignettes.R rnaseq     # just vignettes/rnaseq.*
#
# Then commit BOTH the edited `.qmd.orig` and the regenerated `.qmd`.
# Refresh whenever you edit a vignette's code or want to pick up GEO-side changes.

suppressWarnings(suppressMessages({
    library(knitr)
    pkgload::load_all(".", quiet = TRUE)   # use the working tree, not the install
}))

vig_dir <- "vignettes"
args <- commandArgs(trailingOnly = TRUE)

orig_files <- if (length(args)) {
    file.path(vig_dir, paste0(args, ".qmd.orig"))
} else {
    list.files(vig_dir, pattern = "\\.qmd\\.orig$", full.names = TRUE)
}

missing <- orig_files[!file.exists(orig_files)]
if (length(missing)) {
    stop("no such source(s): ", paste(missing, collapse = ", "))
}

# `#>` output prefix + collapsed code/output to match the house style.
knitr::opts_chunk$set(
    comment = "#>",
    collapse = TRUE,
    message = FALSE,        # keep download chatter / one-time notices out of docs
    fig.path = "figures/"   # only used if a doc emits a plot; commit any output
)

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
for (orig in sort(orig_files)) {
    out <- sub("\\.qmd\\.orig$", ".qmd", orig)
    message("precompiling ", orig, " -> ", out)
    # knit() resolves relative paths (figures) against the working dir; run it
    # from inside vignettes/ so any assets land alongside the vignette.
    setwd(vig_dir)
    knitr::knit(basename(orig), basename(out), quiet = TRUE)
    setwd(old_wd)
}

# The single-cell manifest step can leave empty per-accession download dirs
# (GSE*/GSM*) under vignettes/; remove them so they never reach the tarball.
stray <- list.dirs(vig_dir, recursive = FALSE)
stray <- stray[grepl("/GS[EM][0-9]+$", stray)]
for (d in stray) {
    if (length(list.files(d, all.files = TRUE, no.. = TRUE)) == 0L) unlink(d, recursive = TRUE)
}

message("done. Commit the .qmd.orig sources and the regenerated .qmd files.")
