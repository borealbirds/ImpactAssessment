#!/usr/bin/env Rscript
# analyze_seff.R
# Summarise sacct resource usage (collect_seff.sh -> seff_summary.psv) by
# coalition size, and emit recommended sbatch --mem / --time directives.
#   Rscript analyze_seff.R [path/to/seff_summary.psv]

args <- commandArgs(trailingOnly = TRUE)
psv  <- if (length(args) >= 1) args[1] else "logs/seff_summary.psv"

d <- read.delim(psv, sep = "|", stringsAsFactors = FALSE, colClasses = "character")

# MaxRSS is reported only on the .batch substep; everything else on the main row.
main <- d[grepl("^[0-9]+_[0-9]+$", d$JobID), ]
bat  <- d[grepl("\\.batch$", d$JobID), ]
bat$base <- sub("\\.batch$", "", bat$JobID)

parse_mem <- function(x) {                       # sacct mem string -> MB
  x <- trimws(x); out <- rep(NA_real_, length(x))
  ok <- nzchar(x) & x != "0"
  num  <- suppressWarnings(as.numeric(sub("[A-Za-z]+$", "", x[ok])))
  unit <- toupper(sub("^[0-9.]+", "", x[ok]))
  mult <- ifelse(unit == "K", 1/1024, ifelse(unit == "M", 1,
          ifelse(unit == "G", 1024, ifelse(unit == "T", 1024^2, 1))))
  out[ok] <- num * mult; out
}

main$maxrss_mb <- parse_mem(bat$MaxRSS[match(main$JobID, bat$base)])
main$reqmem_mb <- parse_mem(main$ReqMem)
main$elapsed_s <- suppressWarnings(as.numeric(main$ElapsedRaw))

cid  <- suppressWarnings(as.integer(sub("^coal_", "", main$JobName)))
main <- main[!is.na(cid), ]; cid <- cid[!is.na(cid)]
popcount <- function(n) vapply(n, function(x) sum(as.integer(intToBits(x))), integer(1))
main$cid  <- cid
main$size <- popcount(cid - 1L)

cat(sprintf("parsed %d coalition job-tasks across cids %d..%d\n",
            nrow(main), min(cid), max(cid)))
st <- table(main$State)
cat("states: ", paste(sprintf("%s=%d", names(st), st), collapse = "  "), "\n", sep = "")

# --- per-size summary over COMPLETED tasks --------------------------------
hdr <- sprintf("%-4s %-6s %-8s %9s %9s %9s %9s %8s %8s %8s",
               "size","nDone","OOM/TO","RSSminG","RSSmedG","RSSp95G","RSSpkG",
               "Tmed_h","Tp95_h","Tpk_h")
cat("\n", hdr, "\n", strrep("-", nchar(hdr)), "\n", sep = "")

rec <- list()
for (s in sort(unique(main$size))) {
  rows <- main[main$size == s, ]
  done <- rows[rows$State == "COMPLETED" & is.finite(rows$maxrss_mb), ]
  n_oom <- sum(rows$State == "OUT_OF_MEMORY")
  n_to  <- sum(rows$State == "TIMEOUT")
  if (nrow(done) == 0) {
    cat(sprintf("%-4d %-6d %-8s %9s %9s %9s %9s %8s %8s %8s\n",
                s, 0L, paste0(n_oom,"/",n_to), "-","-","-","-","-","-","-"))
    next
  }
  rss <- done$maxrss_mb / 1024; el <- done$elapsed_s / 3600
  qf <- function(v,p) as.numeric(quantile(v, p, na.rm = TRUE))
  cat(sprintf("%-4d %-6d %-8s %9.1f %9.1f %9.1f %9.1f %8.2f %8.2f %8.2f\n",
              s, nrow(done), paste0(n_oom,"/",n_to),
              min(rss), median(rss), qf(rss,.95), max(rss),
              median(el), qf(el,.95), max(el)))
  rec[[as.character(s)]] <- list(size=s, rss_pk=max(rss), el_pk=max(el),
                                 n_oom=n_oom, n_to=n_to,
                                 reqmem_g=median(done$reqmem_mb,na.rm=TRUE)/1024)
}

# --- recommended sbatch directives ----------------------------------------
roundup <- function(x, step) ceiling(x / step) * step
fmt_t <- function(h) { hh <- ceiling(h)
  if (hh < 24) sprintf("%02d:00:00", hh)
  else sprintf("%d-%02d:00:00", hh %/% 24, hh %% 24) }
cat("\n=== recommended directives (MEM x1.20, TIME x1.30; OOM/TO -> bump) ===\n")
cat(sprintf("%-4s %-7s %-11s %s\n","size","--mem","--time","basis"))
for (r in rec) {
  mem_g <- roundup(r$rss_pk * 1.20, 16)
  if (r$n_oom > 0) mem_g <- max(mem_g, roundup(r$reqmem_g * 1.25, 16))
  hrs <- r$el_pk * 1.30; if (r$n_to > 0) hrs <- hrs * 1.5
  basis <- sprintf("peakRSS %.0fG, peakT %.1fh, n=%s",
                   r$rss_pk, r$el_pk, "")
  basis <- sub(", n=$", "", basis)
  if (r$n_oom > 0) basis <- paste0(basis, sprintf(" [+%d OOM]", r$n_oom))
  if (r$n_to  > 0) basis <- paste0(basis, sprintf(" [+%d TIMEOUT]", r$n_to))
  cat(sprintf("%-4d %-7s %-11s %s\n", r$size, paste0(mem_g,"G"), fmt_t(hrs), basis))
}

# --- per-coalition CSV (within-size variance is large) --------------------
out <- main[main$State == "COMPLETED" & is.finite(main$maxrss_mb),
            c("cid","size","JobID","maxrss_mb","elapsed_s","reqmem_mb","State")]
out$maxrss_g  <- round(out$maxrss_mb/1024, 1)
out$elapsed_h <- round(out$elapsed_s/3600, 2)
agg <- aggregate(cbind(maxrss_g, elapsed_h) ~ cid + size, data = out, FUN = max)
agg <- agg[order(agg$size, agg$cid), ]
csv <- file.path(dirname(psv), "seff_by_coalition.csv")
write.csv(agg, csv, row.names = FALSE)
cat(sprintf("\nwrote per-coalition peaks: %s (%d coalitions)\n", csv, nrow(agg)))
