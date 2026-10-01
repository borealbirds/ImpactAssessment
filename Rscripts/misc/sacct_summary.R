# Summarise a `sacct -P --units=G` dump: per array task, peak memory (MaxRSS of the .batch
# step), wall time and CPU efficiency. Dump on Fir first (one line):
#   sacct -j <job> -P --units=G -o JobID,State,ExitCode,Submit,Start,End,Elapsed,TotalCPU,AllocCPUS,ReqMem,MaxRSS,NodeList > sacct_<job>.txt
# usage: Rscript sacct_summary.R sacct_<job>.txt <job>
args <- commandArgs(TRUE)
d <- read.table(args[1], sep = "|", header = TRUE, quote = "", fill = TRUE,
                stringsAsFactors = FALSE, comment.char = "")
job <- args[2]

secs <- function(x) {                      # [D-]HH:MM:SS, MM:SS.mmm
  vapply(x, function(s) {
    if (is.na(s) || s == "") return(NA_real_)
    days <- 0
    if (grepl("-", s)) { days <- as.numeric(sub("-.*", "", s)); s <- sub(".*-", "", s) }
    p <- as.numeric(strsplit(s, ":")[[1]])
    p <- c(rep(0, 3 - length(p)), p)
    days * 86400 + p[1] * 3600 + p[2] * 60 + p[3]
  }, numeric(1), USE.NAMES = FALSE)
}
gb <- function(x) {                        # "45.23G", "512M", "0"
  x[x == ""] <- NA
  u <- sub("^[0-9.]+", "", x); v <- as.numeric(sub("[A-Za-z]+$", "", x))
  v * c(K = 1 / 1024^2, M = 1 / 1024, G = 1, T = 1024)[ifelse(u == "", "G", u)]
}

top <- d[grepl(paste0("^", job, "_[0-9]+$"), d$JobID), ]
bat <- d[grepl(paste0("^", job, "_[0-9]+\\.batch$"), d$JobID), ]
top$task <- as.integer(sub(".*_", "", top$JobID))
bat$task <- as.integer(sub(".*_([0-9]+)\\.batch$", "\\1", bat$JobID))
t <- merge(top[, c("task", "State", "Elapsed", "TotalCPU", "AllocCPUS", "ReqMem")],
           data.frame(task = bat$task, maxrss_gb = gb(bat$MaxRSS)), by = "task", all.x = TRUE)
t$elapsed_h <- secs(t$Elapsed) / 3600
t$cpu_eff   <- secs(t$TotalCPU) / (secs(t$Elapsed) * as.numeric(t$AllocCPUS))
t <- t[order(t$task), ]

cat(sprintf("job %s: %d tasks; states: %s\n", job, nrow(t),
            paste(names(table(t$State)), table(t$State), sep = "=", collapse = ", ")))
cat("requested:", unique(t$ReqMem), "\n\n")
q <- c(0, .5, .9, .95, .99, 1)
cat("peak memory (GB), quantiles:\n"); print(round(quantile(t$maxrss_gb, q, na.rm = TRUE), 1))
cat("\nwall time (h), quantiles:\n");  print(round(quantile(t$elapsed_h, q, na.rm = TRUE), 2))
cat("\nCPU efficiency, quantiles:\n");  print(round(quantile(t$cpu_eff, q, na.rm = TRUE), 2))
cat(sprintf("\ntotal wall: %.0f task-hours\n", sum(t$elapsed_h, na.rm = TRUE)))
cat("\n10 largest by memory:\n")
print(head(t[order(-t$maxrss_gb), c("task", "State", "maxrss_gb", "elapsed_h", "cpu_eff")], 10),
      row.names = FALSE)

