# Download and prepare daily Fama-French 2x3 bivariate portfolio returns.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE))
} else {
  "."
}

find_compendium_root <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        grepl("^Package: replicateAGCApaper", readLines(file.path(path, "DESCRIPTION"), n = 1L))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the replicateAGCApaper root.", call. = FALSE)
    }
    path <- parent
  }
}

repo_dir <- find_compendium_root(script_dir)
raw_dir <- file.path(repo_dir, "data-raw", "empirics", "ff")
output_dir <- file.path(repo_dir, "data", "empirics", "ff")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

refresh <- "--refresh" %in% commandArgs(TRUE)

sort_specs <- data.frame(
  sort = c("size_bm", "size_op", "size_inv", "size_mom"),
  label = c("Size and Book-to-Market", "Size and Operating Profitability",
            "Size and Investment", "Size and Momentum"),
  url = c(
    "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/6_Portfolios_2x3_daily_CSV.zip",
    "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/6_Portfolios_ME_OP_2x3_daily_CSV.zip",
    "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/6_Portfolios_ME_INV_2x3_daily_CSV.zip",
    "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/6_Portfolios_ME_Prior_12_2_Daily_CSV.zip"
  ),
  stringsAsFactors = FALSE
)

rds_file <- file.path(output_dir, "ff_2x3_sorts_daily.rds")
returns_csv_file <- file.path(output_dir, "ff_2x3_sorts_daily_returns.csv")
losses_csv_file <- file.path(output_dir, "ff_2x3_sorts_daily_losses.csv")
metadata_csv_file <- file.path(output_dir, "ff_2x3_sorts_daily_metadata.csv")

download_file <- function(url, destfile) {
  curl <- Sys.which("curl")
  if (curl == "") {
    stop("The system command 'curl' is required to download the data.",
         call. = FALSE)
  }

  tmp_file <- paste0(destfile, ".tmp")
  tmp_stderr <- paste0(destfile, ".stderr")
  if (file.exists(tmp_file)) {
    unlink(tmp_file)
  }
  if (file.exists(tmp_stderr)) {
    unlink(tmp_stderr)
  }

  result <- system2(
    curl,
    args = c("-L", "-sS", "--fail", shQuote(url)),
    stdout = tmp_file,
    stderr = tmp_stderr
  )
  status <- attr(result, "status")
  if (is.null(status)) {
    status <- result
  }
  curl_stderr <- if (file.exists(tmp_stderr)) {
    paste(readLines(tmp_stderr, warn = FALSE), collapse = "\n")
  } else {
    ""
  }

  if (status != 0L || !file.exists(tmp_file) || file.info(tmp_file)$size == 0L) {
    if (file.exists(tmp_file)) {
      unlink(tmp_file)
    }
    if (file.exists(tmp_stderr)) {
      unlink(tmp_stderr)
    }
    stop(
      "Download failed. curl status: ", status, "\n",
      curl_stderr,
      call. = FALSE
    )
  }

  if (file.exists(tmp_stderr)) {
    unlink(tmp_stderr)
  }
  if (file.exists(destfile)) {
    unlink(destfile)
  }
  if (!file.rename(tmp_file, destfile)) {
    stop("Could not move downloaded file to: ", destfile, call. = FALSE)
  }

  destfile
}

zip_file_is_valid <- function(file) {
  if (!file.exists(file) || file.info(file)$size == 0L) {
    return(FALSE)
  }
  ok <- tryCatch({
    listing <- unzip(file, list = TRUE)
    any(grepl("[.]csv$", listing$Name, ignore.case = TRUE))
  }, error = function(e) FALSE)
  isTRUE(ok)
}

read_zip_csv_lines <- function(zip_file) {
  listing <- unzip(zip_file, list = TRUE)
  csv_name <- listing$Name[grepl("[.]csv$", listing$Name, ignore.case = TRUE)][1L]
  if (is.na(csv_name)) {
    stop("The Fama-French zip archive does not contain a CSV file.",
         call. = FALSE)
  }
  readLines(unz(zip_file, csv_name), warn = FALSE)
}

parse_daily_value_weighted <- function(zip_file, prefix) {
  lines <- read_zip_csv_lines(zip_file)
  marker <- grep("Average Value Weighted Returns --", lines,
                 fixed = TRUE)
  if (length(marker) == 0L) {
    stop("Could not find the value-weighted returns block in ",
         zip_file, ".", call. = FALSE)
  }

  header_idx <- marker[1L] + 1L
  while (header_idx <= length(lines) && !grepl(",", lines[header_idx],
                                              fixed = TRUE)) {
    header_idx <- header_idx + 1L
  }
  if (header_idx > length(lines)) {
    stop("Could not find the portfolio header row in ", zip_file, ".",
         call. = FALSE)
  }

  data_start <- header_idx + 1L
  while (data_start <= length(lines) &&
         !grepl("^\\s*[0-9]{8}\\s*,", lines[data_start])) {
    data_start <- data_start + 1L
  }
  if (data_start > length(lines)) {
    stop("Could not find daily return rows in ", zip_file, ".",
         call. = FALSE)
  }

  data_end <- data_start
  while (data_end <= length(lines) &&
         grepl("^\\s*[0-9]{8}\\s*,", lines[data_end])) {
    data_end <- data_end + 1L
  }
  data_end <- data_end - 1L

  daily_lines <- c(lines[header_idx], lines[data_start:data_end])
  daily <- read.csv(
    text = paste(daily_lines, collapse = "\n"),
    header = TRUE,
    check.names = FALSE,
    strip.white = TRUE,
    na.strings = c("", "-99.99", "-999", "-999.00")
  )
  names(daily)[1L] <- "yyyymmdd"
  names(daily) <- trimws(names(daily))

  daily$yyyymmdd <- as.character(daily$yyyymmdd)
  daily$date <- as.Date(daily$yyyymmdd, format = "%Y%m%d")
  if (any(is.na(daily$date))) {
    stop("Could not parse one or more daily dates in ", zip_file, ".",
         call. = FALSE)
  }

  portfolio_columns <- setdiff(names(daily), c("yyyymmdd", "date"))
  for (column in portfolio_columns) {
    daily[[column]] <- as.numeric(daily[[column]])
  }
  names(daily)[match(portfolio_columns, names(daily))] <-
    paste(prefix, make.names(portfolio_columns), sep = "_")

  daily <- daily[, c("date", paste(prefix, make.names(portfolio_columns),
                                   sep = "_"))]
  daily <- daily[order(daily$date), ]
  rownames(daily) <- NULL
  attr(daily, "raw_portfolio_names") <- portfolio_columns
  daily
}

tables <- vector("list", nrow(sort_specs))
metadata <- vector("list", nrow(sort_specs))

for (i in seq_len(nrow(sort_specs))) {
  spec <- sort_specs[i, ]
  zip_file <- file.path(raw_dir, paste0(basename(spec$url)))
  if (!zip_file_is_valid(zip_file) || refresh) {
    message("Downloading daily Fama-French 2x3 sort: ", spec$label)
    invisible(download_file(spec$url, zip_file))
  } else {
    message("Using existing daily Fama-French zip file: ", zip_file)
  }

  table <- parse_daily_value_weighted(zip_file, spec$sort)
  columns <- setdiff(names(table), "date")
  raw_names <- attr(table, "raw_portfolio_names")
  tables[[i]] <- table
  metadata[[i]] <- data.frame(
    variable = columns,
    sort = spec$sort,
    sort_label = spec$label,
    portfolio = raw_names,
    source_url = spec$url,
    stringsAsFactors = FALSE
  )
}

returns <- Reduce(function(x, y) merge(x, y, by = "date", all = TRUE), tables)
portfolio_columns <- setdiff(names(returns), "date")
losses <- returns
losses[, portfolio_columns] <- -losses[, portfolio_columns, drop = FALSE]

complete_returns <- returns[complete.cases(returns[, portfolio_columns]), ]
complete_losses <- losses[complete.cases(losses[, portfolio_columns]), ]
metadata <- do.call(rbind, metadata)

prepared <- list(
  returns = returns,
  losses = losses,
  complete_returns = complete_returns,
  complete_losses = complete_losses,
  portfolios = metadata,
  portfolio_columns = portfolio_columns,
  source = list(
    provider = "Kenneth R. French Data Library",
    data_set = paste(sort_specs$label, collapse = "; "),
    units = "percent daily return; losses are negative returns",
    downloaded_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    urls = sort_specs$url
  )
)

saveRDS(prepared, rds_file)
write.csv(returns, returns_csv_file, row.names = FALSE)
write.csv(losses, losses_csv_file, row.names = FALSE)
write.csv(metadata, metadata_csv_file, row.names = FALSE)

message("Saved prepared daily Fama-French 2x3 bivariate-sort data:")
message("  ", rds_file)
message("  ", returns_csv_file)
message("  ", losses_csv_file)
message("  ", metadata_csv_file)
message("Complete daily panel: ", nrow(complete_losses), " observations from ",
        min(complete_losses$date), " to ", max(complete_losses$date))
message("Portfolio variables: ", length(portfolio_columns))
