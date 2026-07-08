# Download and prepare daily OSAP value-weighted quintile portfolio returns.

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
raw_dir <- file.path(repo_dir, "data-raw", "empirics", "osap")
output_dir <- file.path(repo_dir, "data", "empirics", "osap")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

refresh <- "--refresh" %in% commandArgs(TRUE)

osap_release <- "2025.10"
osap_data_url <- "https://www.openassetpricing.com/data/"
gdrive_file_id <- "1JAxLf5d7THy9av9_IliGt8I4vR-zkwGo"
gdrive_url <- paste0(
  "https://drive.google.com/uc?export=download&id=",
  gdrive_file_id
)

zip_file <- file.path(raw_dir, "CtsPredictorQuintileVW.zip")
rds_file <- file.path(output_dir, "osap_daily_quintile_vw.rds")
returns_csv_file <- file.path(output_dir,
                              "osap_daily_quintile_vw_returns.csv")
losses_csv_file <- file.path(output_dir,
                             "osap_daily_quintile_vw_losses.csv")
metadata_csv_file <- file.path(output_dir,
                               "osap_daily_quintile_vw_metadata.csv")

signal_specs <- data.frame(
  signal = c("BidAskSpread", "Illiquidity", "DolVol", "zerotrade1M",
             "zerotrade6M", "Price"),
  label = c("Bid-ask spread", "Illiquidity", "Dollar volume",
            "Zero-trading frequency, 1 month",
            "Zero-trading frequency, 6 months", "Price"),
  category = c("liquidity", "liquidity", "trading activity",
               "trading activity", "trading activity", "price level"),
  stringsAsFactors = FALSE
)

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
    any(grepl("^CtsPredictorQuintileVW/.+_ret[.]csv$", listing$Name))
  }, error = function(e) FALSE)
  isTRUE(ok)
}

if (!zip_file_is_valid(zip_file) || refresh) {
  message("Downloading OSAP daily value-weighted quintile portfolios.")
  invisible(download_file(gdrive_url, zip_file))
}

if (!zip_file_is_valid(zip_file)) {
  stop(
    "The downloaded OSAP file is not a valid zip archive. ",
    "Google Drive may have returned a confirmation page. ",
    "Download DailyPortfolios/CtsPredictorQuintileVW.zip from ",
    osap_data_url, " and place it at: ", zip_file,
    call. = FALSE
  )
}

read_signal_returns <- function(signal) {
  member <- file.path("CtsPredictorQuintileVW", paste0(signal, "_ret.csv"))
  table <- read.csv(
    unz(zip_file, member),
    header = TRUE,
    check.names = FALSE,
    strip.white = TRUE,
    na.strings = c("", "NA", "-99.99", "-999", "-999.00")
  )
  if (!"date" %in% names(table)) {
    stop("Missing date column in OSAP file member: ", member, call. = FALSE)
  }
  quintile_columns <- paste0("port0", 1:5)
  missing_columns <- setdiff(quintile_columns, names(table))
  if (length(missing_columns) > 0L) {
    stop("Missing quintile columns in OSAP file member ", member, ": ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  table$date <- as.Date(table$date)
  if (any(is.na(table$date))) {
    stop("Could not parse dates in OSAP file member: ", member,
         call. = FALSE)
  }
  for (column in quintile_columns) {
    table[[column]] <- as.numeric(table[[column]])
  }

  out <- table[, c("date", quintile_columns)]
  names(out) <- c(
    "date",
    paste(tolower(signal), paste0("q", 1:5), sep = "_")
  )
  out <- out[order(out$date), ]
  rownames(out) <- NULL
  attr(out, "raw_portfolio_names") <- quintile_columns
  out
}

tables <- vector("list", nrow(signal_specs))
metadata <- vector("list", nrow(signal_specs))

for (i in seq_len(nrow(signal_specs))) {
  spec <- signal_specs[i, ]
  table <- read_signal_returns(spec$signal)
  columns <- setdiff(names(table), "date")
  raw_names <- attr(table, "raw_portfolio_names")
  tables[[i]] <- table
  metadata[[i]] <- data.frame(
    variable = columns,
    signal = spec$signal,
    signal_label = spec$label,
    category = spec$category,
    portfolio = raw_names,
    source = paste0("OSAP DailyPortfolios/CtsPredictorQuintileVW.zip, ",
                    "release ", osap_release),
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
  signal_specs = signal_specs,
  source = list(
    provider = "Open Source Asset Pricing",
    data_set = "Daily value-weighted continuous-predictor quintile portfolios",
    release = osap_release,
    units = "percent daily return; losses are negative returns",
    url = osap_data_url,
    google_drive_file_id = gdrive_file_id,
    downloaded_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
)

saveRDS(prepared, rds_file)
write.csv(returns, returns_csv_file, row.names = FALSE)
write.csv(losses, losses_csv_file, row.names = FALSE)
write.csv(metadata, metadata_csv_file, row.names = FALSE)

message("Saved prepared OSAP daily value-weighted quintile data:")
message("  ", rds_file)
message("  ", returns_csv_file)
message("  ", losses_csv_file)
message("  ", metadata_csv_file)
message("Complete daily panel: ", nrow(complete_losses), " observations from ",
        min(complete_losses$date), " to ", max(complete_losses$date))
message("Portfolio variables: ", length(portfolio_columns))
