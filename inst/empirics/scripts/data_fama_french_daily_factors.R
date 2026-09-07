# Download and prepare the daily US Fama--French five factors and momentum.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE))
} else {
  "."
}

find_compendium_root <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  repeat {
    description <- file.path(path, "DESCRIPTION")
    if (file.exists(description) &&
        grepl("^Package: replicateAGCApaper",
              readLines(description, n = 1L))) {
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

archive_specs <- data.frame(
  id = c("five_factor", "momentum"),
  archive = c(
    "F-F_Research_Data_5_Factors_2x3_daily_CSV.zip",
    "F-F_Momentum_Factor_daily_CSV.zip"
  ),
  url = c(
    paste0(
      "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/",
      "F-F_Research_Data_5_Factors_2x3_daily_CSV.zip"
    ),
    paste0(
      "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/",
      "F-F_Momentum_Factor_daily_CSV.zip"
    )
  ),
  stringsAsFactors = FALSE
)

download_archive <- function(url, destination) {
  curl <- Sys.which("curl")
  if (curl == "") {
    stop("The system command 'curl' is required to download factor data.",
         call. = FALSE)
  }
  temporary <- paste0(destination, ".tmp")
  status <- system2(
    curl,
    args = c("-L", "--fail", "--silent", "--show-error", url,
             "-o", temporary)
  )
  if (!identical(status, 0L) || !file.exists(temporary) ||
      file.info(temporary)$size == 0L) {
    if (file.exists(temporary)) {
      unlink(temporary)
    }
    stop("Could not download factor archive: ", url, call. = FALSE)
  }
  if (file.exists(destination)) {
    unlink(destination)
  }
  if (!file.rename(temporary, destination)) {
    stop("Could not move factor archive to: ", destination, call. = FALSE)
  }
  invisible(destination)
}

archive_is_valid <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0L) {
    return(FALSE)
  }
  tryCatch({
    listing <- unzip(path, list = TRUE)
    any(grepl("[.]csv$", listing$Name, ignore.case = TRUE))
  }, error = function(e) FALSE)
}

for (i in seq_len(nrow(archive_specs))) {
  destination <- file.path(raw_dir, archive_specs$archive[i])
  if (!archive_is_valid(destination)) {
    message("Downloading ", archive_specs$id[i], " daily factor archive.")
    download_archive(archive_specs$url[i], destination)
  }
}

read_factor_archive <- function(path, header_pattern) {
  listing <- unzip(path, list = TRUE)
  member <- listing$Name[grepl("[.]csv$", listing$Name,
                              ignore.case = TRUE)][1L]
  if (is.na(member)) {
    stop("No CSV member found in: ", path, call. = FALSE)
  }
  lines <- readLines(unz(path, member), warn = FALSE)
  header <- grep(header_pattern, lines)[1L]
  if (is.na(header)) {
    stop("Could not locate factor header in: ", path, call. = FALSE)
  }
  data_rows <- grep("^[0-9]{8},", lines)
  data_rows <- data_rows[data_rows > header]
  if (length(data_rows) == 0L) {
    stop("Could not locate daily factor rows in: ", path, call. = FALSE)
  }
  parsed <- read.csv(
    text = paste(c(lines[header], lines[data_rows]), collapse = "\n"),
    check.names = FALSE,
    na.strings = c("-99.99", "-999", "-999.00"),
    strip.white = TRUE
  )
  names(parsed)[1L] <- "yyyymmdd"
  parsed$yyyymmdd <- sprintf("%08d", as.integer(parsed$yyyymmdd))
  parsed$date <- as.Date(parsed$yyyymmdd, format = "%Y%m%d")
  if (anyNA(parsed$date)) {
    stop("One or more factor dates could not be parsed in: ", path,
         call. = FALSE)
  }
  parsed$yyyymmdd <- NULL
  parsed[, c("date", setdiff(names(parsed), "date")), drop = FALSE]
}

five <- read_factor_archive(
  file.path(raw_dir, archive_specs$archive[archive_specs$id == "five_factor"]),
  "^,Mkt-RF,SMB,HML,RMW,CMA,RF"
)
momentum <- read_factor_archive(
  file.path(raw_dir, archive_specs$archive[archive_specs$id == "momentum"]),
  "^,Mom"
)

names(five)[names(five) == "Mkt-RF"] <- "MKT_RF"
names(momentum)[names(momentum) == "Mom"] <- "MOM"
factors <- merge(five, momentum, by = "date", all = FALSE)
factors <- factors[order(factors$date), ]
rownames(factors) <- NULL

numeric_columns <- setdiff(names(factors), "date")
for (column in numeric_columns) {
  factors[[column]] <- as.numeric(factors[[column]])
}
if (anyNA(factors)) {
  stop("Merged daily factor panel contains missing values.", call. = FALSE)
}

prepared <- list(
  factors = factors,
  units = "percent daily return",
  candidate_factors = c("SMB", "HML", "RMW", "CMA", "MOM"),
  nuisance_factors = "MKT_RF",
  source = list(
    provider = "Kenneth R. French Data Library",
    downloaded_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    urls = setNames(archive_specs$url, archive_specs$id)
  )
)

rds_file <- file.path(output_dir, "ff_six_factors_daily.rds")
csv_file <- file.path(output_dir, "ff_six_factors_daily.csv")
saveRDS(prepared, rds_file)
write.csv(factors, csv_file, row.names = FALSE)

message("Saved daily factor data:")
message("  ", rds_file)
message("  ", csv_file)
message("Factor panel: ", nrow(factors), " observations from ",
        min(factors$date), " through ", max(factors$date), ".")
