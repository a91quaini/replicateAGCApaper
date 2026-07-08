# Author: Alberto Quaini

#' Run the portfolio empirical scripts
#'
#' @param scripts Character vector selecting scripts. Use `"all"` to run both
#'   empirical workflows.
#' @return Invisibly returns `TRUE` when scripts complete.
#' @export
run_empirics <- function(scripts = c("all", "agca_reporting", "tail_functionals")) {
  scripts <- match.arg(scripts, several.ok = TRUE)
  if ("all" %in% scripts) {
    scripts <- c("tail_functionals", "agca_reporting")
  }
  files <- c(
    tail_functionals = replication_file(
      "inst", "empirics", "scripts", "portfolio_tail_functionals_agca.R"
    ),
    agca_reporting = replication_file(
      "inst", "empirics", "scripts", "portfolio_agca_reporting_figures.R"
    )
  )
  lapply(files[scripts], source, local = new.env(parent = globalenv()), chdir = TRUE)
  invisible(TRUE)
}

#' Run the simulation scripts
#'
#' @param scripts Character vector selecting scripts. Use `"all"` to run both
#'   simulation workflows.
#' @return Invisibly returns `TRUE` when scripts complete.
#' @export
run_simulations <- function(scripts = c("all", "simulation_3d", "simulation_10d")) {
  scripts <- match.arg(scripts, several.ok = TRUE)
  if ("all" %in% scripts) {
    scripts <- c("simulation_3d", "simulation_10d")
  }
  files <- c(
    simulation_3d = replication_file(
      "inst", "simulations", "scripts", "standard_simulations_3d.R"
    ),
    simulation_10d = replication_file(
      "inst", "simulations", "scripts", "standard_simulations.R"
    )
  )
  lapply(files[scripts], source, local = new.env(parent = globalenv()), chdir = TRUE)
  invisible(TRUE)
}

#' Run all replication workflows
#'
#' This runs the simulation and empirical scripts with their current defaults.
#' The full workflow can be computationally expensive because it includes
#' bootstrap and finite-sample simulation loops.
#'
#' @return Invisibly returns `TRUE` when scripts complete.
#' @export
run_all <- function() {
  run_simulations()
  run_empirics()
  invisible(TRUE)
}
