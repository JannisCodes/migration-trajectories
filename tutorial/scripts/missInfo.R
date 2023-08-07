#' Extract missing information from datasets
#'
#' This function compares the full and reduced datasets and returns a summary table of the missing data.
#'
#' @param full A data frame for the full dataset.
#' @param reduced A data frame for the reduced dataset.
#' @return A data frame containing information about:
#'   * Total number of rows (nFull, nRed, nDif, nDifPerc).
#'   * Unique participants (pptFull, pptRed, pptDif, pptDifPerc).
#'   * Unique time points (timeFull, timeRed, timeDif, timeDifPerc).
#' @examples
#' # Create example datasets
#' full_data <- data.frame(PID = c(1,2,2,3,3,3), TID = c(1,1,2,1,2,3))
#' reduced_data <- data.frame(PID = c(1,2,3), TID = c(1,1,2))
#' result <- missInfo(full = full_data, reduced = reduced_data)
#' @ImportFrom dplyr mutate select %>%

missInfo <- function(full = NA, reduced = NA){
  missTab <- data.frame(
    nFull = nrow(full),
    nRed = nrow(reduced),
    pptFull = length(unique(full$PID)),
    pptRed = length(unique(reduced$PID)),
    timeFull = length(unique(full$TID)),
    timeRed = length(unique(reduced$TID))
  ) %>%
    mutate(
      nDif = nFull-nRed,
      nDifPerc = nDif/nFull*100,
      pptDif = pptFull-pptRed,
      pptDifPerc = pptDif/pptFull*100,
      timeDif = timeFull-timeRed,
      timeDifPerc = timeDif/timeFull*100
    ) %>%
    select(
      nFull, nRed, nDif, nDifPerc,
      pptFull, pptRed, pptDif, pptDifPerc,
      timeFull, timeRed, timeDif, timeDifPerc
    )
  missTab
}