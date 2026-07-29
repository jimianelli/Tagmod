args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath(".")
archive_dir <- file.path(project_dir, "runs", "results")
generated_dir <- file.path(project_dir, "runs", "generated")
output_csv <- file.path(generated_dir, "mode-validation.csv")
output_md <- file.path(generated_dir, "mode-validation.md")

read_par <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header <- lines[[1]]
  objective <- as.numeric(sub(
    ".*Objective function value = ([^ ]+).*",
    "\\1",
    header
  ))
  gradient <- as.numeric(sub(
    ".*Maximum gradient component = ([^ ]+).*",
    "\\1",
    header
  ))
  labels <- grep("^# .*:$", lines)
  values <- vapply(labels, function(index) {
    as.numeric(strsplit(trimws(lines[[index + 1]]), "[[:space:]]+")[[1]][[1]])
  }, numeric(1))
  names(values) <- sub("^# (.*):$", "\\1", lines[labels])
  list(objective = objective, gradient = gradient, parameters = values)
}

archive_files <- list.files(archive_dir, pattern = "[.]par$", full.names = TRUE)
archive_names <- sub("[.]par$", "", basename(archive_files))

rows <- lapply(seq_along(archive_files), function(index) {
  model <- archive_names[[index]]
  candidates <- list.dirs(generated_dir, recursive = FALSE, full.names = TRUE)
  pattern <- paste0("^", tolower(model), "_mpd_[0-9]{8}t[0-9]{6}z$")
  candidates <- candidates[
    grepl(pattern, tolower(basename(candidates)))
  ]
  if (length(candidates) == 0) {
    return(data.frame(
      model = model,
      status = "missing generated fit",
      objective_difference = NA_real_,
      maximum_parameter_difference = NA_real_,
      archived_gradient = NA_real_,
      generated_gradient = NA_real_
    ))
  }

  generated_path <- file.path(sort(candidates, decreasing = TRUE)[[1]], "tm.par")
  archived <- read_par(archive_files[[index]])
  generated <- read_par(generated_path)
  common <- intersect(names(archived$parameters), names(generated$parameters))
  parameter_difference <- if (length(common)) {
    max(abs(
      archived$parameters[common] - generated$parameters[common]
    ))
  } else {
    NA_real_
  }
  objective_difference <- abs(archived$objective - generated$objective)
  tolerance <- 1e-8
  status <- if (
    is.finite(objective_difference) &&
    objective_difference <= tolerance &&
    is.finite(parameter_difference) &&
    parameter_difference <= tolerance
  ) {
    "reproduced"
  } else {
    "different"
  }

  data.frame(
    model = model,
    status = status,
    objective_difference = objective_difference,
    maximum_parameter_difference = parameter_difference,
    archived_gradient = archived$gradient,
    generated_gradient = generated$gradient
  )
})

validation <- do.call(rbind, rows)
validation <- validation[order(tolower(validation$model)), ]
write.csv(validation, output_csv, row.names = FALSE)

format_number <- function(x) {
  ifelse(is.na(x), "", format(x, scientific = TRUE, digits = 4))
}

markdown <- c(
  "# Archived mode-fit validation",
  "",
  paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "",
  "| Model | Status | Objective difference | Maximum parameter difference | Archived max gradient | Generated max gradient |",
  "|---|---:|---:|---:|---:|---:|",
  apply(validation, 1, function(row) {
    paste0(
      "| ", row[["model"]],
      " | ", row[["status"]],
      " | ", format_number(as.numeric(row[["objective_difference"]])),
      " | ", format_number(as.numeric(row[["maximum_parameter_difference"]])),
      " | ", format_number(as.numeric(row[["archived_gradient"]])),
      " | ", format_number(as.numeric(row[["generated_gradient"]])),
      " |"
    )
  })
)
writeLines(markdown, output_md)

print(validation, row.names = FALSE)
message("Wrote ", output_csv)
message("Wrote ", output_md)

if (any(validation$status != "reproduced")) {
  quit(status = 1)
}
