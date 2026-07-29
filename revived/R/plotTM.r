args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[[1]] else file.path("runs", "results")
output_file <- if (length(args) >= 2) args[[2]] else "atka-mcmc-summary.pdf"

if (!dir.exists(results_dir)) {
  stop("Results directory does not exist: ", results_dir)
}

files <- list.files(
  results_dir,
  pattern = "_mcmc[.]rep$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No *_mcmc.rep files found in ", results_dir)
}

read_mcmc <- function(path) {
  value <- read.table(path, header = TRUE, check.names = FALSE)
  value$Model <- sub("_mcmc[.]rep$", "", basename(path))
  value
}

tables <- lapply(files, read_mcmc)
all_columns <- unique(unlist(lapply(tables, names)))
tables <- lapply(tables, function(value) {
  missing_columns <- setdiff(all_columns, names(value))
  value[missing_columns] <- NA
  value[all_columns]
})
samples <- do.call(rbind, tables)
if (!"Biomass" %in% names(samples)) {
  stop("MCMC files do not contain a Biomass column")
}

summary_table <- do.call(
  rbind,
  lapply(split(samples$Biomass, samples$Model), function(x) {
    c(
      samples = length(x),
      median = median(x, na.rm = TRUE),
      lower_95 = unname(quantile(x, 0.025, na.rm = TRUE)),
      upper_95 = unname(quantile(x, 0.975, na.rm = TRUE))
    )
  })
)
print(summary_table)

grDevices::pdf(output_file, width = 10, height = 7)
model_names <- unique(samples$Model)
colours <- grDevices::hcl.colors(length(model_names), "Dark 3")
densities <- lapply(
  model_names,
  function(model) stats::density(samples$Biomass[samples$Model == model], na.rm = TRUE)
)
x_range <- range(vapply(densities, function(x) range(x$x), numeric(2)))
y_range <- c(0, max(vapply(densities, function(x) max(x$y), numeric(1))))
plot(
  densities[[1]],
  xlim = x_range,
  ylim = y_range,
  col = colours[[1]],
  lwd = 2,
  main = "Atka tag-model posterior biomass",
  xlab = "Biomass",
  ylab = "Density"
)
if (length(densities) > 1) {
  for (index in 2:length(densities)) {
    lines(densities[[index]], col = colours[[index]], lwd = 2)
  }
}
legend("topright", legend = model_names, col = colours, lwd = 2, cex = 0.7)
grDevices::dev.off()

message("Wrote ", output_file)
