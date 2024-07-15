source("SPAGeDi_functions.R")
source("helper_functions.R")

example_fam <- read.table("spatial_dist_example.fam")
example_bim <- read.table("spatial_dist_example.bim")
example_coords <- read.table("spatial_dist_example_coords.txt", header = TRUE)
example_geno <- as.matrix(read.table("spatial_dist_example.geno"))

# Or bed files can be imported using the BEDMatrix package
# example_geno <- as.matrix(BEDMatrix::BEDMatrix("spatial_dist_example.bed"))

plot(example_coords[, 2:3])


pairwise_spatial_dists <- calc_spatial_dist(data.frame(example_coords))

morans_i <- calc_morans_i(
  IDs = example_coords$id,
  genotypes = example_geno,
  spatial_distances = pairwise_spatial_dists
)


num_groups <- 15
distance_classes <- unique(quantile(pairwise_spatial_dists[, "dist"],
  seq(0, 1, length.out = num_groups + 1),
  na.rm = TRUE
))

correlogram_out <- create_correlogram(
  dist_classes = distance_classes,
  spatial_distances = pairwise_spatial_dists[, "dist"],
  MoransI_values = morans_i
)
plot(correlogram_out)



roussets_a <- calc_roussets_a(
  IDs = example_coords$id,
  genotypes = example_geno
)


b <- lm(roussets_a ~ log(pairwise_spatial_dists[, "dist"]))$coefficients[2]


####### Below is still work in progress and will eventually demonstrate the output
####### from these functions is comparable to SPAGeDi - can't promise all, or any of
####### it is currently fully functional

##### Create Spagedi files

assign_genotype <- function(x) {
  switch(as.character(x),
    "0" = "1,1",
    "1" = "1,2",
    "2" = "2,2",
    0
  )
}


spagedi_geno <- apply(example_geno, 2, function(x) sapply(x, assign_genotype))

cmd_file_name <- "cmd_input.txt"
output_dir <- "."
data_file_name <- "spatial_dist_example_spagedi.data"
output_file <- "spatial_dist_example_spagedi.out"


cat("// Line1: Num Individuals, Num Categories(0 if none defined", "Num Spatial Coordinates (-2 for lat & lon), Num Loci, Num Digits, Ploidy\n",
  file = data_file_name
)
line1 <- c(nrow(example_fam), 0, -2, nrow(example_bim), 2, 2)
write.table(rbind(line1), file = data_file_name, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
cat("//Line2: Num of Distance Intervals, maximal distance in each interval\n", file = data_file_name, append = TRUE)
line2 <- c(num_groups, distance_classes)
write.table(rbind(line2), file = data_file_name, append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("//Line3: Name for individuals, Name for coordinates, Name for each locus\n", file = data_file_name, append = TRUE)
line3 <- c("IND", "LAT", "LONG", example_bim$V2)
write.table(rbind(line3), sep = "\t", file = data_file_name, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("//Line4+: Individual Name, Latitude, Longitude, Genotype at locus\n", file = data_file_name, append = TRUE)

write.table(cbind(example_coords[, c("id", "latitude", "longitude")], spagedi_geno), file = data_file_name, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(rbind(c("END")), file = data_file_name, append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#### Create Command File

# Function to create command file to pass to spagedi

# ind_stats: STATISTICS for individual level analyses
# comp_opts: COMPUTATIONAL OPTIONS See Manual - 0 for None
# output_opts: OUTPUT OPTIONS - 0 for None
# spat_format: FORMAT FOR PAIRWISE SPATIAL AND GENETIC DISTANCES


# TODO: This needs fixed as option only appears if categories are provided
# LEVEL OF ANALYSES Spatial analyses carried out INDIVIDUALS (1) POPULATION (2)
# level_of_analysis <- 1


createcmd <- function(cmd_file,
                      output_file,
                      ind_stats = "A", comp_opts = 0, output_opts = 4, spat_format = 3) {
  if (comp_opts == 0) {
    comp_opts <- ""
  }

  if (output_opts == 0) {
    output_opts <- ""
  }

  # cat("\n", file=cmd_file)
  cat(data_file_name, "\n\n", file = cmd_file, sep = "")
  cat(output_file, "\n\n", file = cmd_file, append = TRUE, sep = "")


  if (file.exists(output_file)) {
    message("WARNING SPECIFIED OUTPUT FILE WILL OVERWRITE PREVIOUS\n")
    will_continue <- readline(prompt = "Press [enter] to continue. Any other key then [enter] to stop: ")
    if (will_continue != "") {
      stop("Processing stopped by user input")
    }
    file.remove(output_file)
  }


  cat("\n\n", file = cmd_file, append = TRUE, sep = "")
  # cat(level_of_analysis, '\n\n', file=cmd_file, append=TRUE)
  cat(ind_stats, "\n\n", file = cmd_file, append = TRUE, sep = "")
  cat(comp_opts, "\n", file = cmd_file, append = TRUE, sep = "")
  cat(output_opts, "\n\n", file = cmd_file, append = TRUE, sep = "")
  cat(spat_format, "\n\n\n\n", file = cmd_file, append = TRUE, sep = "")
}

# For Moran's I - note you can combine run multiple analyses together e.g. "13A",
# but this may not be possible if you don't have much RAM or are using a 32bit version


cmd_file <- paste0(output_dir, "\\", cmd_file_name)
spagedi_output_file <- paste0(output_dir, "\\", output_file)

createcmd(cmd_file,
  spagedi_output_file,
  ind_stats = "3",
  comp_opts = 0,
  output_opts = 4,
  spat_format = 3
)


#### Exact implementation may differ depending on OS - this creates a bat file to allow it work on Windows based machines
#### Assumes spagedi program is on PATH and is called as spagedi


run_spagdi <- function(cmd_file) {
  tmp_file <- create_temp_file(".bat")
  cat(paste0("spagedi < ", cmd_file), file = tmp_file)
  system2(tmp_file, stdout = "test.out")
  file.remove(tmp_file)
}

run_spagdi(cmd_file)


extract_table <- function(file, table_name, transpose = FALSE, skip_first = 0) {
  scanned_file <- scan(file, what = list(""), sep = "\n", blank.lines.skip = FALSE)

  start_line <- which(unlist(lapply(scanned_file, function(x) {
    table_name == substr(x, 1, nchar(table_name))
  })))
  blank_lines <- which(unlist(lapply(scanned_file, function(x) {
    grepl("^\\s*$", x, perl = TRUE)
  })))
  end_line <- blank_lines[blank_lines > start_line][[1]] - 1

  table_list <- scanned_file[[1]][start_line:end_line]
  table_content_lines <- unlist(lapply(table_list, function(x) {
    grepl("\\t", x)
  }))

  print("Table Info:")
  for (line in table_list[!table_content_lines]) {
    print(line, quote = FALSE)
  }

  df <- as.data.frame(do.call(rbind, strsplit(table_list[table_content_lines], "\t")))
  if (transpose) {
    df <- as.data.frame(t(df))
  }
  if (skip_first > 0) {
    df <- df[-c(1, skip_first), ]
  }
  # Remove row is all values are blank
  names(df) <- make.names(df[1, ])
  df <- df[-1, ]
  df <- df[apply(df, 1, function(x) !all(x == "")), ]
  df <- df[, apply(df, 2, function(x) !all(x == ""))]
  return(df)
}


computed_stats_df <- extract_table(output_file,
  "VALUES OF THE COMPUTED STATISTICS",
  transpose = TRUE
)

morans_i_spagedi <- extract_table(output_file,
  "Pairwise RELATIONSHIP coefficients (Moran's I for individual allele freq)",
  skip_first = 1
)
row.names(morans_i_spagedi) <- morans_i_spagedi[, 1]
morans_i_spagedi <- morans_i_spagedi[, -1]


plot(unlist(morans_i_spagedi[1, 2:16]), correlogram_out)


### For Rousset's a

createcmd(cmd_file,
  spagedi_output_file,
  ind_stats = "A",
  comp_opts = 0,
  output_opts = 4,
  spat_format = 3
)


run_spagdi(cmd_file)


roussets_a_spagedi <- extract_table(output_file,
  "Pairwise Rousset's distance between individuals ('a' in Rousset, 2000)",
  skip_first = 1
)

row.names(roussets_a_spagedi) <- roussets_a_spagedi[, 1]
roussets_a_spagedi <- roussets_a_spagedi[, -1]

cat(unname(b), " is (almost) the same as ", roussets_a_spagedi[1, "b.log.slope.log.dist."])
