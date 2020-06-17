library(stringr)

assert <- function(cond, msg) {
    if (!cond)
        stop(paste("assertion failed:", msg))
}

parseArgs <- function() {
    args <- commandArgs(trailingOnly=TRUE)

    if (length(args) != 2) {
      cat("usage: DATA_DIR PLOT_DIR\n")
      quit()
    }

    list("data_dir" = args[1], "plot_dir" = args[2])
}

readData_ <- function(data_dir, columns) {
    d_files <- list.files(path=data_dir, pattern=".txt", full.names=TRUE)
    assert(length(d_files) > 0, "data directory not empty")

    d_all = NULL

    invisible(lapply(d_files, function(d_file) {
        d_file_pattern =
          paste(".*/",
                paste(rep("(.*?)_", length(columns) - 1), collapse=''),
                "(.*).txt",
                sep='')

        # parse file name
        m <- str_match(d_file, d_file_pattern)
        assert(length(m) == length(columns) + 1, "valid data file name")

        # read data into data frame
        d <- read.csv(d_file)

        dim_d = dim(d)
        assert(length(dim_d) == 2 && dim_d[2], "data dimensions correct")

        d$runtime = as.numeric(d$runtime)

        for (i in 1:length(columns))
          d[[columns[i]]] = m[i + 1]

        # concatenate to main data frame
        if (is.null(d_all))
          d_all <<- d
        else
          d_all <<- rbind(d_all, d)
    }))

    d_all
}

readData <- function(data_dir, columns) {
    cat(paste("reading data from directory '", data_dir, "'\n", sep=''))
    d <- readData_(data_dir, columns)
}
