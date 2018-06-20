#!/usr/bin/env Rscript

# Original author:  Leonardo Tenori (2016)
# Refactored by Luca Pireddu (2018)

# Usage: Rscript integrals.R --left=5.257 --right=5.225 --where=5.243 --plotfile somefile.pdf . metabolites.txt
# Usage: Rscript integrals.R --left=5.257 --right=5.225 --where=5.243 --plotfile somefile.pdf metabolites.txt reference_dir test_dir

# `kind` correspnods to the subdirectory of the spectrum directory.  I saw: { 1, 3, 4, 98888 }
# Used by load_spectrum to determine from where to load the data
kind=3

library("caTools")
library("optparse") # for python-style option parsing

get_spectra_dirs <- function(data_root) {
  d <- dir(data_root, full.names=TRUE)
  d[ dir.exists(d) ]
}

parse_spectrum_parameters <- function(procs_dir) {
    procs_file <- file.path(procs_dir, "procs")
    one_r_file <- file.path(procs_dir, "1r")

    Procs <- scan(file=procs_file, what="list")

    BYTORD <- Procs[ which(Procs == "##$BYTORDP=") + 1 ]
    ENDIAN <- if (BYTORD=="0") "little" else "big"

    NC_proc <- as.numeric( Procs[ which(Procs == "##$NC_proc=") + 1 ] )
    OFFSET  <- as.numeric( Procs[ which(Procs == "##$OFFSET=")  + 1 ] )
    SW_p    <- as.numeric( Procs[ which(Procs == "##$SW_p=")    + 1 ] )
    SF      <- as.numeric( Procs[ which(Procs == "##$SF=")      + 1 ] )
    SI      <- as.numeric( Procs[ which(Procs == "##$SI=")      + 1 ] )

    con <- file(one_r_file, open="rb")
    y <- readBin(con, what="int", n=SI, endian=ENDIAN)
    close(con)
    x <- as.numeric( seq(OFFSET, OFFSET - SW_p/SF, length=SI) )
    y <- as.numeric( y/(2^(-NC_proc)) )

    retval <- list(
        x=x,
        y=y,
        NC_proc= NC_proc,
        OFFSET = OFFSET,
        SW_p   = SW_p,
        SF     = SF,
        SI     = SI
    )
    return (retval)
}

# returns a named list
load_spectrum <- function(spectrum_dir) {
    procs_dir <- file.path(spectrum_dir, as.character(kind), "pdata", "1")

    ## Read in parameters, needed for data processing
    params <- parse_spectrum_parameters(procs_dir)
    params$path <- spectrum_dir

    return (params)
}

align_spectrum <- function(spectrum_object, left, right, where) {
    # align spectra to a reference signal

    v <- which(left >= spectrum_object$x & right <= spectrum_object$x)
    s <- which.max(spectrum_object$y[v])
    reference <- spectrum_object$x[v][s]
    shift <- where - reference
    new_offset <- spectrum_object$OFFSET + shift

    return (new_offset)
}

align_all <- function(spectrum_objects, left, right, where) {
    new_offsets <- sapply(spectrum_objects, align_spectrum, left=left, right=right, where=where)

    # This function returns a named numeric array
    return (new_offsets)
}

integrate_all <- function(spectrum_objects, aligned_offsets, metabolites_table) {
    # integrate spectra
    require(caTools)
    if (length(spectrum_objects) != length(aligned_offsets)) {
        stop("number of aligned offsets doesn't match number of spectra")
    }

    inte = matrix(nrow=length(spectrum_objects), ncol=nrow(metabolites_table))
    colnames(inte) = metabolites_table[,1]
    xx = list(list(),list()); yy=list(list(),list())

    for (idx in seq_along(spectrum_objects)) {
        spectrum <- spectrum_objects[[idx]]

        # normalize
        v <- which(spectrum$x >= 5.5 | spectrum$x <= 4.35)
        v <- sum(abs(spectrum$y[v]))
        y <- 1e6*spectrum$y/v

        for (ff in 1:nrow(metabolites_table)) {

            left = metabolites_table[ff,2]; right = metabolites_table[ff,3]; where=metabolites_table[ff,4];
            l1 = left + 0.0005;  l2 = left + 0.0001;
            r1 = right - 0.0001; r2 = right - 0.0005;

            OFFSET <- aligned_offsets[idx]
            if (metabolites_table[ff,6]==1){

                ## align
                v <- which(left >= spectrum$x & right <= spectrum$x)
                s <- which.max(y[v])
                reference <- spectrum$x[v][s]
                shift <- where - reference
                OFFSET <- aligned_offsets[idx] + shift
            }

            if (metabolites_table[ff,5]==1){
                ### raddrizza
                x <- as.numeric(seq(OFFSET, OFFSET - spectrum$SW_p/spectrum$SF, length=spectrum$SI))
                v <- which(left>=x & right<=x)

                A = min(as.numeric( y[which(l1 >= x & l2 <= x)] ))
                B = min(as.numeric( y[which(r1 >= x & r2 <= x)] ))

                left1 = mean(l1, l2)
                right1 = mean(r1, r2)
                v = which(left1 >= x & right1 <= x)
                xxxx = c(left1, right1)
                yyyy = c(A, B)
                f = coef(lm(yyyy ~ xxxx))
                baseline = as.numeric(f[1]) + x[v]*as.numeric(f[2])
                y[v] = y[v] - baseline
            }
            ## integral
            x <- as.numeric(seq (OFFSET, OFFSET - spectrum$SW_p/spectrum$SF, length=spectrum$SI))
            v <- which(left >= x & right <= x)

            inte[idx,ff] <- trapz(rev(x[v]), y[v])

            xx[[idx]][[ff]] = rev(x[v])
            yy[[idx]][[ff]] = rev(y[v])
            x <- as.numeric( seq(aligned_offsets[idx], aligned_offsets[idx] - spectrum$SW_p/spectrum$SF, length=spectrum$SI) )
        }
    }

    return (list(inte=inte, x=xx, y=yy))
}

print_report <- function(integration, metabolites_table) {
    for (h in 1:nrow(metabolites_table) ) {
        ss = paste(metabolites_table[h,1], "ok")
        pp = paste(metabolites_table[h,1], "bad")

        reference = integration[1,h]
        test = integration[2,h]

        tt = test * 100 / reference
        ttt = abs(100 - tt)
        output <- ifelse(ttt < 20, ss, pp)
        cat(output, "\tdiff\t", round(ttt,2), "\t%\n")
    }
}


plot_metabolites <- function(metabolites, yy, xx) {
    for (h in 1:nrow(metabolites) ) {
        left = metabolites[h,2]; right = metabolites[h,3]; where=metabolites[h,4]

        kk = c(max(yy[[1]][[h]]), max(yy[[2]][[h]]))
        dd = length(xx[[1]][[h]])

        plot(xx[[1]][[h]], rev(yy[[1]][[h]]),
             ylim=c(0,(max(kk) + max(kk) / 5)),
             type="l", xaxt="n",
             xlab="ppm", ylab="Relative intensity")

        axis(1,
             at=c(xx[[1]][[h]][1], xx[[1]][[h]][dd]),
             labels=c(round(rev(xx[[1]][[h]])[1], 3), round(rev(xx[[1]][[h]])[dd], 3)),
             las=1, cex.axis=1)
        points(xx[[2]][[h]],rev(yy[[2]][[h]]), type="l",col=2)

        abline(v=left, lty=2, col="blue")
        abline(v=right, lty=2, col="blue")

        title(metabolites[h,1])
    }
}


do_arg_parsing <- function() {
    option_list <- list(
        make_option("--left", type="numeric", default=5.257, help="Left", metavar="N"),
        make_option("--right", type="numeric", default=5.225, help="Right", metavar="N"),
        make_option("--where", type="numeric", default=5.243, help="Where", metavar="N"),
        make_option("--plotfile", type="character", default="plot.pdf", help="Plot file", metavar="FILE.pdf")
        )
    epilogue_help <- paste("\tMETABOLITES_TABLE\n\t\tMetabolites table in tab-separated format\n",
													 "\tREFERENCE_SPECTRUM\n\t\tReference dataset (Bruker NMR)\n\n",
													 "\tTEST_SPECTRUM\n\t\tTest dataset (Bruker NMR)\n\n"
													)

    parser <- OptionParser(
                  usage = "%prog [options] SPECTRA_DIR METABOLITES_TABLE",
                  option_list=option_list,
                  epilogue=epilogue_help)

    args <- parse_args(parser, positional_arguments = 3)

    metab_file <- args$args[1]
    reference_dir <- args$args[2]
    test_dir <- args$args[3]
    if (!dir.exists(reference_dir) || file.access(reference_dir, 1 | 4) < 0) { # require r-x perms
        stop(sprintf("The reference dataset directory %s either doesn't exist or we don't have access permission", reference_dir))
    }
		else if (!dir.exists(test_dir) || file.access(test_dir, 1 | 4) < 0) { # require r-x perms
        stop(sprintf("The reference dataset directory %s either doesn't exist or we don't have access permission", test_dir))
    }
    else if (file.access(metab_file, 4) < 0) {
        stop(sprintf("Can't read metabolites table file %s", metab_file))
    }

    args$reference_dir <- reference_dir
    args$test_dir <- test_dir
    args$metabolites_file <- metab_file
    return (args)
}

main <- function(path_to_spectra) {
    args <- do_arg_parsing()

    message(sprintf("Loading metabolites table from %s", args$metabolites_file))
    metabolites_t <- read.table(args$metabolites_file)

		# order is important.  Later functions assume the reference is the first in the list
    spectra_dirs <- c(args$reference_dir, args$test_dir)

    message(sprintf("Loading %d spectra...", length(spectra_dirs)))
    spectrum_objects <- lapply(spectra_dirs, load_spectrum)
    message(paste("Loaded", length(spectrum_objects), "spectrum objects"))

    message("Aligning")
    new_offsets <- align_all(spectrum_objects, args$options$left, args$options$right, args$options$where)
    message("Integrating")
    integr_results <- integrate_all(spectrum_objects, new_offsets, metabolites_t)
    message("Done.  Printing output")

    print_report(integr_results$inte, metabolites_t)

    pdf(args$options$plotfile)
    plot_metabolites(metabolites_t, integr_results$y, integr_results$x)
    ignored_ret_device <- dev.off()
}

main(spectra)
