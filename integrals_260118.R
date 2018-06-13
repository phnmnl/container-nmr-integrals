
# Original author:  Leonardo Tenori (2016)
# Refactored by Luca Pireddu (2018)

#datadir = "G:\\CONTROLLO QUALITA\\Long_Term_Storage_SIERO\\nmr" #percorso di dove si salvano gli spettri nmr

# `kind` corrisponde ad una sottocartella dello spettro.  Ho visto { 1, 3, 4, 98888 }
# Used by load_spectrum to determine from where to load the data
kind=3

library("caTools")

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


spectra <- "/home/pireddu/Projects/phenomenal/cerm/Script/"
metabolites <- "/home/pireddu/Projects/phenomenal/cerm/Script/metabolites.txt"

main <- function(path_to_spectra) {
    left <- 5.257
    right <- 5.225
    where <- 5.243

    # percorso del file metabolites.txt
    metabolites_t <- read.table(metabolites)

    spectra_dirs <- get_spectra_dirs(path_to_spectra)
    message("Loading spectra...")
    spectrum_objects <- lapply(spectra_dirs, load_spectrum)
    message(paste("Loaded", length(spectrum_objects), "spectrum objects"))

    new_offsets <- align_all(spectrum_objects, left, right, where)
    integr_results <- integrate_all(spectrum_objects, new_offsets, metabolites_t)

    print_report(integr_results$inte, metabolites_t)

    pdf("plot.pdf")
    plot_metabolites(metabolites_t, integr_results$y, integr_results$x)
    ignored_ret_device <- dev.off()
}

main(spectra)
