
# Original author:  Leonardo Tenori (2016)
# Refactored by Luca Pireddu (2018)

#datadir = "G:\\CONTROLLO QUALITA\\Long_Term_Storage_SIERO\\nmr" #percorso di dove si salvano gli spettri nmr

# `kind` corrisponde ad una sottocartella dello spettro.  Ho visto { 1, 3, 4, 98888 }
kind=3

library("caTools")

get_spectra_dirs <- function(data_root) list.dirs(data_root, full.names=TRUE, recursive=FALSE)

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


align_spectrum <- function(spectrum_dir, kind, left, right, where) {
    # align spectra to a reference signal

    procs_dir <- file.path(spectrum_dir, as.character(kind), "pdata", "1")

    ## Read in parameters, needed for data processing
    params <- parse_spectrum_parameters(procs_dir)
    
    ## align
    v <- which(left >= params$x & right <= params$x)
    s <- which.max(params$y[v])
    reference <- params$x[v][s]
    shift <- where - reference
    OFFSETNEW <- params$OFFSET + shift

    ## Ovewrite the procs file setting the new offset
    procs <- file.path(procs_dir, "procs")
    x <- readLines(procs)
    y <- gsub( paste("OFFSET=", params$OFFSET), paste("OFFSET=", OFFSETNEW), x )
    writeLines(y, procs)
}

align_all <- function(datadir, kind, left, right, where) {
    spectra_dirs <- get_spectra_dirs(datadir)

    for (dir in spectra_dirs) {
        align_spectrum(dir, kind, left, right, where)
    }
}

integrate_all <- function(datadir, kind, input) {
    # integrate spectra 
    require(caTools)

    SpectraDirs <- get_spectra_dirs(datadir)
    
    inte = matrix(nrow=length(SpectraDirs), ncol=nrow(input))
    colnames(inte) = input[,1]
    xx = list(list(),list()); yy=list(list(),list())
    
    i = 0
    for (dir in SpectraDirs) {
        procs_dir <- file.path(dir, as.character(kind), "pdata", "1")
      
        ## Read in parameters, needed for data processing
        params <- parse_spectrum_parameters(procs_dir)
  
        # normalize
        v <- which(params$x >= 5.5 | params$x <= 4.35)
        v <- sum(abs(params$y[v]))
        y <- 1e6*params$y/v
        i <- i+1
  
        for (ff in 1:nrow(input)) {
  
            left = input[ff,2]; right = input[ff,3]; where=input[ff,4]; 
            l1 = left + 0.0005;  l2 = left + 0.0001;
            r1 = right - 0.0001; r2 = right - 0.0005;
  
            OFFSETNEW <- params$OFFSET
            if (input[ff,6]==1){
  
                ## align
                v <- which(left >= params$x & right <= params$x)
                s <- which.max(y[v])
                reference <- params$x[v][s]
                shift <- where - reference
                OFFSETNEW <- params$OFFSET + shift
            }
  
            if (input[ff,5]==1){
                ### raddrizza
                x <- as.numeric(seq(OFFSETNEW, OFFSETNEW - params$SW_p/params$SF, length=params$SI))
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
            x <- as.numeric(seq (OFFSETNEW, OFFSETNEW - params$SW_p/params$SF, length=params$SI))
            v <- which(left >= x & right <= x)
  
            inte[i,ff] <- trapz(rev(x[v]), y[v])
  
            xx[[i]][[ff]] = rev(x[v])
            yy[[i]][[ff]] = rev(y[v])
            x <- as.numeric( seq(params$OFFSET, params$OFFSET - params$SW_p/params$SF, length=params$SI) )
        }
    }

    # report
    for (h in 1:nrow(input) ) {
        ss = paste(input[h,1], "ok") 
        pp = paste(input[h,1], "bad")

        reference = inte[1,h]
        test = inte[2,h]

        tt = test * 100 / reference
        ttt = abs(100 - tt)
        ifelse(ttt < 20, print(ss), print(pp))
        print(paste("diff", round(ttt,2), "%"))
    }

    return (list(inte=inte, x=xx, y=yy))
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
    kind <- 3
    left <- 5.257
    right <- 5.225
    where <- 5.243

    # percorso del file metabolites.txt
    metabolites_t <- read.table(metabolites)

    align_all(path_to_spectra, kind, left, right, where)
    res <- integrate_all(path_to_spectra, kind, metabolites_t)

    pdf("plot.pdf")
    plot_metabolites(metabolites_t, res$y, res$x)
    dev.off()
}

main(spectra)
