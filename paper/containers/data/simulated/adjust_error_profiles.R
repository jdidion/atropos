# This script adjusts an ART error profile to increase the cumulative error
# rate to a given value.

## Optimization of multiplier of counts vector to generate
## a distribution with (approximately) a given error rate
q2p <- function(q) 10^(-q/10)
cum_err <- function(p, counts) {
    s <- sum(counts)
    sum(p * counts) / s
}
eq <- function(a, b, tol) {
    return (abs(a-b) <= tol)
}
new_err <- function(d, p, counts, target.err) {
    total <- sum(counts)
    new.counts <- counts + (total * d)
    abs(cum_err(p, (new.counts / sum(new.counts)) * total) - target.err)
}
find_new_counts <- function(p, counts, target.err, rescale=TRUE, d0=target.err*0.1) {
    res <- optim(par=d0, fn=new_err, p=p, counts=counts, target.err=target.err, method = "BFGS")
    new.counts <- counts + (sum(counts) * res$par)
    if (rescale) {
        m <- sum(counts) / sum(new.counts)
        new.counts <- new.counts * m
    }
    list(counts=new.counts, res=res)
}

profile.error <- function(path) {
    if (endsWith(path, ".gz")) {
        lines <- readLines(gzfile(path))
    }
    else {
        lines <- readLines(path)
    }
    comments <- grep("^#", lines, perl=T)
    if (length(comments) > 0) {
        lines <- lines[-comments]
    }
    num <- NULL
    denom <- NULL
    for (i in seq(1, length(lines), 2)) {
        q <- unlist(strsplit(lines[i], "\t"))
        if (length(q) < 3) {
            next
        }
        h <- q[1:2]
        q <- as.integer(q[3:length(q)])
        p <- sapply(q, q2p)
        counts <- unlist(strsplit(lines[i+1], "\t"))
        stopifnot(counts[1] == h[1] && counts[2] == h[2])
        counts <- as.numeric(counts[3:length(counts)])
        stopifnot(length(counts) == length(p))
        num <- c(num, p*counts)
        denom <- c(denom, counts)
    }
    as.numeric(sum(num)) / as.numeric(sum(denom))
}

## Generate a new error profile
adjust.profile <- function(orig.path, new.path, target.err, rescale=TRUE) {
    if (endsWith(orig.path, ".gz")) {
        lines <- readLines(gzfile(orig.path))
    }
    else {
        lines <- readLines(orig.path)
    }
    comments <- grep("^#", lines, perl=T)
    if (length(comments) > 0) {
        lines <- lines[-comments]
    }
    n.lines <- length(lines)
    read.length <- n.lines / 12
    new.lines <- lines
    process.lines <- function(i) {
        q <- unlist(strsplit(lines[i], "\t"))
        h <- q[1:2]
        q <- as.integer(q[3:length(q)])
        p <- sapply(q, q2p)
        counts <- unlist(strsplit(lines[i+1], "\t"))
        stopifnot(counts[1] == h[1] && counts[2] == h[2])
        counts <- as.numeric(counts[3:length(counts)])
        stopifnot(length(counts) == length(p))
        new.counts <- find_new_counts(p, counts, target.err, rescale)
        list(h=h, q=q, p=p, counts=new.counts$counts, line=paste0(c(h, round(new.counts$counts)), collapse="\t"))
    }
    new.N <- rep(NA, read.length)
    for (ii in 1:read.length) {
        i <- (ii*2)-1
        result <- process.lines(i)
        new.lines[i+1] <- result$line
        w <- which(result$q == 2)
        if (length(w) > 0) {
            new.N[ii] <- result$counts[w[1]]
        }
    }
    end <- read.length * 2 * 5
    for (i in seq((read.length * 2) + 1, end, 2)) {
        result <- process.lines(i)
        new.lines[i+1] <- result$line
    }
    for (i in 1:read.length) {
        if (!is.na(new.N[i])) {
            new.lines[end+(i*2)] <- paste0(c('N', i-1, new.N[i]), collapse="\t")
        }
    }
    if (endsWith(new.path, ".gz")) {
        out <- gzfile(new.path, "w")
    }
    else {
        out <- file(new.path, "w")
    }
    writeLines(new.lines, out)
    close(out)
}

parser <- ArgumentParser()
parser$add_argument(
    "-g", "--generate", action="store_true", default=FALSE,
    help="Generate new profile.")
parser$add_argument(
    "-c", "--calculate", action="store_true", default=FALSE,
    help="Calculate error rate on existing profile.")
parser$add_argument(
    "-p", "--profile", help="The source ART error profile.")
parser$add_argument(
    "-e", "--error-rate", type="double", action="store_true", default=FALSE,
    help="Calculate error rate of existing profile.")
parser$add_argument(
    "-o", "--output", default=NULL, 
    help="New profile file.")
parser$add_argument(
    "-t", "--table", default=NULL, 
    help="File where data will be written in tab-delimited format.")
args <- parser$parse_args()

if (args$generate) {
    # for some reason, the acutal error rate in simulated reads is about twice 
    # what is expected based on the empirical distribution of quality scores
    error.rate <- as.numeric(args$error_rate) / 2
    adjust.profile(args$profile, args$output, error.rate)
    actual.error <- profile.error(args$output)
    print(paste("The acutal error rate of the new profile is", actual.error))
} else if (args$error) {
    print(paste("The acutal error rate of the profile is", round(profile.error(args$profile), 4)))
}

if (args$table) {
    if (!file.exists(args$table)) {
        write(
            paste("Profile", "RequestedErrorRate", "ActualErrorRate", sep="\t"),
            file=args$table)
    }
    write(
        paste(basename(args$profile), error.rate, actual.error, sep="\t"),
        file=args$table, append=TRUE)
}