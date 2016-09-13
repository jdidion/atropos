library(reshpae2)
library(ggplot2)
library(readr)

compare_mapping <- function(read, prog) {
    untrimmed <- (read[1, 'read1_quality']+1) * (read[1, 'read2_quality']+1)
    trimmed <- (read[prog, 'read1_quality']+1) * (read[prog, 'read2_quality']+1)
    if (untrimmed == trimmed) {
        'same_quality'
    }
    else if (untrimmed > trimmed) {
        'worse_quality'
    }
    else {
        'better_quality'
    }
}
outcomes <- do.call(rbind, lapply(seq(1, nrow(tab), num.progs), function(i) {
    if (((i-1)/num.progs) %% 1000 == 0) print((i-1)/num.progs)
    read <- tab[i:(i+num.progs-1),]
    #read <- read[match(progs, read[,1]),]
    sapply(2:num.progs, function(prog) {
        if (read[1, 'skipped']) {
            'skipped'
        }
        else if (read[prog, 'discarded']) {
            'discarded'
        }
        else if (any(read[c(1,prog), 'split'])) {
            'split'
        }
        else {
            # untrimmed and trimmed have different validity
            if (read[1, 'valid'] != read[prog, 'valid']) {
                if (read[1, 'valid']) {
                    'valid->invalid'
                }
                else {
                    'invalid->valid'
                }
            }
            else if (!read[1, 'valid']) {
                'both_invalid'
            }
            # untrimmed is mapped
            else if (read[1, 'proper']) {
                if (read[prog, 'proper']) {
                    compare_mapping(read, prog)
                }
                else if (read[prog, 'read1_mapped'] || read[prog, 'read2_mapped']) {
                    'proper->partial'
                }
                else {
                    'proper->unmapped'
                }
            }
            else if (read[1, 'read1_mapped'] || read[1, 'read2_mapped']) {
                if (read[prog, 'proper']) {
                    'partial->proper'
                }
                else {
                    compare_mapping(read, prog)
                }
            }
            # untrimmed is unmapped
            else if (read[prog, 'proper']) {
                'unmapped->proper'
            }
            else if (read[prog, 'read1_mapped'] || read[prog, 'read2_mapped']) {
                'unmapped->partial'
            }
            else {
                'unmapped->unmapped'
            }
        }
    })
}))

process <- function(tab) {
    progs <- c('untrimmed', sort(setdiff(unique(tab$prog), 'untrimmed')))
    N <- max(tab$read_idx)
    num.progs <- length(progs)
    for (i in c(4:10,19)) {
        tab[,i] <- ifelse(tab[,i]=='True', TRUE, FALSE)
    }
    sample_tabs <- lapply(1:num.progs, function(i) tab[seq(i,nrow(tab),num.progs),])
    names(sample_tabs) <- unlist(lapply(sample_tabs, function(x) x[1,1]))
    quals <- rbind(
        data.frame(read_id=1:N, read=1, do.call(cbind, lapply(sample_tabs, function(x) x$read1_quality))),
        data.frame(read_id=1:N, read=2, do.call(cbind, lapply(sample_tabs, function(x) x$read2_quality)))
    )
    quals$maxq <- apply(quals[,3:ncol(quals)], 1, max)
    pts <- melt(as.data.frame(t(sapply(c(1,seq(5,60,5)), function(th) {
        c(
            th=th,
            x=sum(quals$maxq >= th),
            sapply(progs, function(p) {
                sum(quals[,p] >= th)
            })
        )
    }))), id.vars=c('th', 'x'), measure.vars=progs, variable.name='prog', value.name='y')
    list(progs=progs, quals=quals, pts=pts)
}

plot.data <- function(data) {
    ggplot(data$pts, aes(x=x, y=y, colour=prog, shape=prog)) + geom_line() + geom_point()
}

wgs <- read_tsv('wgs_results.txt')
wgs.data <- process(wgs)

wgbs <- read_tsv('wgbs_results.txt')
wgbs.data <- process(wgbs)

#sum(quals[,7] != -1 & quals[,3] != -1 & quals[,7] > quals[,3])
#sum(quals[,7] != -1 & quals[,3] != -1 & quals[,7] < quals[,3])
#sum(quals[,10] != -1 & quals[,3] != -1 & quals[,10] > quals[,3])
#sum(quals[,10] != -1 & quals[,3] != -1 & quals[,10] < quals[,3])

#w<-apply(quals[,c(3,10)], 1, function(x) any(x==-1))
#q<-quals[!w,]
#q[q<0] <- 0
#apply(q[,3:10],2,mean)
