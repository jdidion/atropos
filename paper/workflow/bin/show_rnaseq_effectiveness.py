#!/usr/bin/env python
import argparse
import pandas as pd
from compute_real_effectiveness import HEADER

def main():
    parser.add_argument("-i", "--input", default="-")
    parser.add_argument("-o", "--output")
    parser.add_argument(
        "-f", "--formats", choices=('txt', 'svg', 'pickle'), nargs="+",
        default=['tex'])
    parser.add_argument("--exclude-discarded", action="store_true", default=False)
    args = parser.parse_args()
    
    
    process.rna <- function(rna) {
        r1 <- rna[,c('prog','read1_in_region','read1_quality')]
        colnames(r1)[2:3] <- c('in_region', 'quality')
        r2 <- rna[,c('prog','read2_in_region','read2_quality')]
        colnames(r2)[2:3] <- c('in_region', 'quality')
        rna.reads <- rbind(r1,r2)
        rna.tab <- dlply(rna.reads, 'prog', table)
        rna.tab.tidy <- lapply(rna.tab[1:3], melt)
        rna.tab.tidy <- do.call(rbind, lapply(names(rna.tab.tidy), function(n) {
            x<-rna.tab.tidy[[n]]
            x$prog <- n
            x
        }))
        for (i in 1:3) rna.tab.tidy[[i]] <- rna.tab.tidy[[i]][2:3, 2:5]
        for (i in 1:3) rna.tab.tidy[[i]] <- (rna.tab.tidy[[i]][2,] - rna.tab.tidy[[4]][2,])
        rna.tab.tidy <- melt(do.call(rbind, rna.tab.tidy[1:3]))
        colnames(rna.tab.tidy) <- c('prog','MAPQ','delta')
        rna.tab.tidy$MAPQ <- factor(rna.tab.tidy$MAPQ)
        rna.tab.tidy
    }
    
    plot.rna.data <- function(rna.data) {
        ggplot(rna.data, aes(x=MAPQ, y=log10(delta), col=Program, group=Program)) +
            geom_line() +
            geom_point() +
            xlab('Mapping Quality Score (MAPQ) Cutoff') +
            ylab('Log10(Difference versus\nUntrimmed Reads)')
    }

    rna <- read_tsv('rnaseq_results.txt')
    rna.data <- process.rna(rna)
    pdf('Figure3.pdf')
    plot.rna.data(rna.data)
    dev.off()

if __name__ == "__main__":
    main()