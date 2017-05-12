#!/usr/bin/env python

def main():
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