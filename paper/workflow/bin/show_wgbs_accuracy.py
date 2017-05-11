#!/usr/bin/env python

def main():
    plot.wgbs.data <- function(data, q.labels=NULL) {
        pts <- data$pts
        progs <- data$progs
        pts$prog <- as.character(pts$prog)
        colnames(pts) <- c("MAPQ", "x", "Q", "y", "Delta")
        if (!is.null(q.labels)) {
            for (i in 1:length(q.labels)) {
                pts[pts$Q == progs[i], 'Q'] <- names(q.labels)[i]
            }
            n <- names(sort(q.labels))
        }
        else {
            n <- sort(unique(pts$Q))
        }
        pts$Q <- factor(pts$Q, levels=n)
        pts <- pts[order(pts$Q),]
        ggplot(pts[pts$Q != 'untrimmed',], aes(x=MAPQ, y=Delta, colour=Q)) +
            geom_line() + geom_point() +
            labs(x="Mapping Quality Score (MAPQ) Cutoff", y="Difference versus Untrimmed Reads")
    }

    wgbs <- read_tsv('wgbs_results.txt')
    wgbs.data <- process.wgbs(wgbs, exclude.discarded=TRUE)
    pdf('Figure2.pdf')
    plot.wgbs.data(wgbs.data)
    dev.off()

if __name__ == "__main__":
    main()


# plot.data <- function(data, prog.labels) {
#     pts <- data$pts
#     progs <- data$progs
#     pts$prog <- as.character(pts$prog)
#     colnames(pts) <- c("MAPQ", "x", "Program", "y", "Delta")
#     for (i in 1:length(prog.labels)) {
#         pts[pts$Program == progs[i], 'Program'] <- prog.labels[i]
#     }
#     pts$Q <- 0
#     pts[grep(pts$Program, pattern = 'Q20'), 'Q'] <- 20
#     pts$Q <- factor(pts$Q)
#     ggplot(pts[pts$Program != 'Untrimmed',], aes(x=MAPQ, y=Delta, colour=Program, shape=Q)) +
#         geom_line() + geom_point() +
#         labs(x="Mapping Quality Score (MAPQ) Cutoff", y="Difference versus Untrimmed Reads")
# }