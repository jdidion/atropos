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
    
    with fileopen(args.input, 'rt') as inp:
        table = pd.read_csv(inp, sep="\t", names=HEADER)
    
    if 'txt' in args.formats:
        with fileopen(args.output + '.txt', 'wt') as out:
            table.to_csv(out, sep="\t", index=False)
        
    if 'pickle' in args.formats:
        import pickle
        with fileopen(args.output + '.pickle', 'wb') as out:
            pickle.dump(table, out)
    
    if 'svg' in args.formats:
        # sort programs but put untrimmed first
        progs = set(tab.prog)
        progs.remove('untrimmed')
        ['untrimmed'] + sorted(progs)
        
        num_reads = max(tab.read_idx)
        
        quals = table.pivot(index='read_idx', columns='prog', values='read1_quality').\
            append(table.pivot(index='read_idx', columns='prog', values='read2_quality')).\
            drop('')
        
        if args.exclude_discarded:
            tab3.apply(lambda x: any(x==-1), 1)
    
        if (exclude.discarded) {
            w <- apply(quals[,3:ncol(quals)], 1, function(x) any(x==-1))
            quals <- quals[!w,]
        }
        quals$maxq <- apply(quals[,3:ncol(quals)], 1, max)
        th_values <- c(1,seq(5,60,5))
        pts <- melt(as.data.frame(t(sapply(th_values, function(th) {
            c(
                th=th,
                x=sum(quals$maxq >= th),
                sapply(progs, function(p) {
                    sum(quals[,p] >= th)
                })
            )
        }))), id.vars=c('th', 'x'), measure.vars=progs, variable.name='prog', value.name='y')
        pts$delta <- NA
        for (th in th_values) {
            for (prog in progs) {
                pts[pts$prog == prog & pts$th == th, 'delta'] <- pts[pts$prog == prog & pts$th == th, 'y'] - pts[pts$prog == 'untrimmed' & pts$th == th, 'y']
            }
        }
        retval<-list(progs=progs, quals=quals, pts=pts)
        retval
    }

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

if __name__ == "__main__":
    main()


# process.wgbs <- function(tab, ) {
#     progs <- c('untrimmed', sort(setdiff(unique(tab$prog), 'untrimmed')))
#     N <- max(tab$read_idx)
#     num.progs <- length(progs)
#     for (i in c(4:10, ifelse(has.regions, 20, 19))) {
#         tab[,i] <- ifelse(tab[,i]=='True', TRUE, FALSE)
#     }
#     sample_tabs <- lapply(1:num.progs, function(i) tab[seq(i,nrow(tab),num.progs),])
#     names(sample_tabs) <- unlist(lapply(sample_tabs, function(x) x[1,1]))
#     quals <- rbind(
#         data.frame(read_id=1:N, read=1, do.call(cbind, lapply(sample_tabs, function(x) x$read1_quality))),
#         data.frame(read_id=1:N, read=2, do.call(cbind, lapply(sample_tabs, function(x) x$read2_quality)))
#     )
#     if (has.regions) {
#         in_region <- rbind(
#             data.frame(read_id=1:N, read=1, do.call(cbind, lapply(sample_tabs, function(x) x$read1_in_region))),
#             data.frame(read_id=1:N, read=2, do.call(cbind, lapply(sample_tabs, function(x) x$read2_in_region)))
#         )
#     }
#     if (exclude.discarded) {
#         w <- apply(quals[,3:ncol(quals)], 1, function(x) any(x==-1))
#         if (has.regions) {
#             w <- w | apply(in_region[,3:ncol(in_region)], 1, function(x) any(x==-1))
#             in_region <- in_region[!w,]
#         }
#         quals <- quals[!w,]
#     }
#     quals$maxq <- apply(quals[,3:ncol(quals)], 1, max)
#     th_values <- c(1,seq(5,60,5))
#     pts <- melt(as.data.frame(t(sapply(th_values, function(th) {
#         c(
#             th=th,
#             x=sum(quals$maxq >= th),
#             sapply(progs, function(p) {
#                 sum(quals[,p] >= th)
#             })
#         )
#     }))), id.vars=c('th', 'x'), measure.vars=progs, variable.name='prog', value.name='y')
#     pts$delta <- NA
#     for (th in th_values) {
#         for (prog in progs) {
#             pts[pts$prog == prog & pts$th == th, 'delta'] <- pts[pts$prog == prog & pts$th == th, 'y'] - pts[pts$prog == 'untrimmed' & pts$th == th, 'y']
#         }
#     }
#     retval<-list(progs=progs, quals=quals, pts=pts)
#     if (has.regions) {
#         in_region$maxi <- apply(in_region[,3:ncol(in_region)], 1, max)
#         pts2 <- melt(as.data.frame(t(sapply(th_values, function(th) {
#             c(
#                 th=th,
#                 x=sum(in_region[quals$maxq >= th, 'max']==1),
#                 sapply(progs, function(p) {
#                     sum(in_region[quals[,p] >= th, p]==1)
#                 })
#             )
#         }))), id.vars=c('th', 'x'), measure.vars=progs, variable.name='prog', value.name='y')
#         pts2$delta <- NA
#         for (th in th_values) {
#             for (prog in progs) {
#                 pts2[pts2$prog == prog & pts2$th == th, 'delta'] <- pts2[pts2$prog == prog & pts2$th == th, 'y'] - pts2[pts2$prog == 'untrimmed' & pts2$th == th, 'y']
#             }
#         }
#         retval <- c(retval, list(reg=in_region, pts2=pts2))
#     }
#     retval
# }

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