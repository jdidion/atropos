#!/usr/bin/env python

def main():
    process.wgbs <- function(tab, exclude.discarded=TRUE, has.regions=FALSE) {
        progs <- c('untrimmed', sort(setdiff(unique(tab$prog), 'untrimmed')))
        N <- max(tab$read_idx)
        num.progs <- length(progs)
        for (i in c(4:10, ifelse(has.regions, 20, 19))) {
            tab[,i] <- ifelse(tab[,i]=='True', TRUE, FALSE)
        }
        sample_tabs <- lapply(1:num.progs, function(i) tab[seq(i,nrow(tab),num.progs),])
        names(sample_tabs) <- unlist(lapply(sample_tabs, function(x) x[1,1]))
        quals <- rbind(
            data.frame(read_id=1:N, read=1, do.call(cbind, lapply(sample_tabs, function(x) x$read1_quality))),
            data.frame(read_id=1:N, read=2, do.call(cbind, lapply(sample_tabs, function(x) x$read2_quality)))
        )
        if (has.regions) {
            in_region <- rbind(
                data.frame(read_id=1:N, read=1, do.call(cbind, lapply(sample_tabs, function(x) x$read1_in_region))),
                data.frame(read_id=1:N, read=2, do.call(cbind, lapply(sample_tabs, function(x) x$read2_in_region)))
            )
        }
        if (exclude.discarded) {
            w <- apply(quals[,3:ncol(quals)], 1, function(x) any(x==-1))
            if (has.regions) {
                w <- w | apply(in_region[,3:ncol(in_region)], 1, function(x) any(x==-1))
                in_region <- in_region[!w,]
            }
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
        if (has.regions) {
            in_region$maxi <- apply(in_region[,3:ncol(in_region)], 1, max)
            pts2 <- melt(as.data.frame(t(sapply(th_values, function(th) {
                c(
                    th=th,
                    x=sum(in_region[quals$maxq >= th, 'max']==1),
                    sapply(progs, function(p) {
                        sum(in_region[quals[,p] >= th, p]==1)
                    })
                )
            }))), id.vars=c('th', 'x'), measure.vars=progs, variable.name='prog', value.name='y')
            pts2$delta <- NA
            for (th in th_values) {
                for (prog in progs) {
                    pts2[pts2$prog == prog & pts2$th == th, 'delta'] <- pts2[pts2$prog == prog & pts2$th == th, 'y'] - pts2[pts2$prog == 'untrimmed' & pts2$th == th, 'y']
                }
            }
            retval <- c(retval, list(reg=in_region, pts2=pts2))
        }
        retval
    }

if __name__ == "__main__":
    main()