plotOne <- function(qname, target, k=10000L){

    if(!qname %in% query$gff$seqid){
      warning(sprintf("'%s' is not present in the query GFF, cannot plot.", qname))
      return(NULL)
    }

    q <- query$gff[query$gff$seqid == qname]

    qids <- findOverlaps(
      q,
      target$syn$query,
      maxgap=k,
      ignore.strand=TRUE
    ) %$% subjectHits

    range2df <- function(r, side='query', group='syn'){
      df <- ranges(r) %>% as.data.frame
      df$chr <- r %>% seqnames %>% as.character
      df$names <- NULL # needed only draw from query$gff
      df$group <- group
      df$side <- side
      df$strand <- strand(r) %>% as.character
      df$id <- 1:nrow(df)
      df
    }

    bindAndPlot <- function(chrname, d) {
        td <- d[[chrname]]
        synids <- subset(td, group == "syn") %$% id

        td <- do.call(
            rbind,
            list(
                td,
                range2df(q, 'query', 'query'),
                range2df(qints, 'query', 'syn') %>% subset(id %in% synids)
            )
        )

        tmin <- subset(td, side=='target' & group=='syn') %$% start %>% min
        tmax <- subset(td, side=='target' & group=='syn') %$% end %>% max
        qmin <- subset(td, side=='query' & group=='syn') %$% start %>% min
        qmax <- subset(td, side=='query' & group=='syn') %$% end %>% max


        tout <- target$syn$target %>% split(target$syn$target %>% seqnames)
        tout <- tout[[chrname]]
        tout <- tout[start(tout) >= tmin & end(tout) <= tmax]

        qchr <- td[td$side=='query', ]$chr[1]
        qout <- target$syn$query %>% split(target$syn$query %>% seqnames)
        qout <- qout[[qchr]]
        qout <- qout[start(qout) >= qmin & end(qout) <= qmax]

        td <- do.call(
            rbind,
            list(
                td,
                range2df(tout, 'target', 'out'),
                range2df(qout, 'query', 'out')
            )
        )

        td <- group_by(td, side, chr) %>%
            mutate(start=start - min(start) + 1, end=start + width - 1) %>%
            ungroup

        td$y <- 0
        td$y[td$side == 'query'  & td$group == 'syn'] <- -1
        td$y[td$side == 'target' & td$group == 'syn'] <-  1 
        td$y[td$side == 'target' & td$group == 'si'] <-  1.3
        td$y <- td$y + seq(from=-0.10, to=0.10, length.out=5)[(1:nrow(td)) %% 5 + 1]
        td$y[td$side == 'query'  & td$group == 'query'] <- -1.2
        if(sum(td$group == 'si') < 2){
          td$y[td$side == 'target' & td$group == 'si'] <-  1.2
        }
        td$y[td$side == 'query'  & td$group == 'out'] <- -1.2

        if(sum(td$group == 'si') > 0){
          td$y[td$side == 'target' & td$group == 'out'] <- max(td[td$group == 'si', 'y']) + 0.1
        } else {
          td$y[td$side == 'target' & td$group == 'out'] <- 1.2
        }

        intlines <- data.frame(
            qx=subset(td, side=='query' & group=='syn') %>% with((start + end) / 2),
            qy=subset(td, side=='query' & group=='syn') %$% y,
            tx=subset(td, side=='target' & group=='syn') %>% with((start + end) / 2),
            ty=subset(td, side=='target' & group=='syn') %$% y
        )

        tgenome=attr(target$si$target, "seqinfo")@genome[1]
        qgenome=attr(target$si$query, "seqinfo")@genome[1]

        ggplot() +
            ggtitle(qname) +
            xlab(
              sprintf(
                "top  %s (%s : %d)\nbottom  %s (%s : %d)",
                tgenome, chrname, tmin, qgenome, qchr, qmin
              )
            ) +
            geom_segment(
                data=subset(td, group=='query'),
                aes(
                    x=start,
                    xend=end,
                    y=y,
                    yend=y
                ),
                color='blue',
                size=4,
            ) +
            geom_segment(
                data=subset(td, group=='out'),
                aes(
                    x=start,
                    xend=end,
                    y=y,
                    yend=y
                ),
                size=2,
                alpha=0.5,
                color='black'
            ) +
            geom_segment(
                data=subset(td, group=='syn'),
                aes(
                    x=start,
                    xend=end,
                    y=y,
                    yend=y,
                    color=strand
                ),
                alpha=0.5,
                size=2,
            ) +
            geom_segment(
                data=subset(td, group=='si'),
                aes(
                    x=start,
                    xend=end,
                    y=y,
                    yend=y
                ),
                size=3,
                alpha=0.5,
                color='orange'
            ) +
            geom_segment(
                data=intlines,
                aes(
                    x=qx,
                    xend=tx,
                    y=qy,
                    yend=ty
                ),
                alpha=0.2
            )
    }

    qints <- target$syn$query[qids]
    tints <- target$syn$target[qints$over]
    qsi   <- target$si$query[which(target$si$query$seqid == qname)]
    tsi   <- target$si$target[qsi$id]

    d <- do.call(
        rbind,
        list(
            range2df(tints, 'target', 'syn'),
            range2df(tsi, 'target', 'si')
        )
    )
    d <- split(d, d$chr)
    gplots <- lapply(names(d), bindAndPlot, d)

    do.call(grid.arrange, c(gplots, ncol=1))
}
