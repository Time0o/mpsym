library(colorspace)
library(ggplot2)
library(gridExtra)
library(scales)
library(tikzDevice)

plotFile_ <- function(plot_dir, plot_name, plot_type) {
    paste(file.path(plot_dir, tolower(plot_name)), '_', plot_type, '.tex', sep='')
}

plot_ <- function(d, z, titles, f, plot_file, width=6, height=4) {
    tikz(file=plot_file, width=width, height=height)

    if (is.vector(d)) {
        plots = list()
        offs = 0
        for (i in 1:length(d)) {
            plots[[i]] = f(d[[i]], titles[[i]], offs)
            offs = offs + length(unique(d[[i]][[z]]))
        }

        grid.arrange(grobs=plots, nrow=1, ncol=length(plots))

    } else {
        plot = f(d, titles, 1)

        print(plot)
    }

    invisible(dev.off())
}

default_colors_ <- function() {
  hues = seq(15, 375, length = 8)
  hcl(h = hues, l = 65, c = 100)
}

createBoxplot <- function(d, x, y, ylabel, z, titles, plot_dir, plot_name) {
    plot_(
        d,
        z,
        titles,
        function (d_, title_, offs) {
            ggplot(na.omit(d_), aes_string(x=x, y=y, fill=z)) +
                stat_boxplot(geom='errorbar') +
                geom_boxplot(outlier.size=1) +
                scale_y_continuous(trans='log10', labels=scientific) +
                theme_bw() +
                theme(axis.title.x=element_blank(),
                      axis.title.y=element_text(margin=margin(r=5)),
                      plot.title=element_text(hjust=0.5),
                      legend.title=element_blank(),
                      legend.position='top',
                      legend.direction='vertical',
                      legend.margin=margin(t=-0.25,
                                           unit='cm')) +
                labs(y=ylabel, title=title_)
        },
        plotFile_(plot_dir, plot_name, 'boxplot')
    )
}

createLinePlot <- function(
    d, x, xlabel, y, ylabel, z, titles, annotate, plot_dir, plot_name) {

    plot_(
        d,
        z,
        titles,
        function(d_, title_, offs) {
            # x limits
            xmin = min(d_[[x]])
            xmax = max(d_[[x]])

            xrange = xmax - xmin
            xmin = xmin - xrange / 10
            xmax = xmax + xrange / 10

            # y limits
            ymin = min(d_[[y]])
            ymax = max(d_[[y]])

            logy = abs(ymax / ymin) > 100

            if (logy)
                ymin = 10^floor(log10(min(d_[[y]])))
            else
                ymin = min(d_[[y]]) - (max(d_[[y]]) - min(d_[[y]])) / 5

            nlines = length(unique(d_[[z]]))

            # colors
            if (offs == 0)
                colors = default_colors_()[1:nlines]
            else
                colors = default_colors_()[(offs + 1):(nlines + offs)]

            # legend dimensions
            lcols = if (nlines %% 2 ==  0) 2 else 3
            lrows = floor(nlines / lcols)

            plot = ggplot(na.omit(d_),
                          aes_string(x=x, y=y, color=z, xmin=xmin, xmax=xmax, ymin=ymin)) +
                geom_line() +
                geom_point(aes_string(shape=z)) +
                scale_x_continuous(breaks=unique(d_[[x]])) +
                scale_y_continuous(trans=(if (logy) 'log10' else 'identity'),
                                   labels=scientific) +
                theme_bw() +
                theme(axis.title.y=element_text(margin=margin(r=5)),
                      plot.title=element_text(hjust=0.5),
                      legend.title=element_blank(),
                      legend.position='top',
                      legend.direction='vertical',
                      legend.margin=margin(t=-0.25,
                                           b=max((3 - lrows) * 0.5, 0),
                                           unit='cm')) +
                guides(color=guide_legend(ncol=lcols,byrow=TRUE)) +
                labs(x=xlabel, y=ylabel, title=title_) +
                scale_color_manual(values=colors)

            plot = annotate(plot, d_)
        },
        plotFile_(plot_dir, plot_name, 'lineplot')
    )
}
