#' Create subheading for each celltype and its corresponding plot
#'
#' For use in RMarkdown code chunks with option `results="asis"`.
#'
#' @param plotlist List of plots. Can be nested.
#' @param use_names Names to be used in the subheading. If not provided, will
#'   use `names(plotlist)`.
#' @param header_level Level of the subheading to use. Default: 3.
#'
#' @export
print_plots <- function(plotlist, use_names = NULL, header_level = 3) {
    stopifnot("Expects list as input." = is.list(plotlist))

    if (!is.null(use_names)) {
        names(plotlist) <- use_names
    }

    ## Helper to recursively traverse through list and print out gg objects
    print_or_walk <- function(x) {
        if (!inherits(x, "gg")) {
            purrr::walk(x, print_or_walk)
        } else {
            print(x)
        }
    }

    purrr::iwalk(plotlist, function(x, i) {
        cat("\n\n", paste(rep("#", header_level), collapse = ""), " ",
            i, "\n\n", sep = ""
        )
        print_or_walk(x)
    })
}


## Boxplot of runtimes per method
#' @export
plot_run_times <- function(x, ...) {
    ggplot(x, aes(method, time, col = method)) +
        geom_jitter(...) +
        labs(x = "") +
        coord_flip()
}

## Plot p-value histogram with ggplot, faceted by method replicate
#' @export
pval_hist <- function(x, binwidth = 0.05, ...,
                      stratified = FALSE, n_groups = 5) {
    methods <- unique(x$method)
    labs <- stringr::str_replace(methods, "nonClustered", "nonClustered\n")
    names(labs) <- methods

    p <- ggplot(x, aes(PValue, fill = method)) +
        geom_histogram(
            col = "black", binwidth = binwidth, boundary = 0,
            show.legend = FALSE, ...
        ) +
        scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
            theme_minimal() +
            theme(
                strip.text = element_text(size = 10, face = "bold"),
                strip.text.y = element_text(angle = 0, hjust = 0)
            )

    if (stratified) {
        p +
            facet_grid(
                rows = vars(method),
                cols = vars(cut_number(avg_count, n = n_groups)),
                scales = "free",
                labeller = labeller(method = labs)
            )
    } else {
        p +
            facet_grid(
                rows = vars(method), cols = vars(replicate), scales = "free",
                labeller = labeller(method = labs)
            )
    }
}

#' @export
pval_hist_strat <- function(..., n_groups = 5) {
    pval_hist(..., stratified = TRUE, n_groups = n_groups)
}


#' @export
#' @importFrom purrr map
#' @importFrom iCOBRA fdrtpr plot_fdrtprcurve
#' @importFrom ggplot2 coord_cartesian
plot_zoomed_fdrtprcurve <- function(cobraplot, ...) {
    p <- fdrtpr(cobraplot)
    x_lims <- c(min(p$FDR), max(p$FDR))
    y_lims <- c(min(p$TPR), max(p$TPR))
    plot_fdrtprcurve(cobraplot, ...) +
        coord_cartesian(xlim = x_lims, ylim = y_lims)
}



# iCOBRA plots ------------------------------------------------------------

#' @export
plot_fdrtpr_points <- function(x, thr = c(0.01, 0.05, 0.1)) {
    ggplot(x, aes(FDR, TPR, color = method, fill = method)) +
        geom_vline(xintercept = thr, linetype = "dashed", alpha = 0.5) +
        geom_path(size = 0.8, alpha = 0.6, show.legend = FALSE) +
        geom_point(size = 4, stroke = 0.8, alpha = 0.8, shape = 21, col = "black")
}


#' @export
plot_fdr_control <- function(cobra_data, use_thr, shape_by = NULL) {
    dat_text <- data.frame(
        label = paste("alpha == ", use_thr),
        method = levels(cobra_data$method)[1],
        FDR = use_thr, TPR = 0
    )

    base_jitter <- function(...) {
        geom_jitter(
            width = 0.2, height = 0,
            alpha = 0.8, stroke = 0.8,
            ...
        )
    }
    if (!is.null(shape_by)) {
        shape_by <- sym(shape_by)
        add_jitter <- base_jitter(mapping = aes(shape = !!enquo(shape_by)))
    } else {
        add_jitter <- base_jitter(shape = 21)
    }

    ggplot(cobra_data, aes(method, FDR, fill = method, size = TPR)) +
        add_jitter +
        geom_hline(yintercept = use_thr, linetype = "dashed", alpha = 0.7) +
        geom_text(
            data = dat_text, aes(label = label), size = 6,
            alpha = 0.5, hjust = 1, vjust = -0.5, parse = TRUE
        ) +
        labs(x = NULL, y = "FDP")
}
