#    This file is part of GeneFunnel-benchmarks.
#    Copyright (C) 2021  Emir Turkes, UK DRI at UCL
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Emir Turkes can be contacted at emir.turkes@eturkes.com

# This file holds common functions and methods.

#' ggplot2 function providing custom aesthetics and automatic placement of
#' categorical labels. For continuous data, a colorbar is implemented.
#'
#' @param data SingleCellExperiment or Seurat object.
#' @param x,y Dimensionality reduction coordinates.
#' @param color Column metadata to color points by.
#' @param type \code{"cat"} is categorical, \code{"cont"} is continuous,
#' \code{"NULL"} is generic.
#' @examples
#' red_dim_plot(
#'     data = sce, x = "tsne1", y = "tsne2", color = "cluster", type = "cat"
#' )
#' red_dim_plot(
#'     data = seurat, x = "umap1", y = "umap2", color = "nUMI", type = "cont"
#')
#'
red_dim_plot <- function(data, x, y, color, type = NULL) {

    if ((class(data))[1] == "SingleCellExperiment") {
        gg_df <- data.frame(colData(data)[ , c(x, y, color)])
    } else if ((class(data))[1] == "Seurat") {
        gg_df <- data.frame(data[[x]], data[[y]], data[[color]])
    }
    rownames(gg_df) <- NULL
    gg_df[[color]] <- factor(gg_df[[color]])

    gg <- ggplot(gg_df, aes_string(x, y, col = color)) +
        geom_point(
            alpha = 0.35, stroke = 0.05, shape = 21, aes_string(fill = color)
        ) +
        theme_classic() +
        theme(
            legend.position = "right", plot.title = element_text(hjust = 0.5),
            legend.title = element_blank()
        ) +
        guides(color = guide_legend(override.aes = list(alpha = 1)))

    if (is.null(type)) {
        return(gg)

    } else if (type == "cat") {
        label_df <- gg_df %>%
            group_by_at(color) %>% summarise_at(vars(x:y), median)
        label_df <- cbind(label_df[[1]], label_df)
        names(label_df) <- c("label", color, x, y)
        gg <- gg +
            geom_label_repel(
                data = label_df, aes(label = label), show.legend = FALSE
            )

    } else if (type == "cont") {
        if ((class(data))[1] == "SingleCellExperiment") {
            gg_df <- data.frame(colData(data)[ , c(x, y, color)])
        } else if ((class(data))[1] == "Seurat") {
            gg_df <- data.frame(data[[x]], data[[y]], data[[color]])
        }
        rownames(gg_df) <- NULL

        gg <- ggplot(gg_df, aes_string(x, y)) +
            geom_point(alpha = 0.35, stroke = 0.05, aes_string(color = color)) +
            theme_classic() +
            theme(
                legend.position = "right",
                plot.title = element_text(hjust = 0.5),
                legend.title = element_blank()
            ) +
            scale_color_viridis()
    }
    gg
}

#' Adds download buttons and horizontal scrolling to \code{"DT::datatable"}.
#'
#' @param dt A data.table object.
#' @examples
#' datatable_download(dt = data_table)
#'
datatable_download <- function(dt) {

    datatable(
        dt,
        list(
            scrollX = TRUE,
            dom = "Bfrtip",
            buttons = list(
                "copy", "print",
                list(
                    extend = "collection", buttons = c("csv", "excel", "pdf"),
                    text = "Download"
                )
            ),
            pageLength = 6
        ),
        extensions = "Buttons"
    )
}

#' Adds download buttons, horizontal scrolling, exponential values to
#' \code{"DT::datatable"}.
#'
#' @param dt A data.table object.
#' @examples
#' datatable_download_exp(dt = data_table)
#'
datatable_download_exp <- function(dt) {

    datatable(
        dt,
        list(
            scrollX = TRUE,
            dom = "Bfrtip",
            buttons = list(
                "copy", "print",
                list(
                    extend = "collection",
                    buttons = c("csv", "excel", "pdf"), text = "Download")
            ),
            rowCallback = JS(
                "function(row, data) {",
                "for (i = 1; i < data.length; i++) {",
                "if (data[i]>=1000 | data[i]<1000) {",
                "$('td:eq('+i+')', row).html(data[i].toExponential(2));}}}"
            ),
            pageLength = 8
        ),
        extensions = "Buttons",
    )
}
