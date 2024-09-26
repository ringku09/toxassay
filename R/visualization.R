#' #----------------- setup
#' library(grid)
#' library(gridExtra)
#' library(ggthemes)
#' #library(extrafont)
#' library(latex2exp)
#' library(scales)
#' # font_import()
#' # loadfonts(device = "win", quiet = TRUE)
#' # fonts()
#'
#' common_legend<-function(gplot){
#'   tmp <- ggplot_gtable(ggplot_build(gplot))
#'   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#'   legend <- tmp$grobs[[leg]]
#'   return(legend)}
#'
#'
#' theme_Publication <- function(base_size=14,
#'                               base_family="Helvetica",
#'                               title_size = rel(1),
#'                               axis_size = rel(1),
#'                               lab_size = rel(0.8),
#'                               leg_size = rel(0.8),
#'                               strip_size = rel(0.8),
#'                               y_lab = element_text(angle=90,vjust =2),
#'                               x_lab = element_text(vjust = -0.2),
#'                               tag_size = rel(0.7),
#'                               top_mar = 5,
#'                               bot_mar = 5,
#'                               ...) {
#'   (theme_foundation(base_size=base_size, base_family=base_family)
#'    + theme(plot.title = element_text(face = "bold",
#'                                      size = rel(title_size), hjust = 0.5),
#'            text = element_text(),
#'            panel.background = element_rect(colour = NA),
#'            plot.background = element_rect(colour = NA),
#'            plot.tag = element_text(face = "bold", size = tag_size),
#'            panel.border = element_rect(colour = NA),
#'            axis.title = element_text(face = "bold",size = axis_size),
#'            axis.title.y = y_lab,
#'            axis.title.x = x_lab,
#'            axis.text = element_text(size = lab_size),
#'            axis.line = element_line(colour="black"),
#'            axis.ticks = element_line(),
#'            panel.grid.major = element_line(colour="#f0f0f0"),
#'            panel.grid.minor = element_blank(),
#'            legend.key = element_rect(colour = NA),
#'            legend.position = "bottom",
#'            legend.direction = "horizontal",
#'            legend.text = element_text(size= leg_size),
#'            legend.key.size= unit(0.4, "cm"),
#'            legend.margin = margin(0,unit = "cm"),
#'            legend.title = element_text(face="italic",size = leg_size),
#'            plot.margin=unit(c(top_mar,3,bot_mar,3),"mm"),
#'            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#'            strip.text = element_text(face="bold", size= strip_size)
#'    ))
#'
#' }
#'
#' scale_fill_Publication <- function(...){
#'   library(scales)
#'   discrete_scale("fill",
#'                  "Publication",
#'                  manual_pal(values = c("#1F78B4", "#33A02C", "#E31A1C","#6A3D9A", "#FF00FF",
#'                                        "#8B4513", "#556B2F", "yellow3","#FF9F00", "paleturquoise4")),
#'                  ...)
#'
#' }
#'
#' scale_colour_Publication <- function(...){
#'   library(scales)
#'   discrete_scale("colour",
#'                  "Publication",
#'                  manual_pal(values = c("#fdb462",
#'                                        "#386cb0",
#'                                        "#ef3b2c",
#'                                        "#7fc97f",
#'                                        "#662506",
#'                                        "#a6cee3",
#'                                        "#fb9a99",
#'                                        "#984ea3",
#'                                        "#ffff33")),
#'                  ...)
#'
#' }
#'
#' scale_colour_Publication2 <- function(...){
#'   library(scales)
#'   discrete_scale("colour",
#'                  "Publication",
#'                  manual_pal(values = c("#7fc97f",
#'                                        "#662506",
#'                                        "#a6cee3",
#'                                        "#fb9a99",
#'                                        "#ffff33",
#'                                        "#984ea3",
#'                                        "#ef3b2c",
#'                                        "#386cb0",
#'                                        "#fdb462")),
#'                  ...)
#'
#' }
#' #--------------------------end setup
#'
#'
#'
#' #             ppi plot
#' #' plot default string network
#' #'
#' #' @inheritParams get_netdata
#' #'
#' #' @return A ppi network plot
#' #' @export
#' #'
#' #' @examples
#' #' \donttest{
#' #' plot_stringnet(res3[1:50, ],
#' #'             gene_col = "gene_name",
#' #'             organism = "rat",
#' #'             score_threshold = 200)
#' #' }
#' plot_stringnet <- function(gene_df,
#'                            gene_col = "gene_name",
#'                            organism = "rat",
#'                            score_threshold = 200) {
#'   string_db <- setup_stringdb(organism = organism, score_threshold = score_threshold)
#'   gene_map <- string_db$map(as.data.frame(gene_df), {{gene_col}}, removeUnmappedRows = TRUE)
#'   gene_map <- gene_map[!duplicated(gene_map$STRING_id),]
#'   string_db$plot_network(gene_map$STRING_id)
#' }
#'
#'
#' #' Plot protein protein interaction network
#' #'
#' #' The function `plot_network()` is used to plot protein protein interaction (PPI) network.
#' #'
#' #' @param net_data A network data of class `tgxtool`. The network data must be a list which contain
#' #' two data frames `vertices` and `edges`. However,  the `vertices` must have column `p.values` and
#' #' `community`.
#' #' @param plot_layout
#' #' @param ... Optional arguments passed to [igraph::plot.igraph()].
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' plot_network <- function(net_data, plot_layout = "fr", rm_noedge = TRUE,arc = 0.5,
#'                          tag =  NULL, tag_size = rel(1), show_legend = TRUE)
#' {
#'   net <- igraph::graph_from_data_frame(d = net_data$edges, vertices = net_data$vertices, directed = TRUE)
#'   if (rm_noedge) {
#'     net <- igraph::delete.vertices(net, igraph::degree(net)==0)
#'   }
#'   col_ver <- c("#1F78B4", "#33A02C", "#E31A1C","#6A3D9A", "#FF00FF",
#'                "#8B4513", "#556B2F", "yellow3","#FF9F00", "paleturquoise4")
#'   # col_edg <- c("deepskyblue", "lightgreen", "#FB9A99","#CAB2D6", "#DDA0DD",
#'   #              "#CD853F", "#CAFF70", "#FFFF99", "#FDBF6F", "paleturquoise3")
#'   col_edg <- edge_color(col_ver)
#'   igraph::V(net)$color <- col_ver[igraph::V(net)$community]
#'   if(all(grepl(".*e-", igraph::V(net)$p_value))){
#'     igraph::V(net)$size <- ceiling(scales::rescale(
#'       as.numeric(gsub(".*e-","", igraph::V(net)$p_value)), to = c(3, 10)))
#'   } else {
#'     igraph::V(net)$size <- ceiling(scales::rescale(1/igraph::V(net)$p_value, to = c(3, 10)))
#'   }
#'   # as.numeric(gsub(".*e-","", igraph::V(net)$p.values))
#'   igraph::E(net)$width <- scales::rescale(igraph::E(net)$combined_score, to = c(0.1, 1))
#'   igraph::E(net)$color <- col_edg[igraph::V(net)$community][match(net_data$edges$from, igraph::V(net)$name)]
#'   pp <- ggraph::ggraph(net, layout = plot_layout) +
#'     ggraph::geom_edge_arc(aes(colour = color, width= width),strength = arc, alpha = 0.5) +
#'     ggraph::geom_node_point(size = igraph::V(net)$size,
#'                             color = igraph::V(net)$color) +
#'     ggraph::geom_node_text(aes(label = igraph::V(net)$name), size = 5, repel = TRUE) +
#'     labs(tag = tag)+
#'     theme(#text = element_text(),
#'       plot.tag = element_text(face = "bold", size = tag_size),
#'       panel.background = element_blank(),
#'       plot.background = element_rect(colour = NA),
#'       legend.position = "none"
#'     )
#'   class_name <- factor(igraph::V(net)$gene_class, levels = levels(net_data$vertices$gene_class))
#'   class_name <- droplevels(class_name)
#'   n_net <- length(levels(class_name))
#'   net_col <- col_ver[sort(unique(igraph::V(net)$community))]
#'   df_leg <- data.frame(x=1:n_net,y=1:n_net,
#'                        net_name= factor(levels(class_name), levels = levels(class_name)))
#'   leg_plot <- ggplot(df_leg,aes(x=x,y=y,color=net_name), size = rel(1))+
#'     geom_point()+
#'     scale_color_manual(values = net_col)+
#'     theme(
#'       legend.title = element_blank(),
#'       legend.position = "bottom",
#'       legend.direction = "horizontal",
#'       legend.text = element_text(size= 14),
#'       legend.margin = margin(0,unit = "mm"),
#'       legend.box.margin= margin(0,unit = "mm"),
#'       legend.spacing.y = unit(-0.3, 'mm'),
#'       legend.spacing.x = unit(1, 'mm'),
#'       legend.key=element_blank())+
#'     guides(color = guide_legend(override.aes = list(size = 6),
#'                                 nrow = dplyr::case_when(n_net <= 3~ 1,
#'                                                         n_net > 3 & n_net <= 6 ~ 2,
#'                                                         n_net > 6 ~ 3), byrow = TRUE))
#'   leg_net <- common_legend(leg_plot)
#'
#'   temp_df <- data.frame(x=rnorm(10),y=rnorm(10),z=runif(10,0.1,1),
#'                         g=sample(letters[1:2],10,replace = TRUE))
#'   leg_edgwt <- ggplot(temp_df) +
#'     geom_line(
#'       aes(x = x, y = y, linewidth = z, group = g),
#'       lineend = "round"
#'     ) +
#'     scale_linewidth("Weight: ", range = c(0.1, 2))+
#'     theme(
#'       legend.position = "bottom",
#'       legend.direction = "horizontal",
#'       legend.title = element_text(size= 14),
#'       legend.text = element_text(size= 12),
#'       legend.margin = margin(0,unit = "mm"),
#'       legend.box.margin= margin(0,unit = "mm"),
#'       legend.spacing.y = unit(-0.3, 'mm'),
#'       legend.spacing.x = unit(1, 'mm'),
#'       legend.key=element_blank())
#'   leg_wtsize <- common_legend(leg_edgwt)
#'
#'   net_logp <- ceiling(-log10(igraph::V(net)$p_value))
#'   net_size <- igraph::V(net)$size
#'   df_nsize <- data.frame(
#'     x= ceiling(c(min(net_logp), median(net_logp), max(net_logp))),
#'     y= 1:3,
#'     z= ceiling(c(min(net_size), median(net_size), max(net_size))))
#'   leg_pnet <- ggplot(df_nsize,aes(x=x,y=y,size=z), color = "gray50")+
#'     geom_point(aes(size=x))+
#'     labs(size = expression("-log"[10]*" (p-value): "))+
#'     theme(
#'       legend.position = "bottom",
#'       legend.direction = "horizontal",
#'       legend.title = element_text(size= 14),
#'       legend.text = element_text(size= 12),
#'       legend.margin = margin(0,unit = "mm"),
#'       legend.box.margin= margin(0,unit = "mm"),
#'       legend.spacing.y = unit(-0.3, 'mm'),
#'       legend.spacing.x = unit(1, 'mm'),
#'       legend.key=element_blank())
#'   leg_nsize <-common_legend(leg_pnet)
#'   if (show_legend) {
#'     ppp <- grid.arrange(pp,
#'                         gridExtra::arrangeGrob(leg_net,
#'                                                gridExtra::arrangeGrob(grobs = list(leg_nsize, leg_wtsize),
#'                                                                       heights = c(5,5), nrow = 2),
#'                                                ncol=2,widths=c(6, 4)),
#'                         nrow=2,heights=c(10, 1.5))
#'   } else {
#'     ppp <- pp
#'   }
#'   return(ppp)
#' }
#'
#' # End ppi plot
#'
#' #---------------------  function for enrich plot ----------------------
#' plot_enrich <- function(enrc_data, plot_layout = "fr",arc =0.3, gene_size = 5,
#'                         alpha = 0.4, title = NULL, title_size = rel(1),
#'                         tag = NULL, tag_size =rel(1)) {
#'   enrc_net <- igraph::graph_from_data_frame(
#'     d = enrc_data$edges,
#'     vertices = enrc_data$vertices,
#'     directed = TRUE
#'   )
#'   n_path <- sum(enrc_data$vertices$nodes == "pathway")
#'   n_gene <- sum(enrc_data$vertices$nodes == "gene")
#'   col_ver <- node_color(n_path)
#'   col_edg <- edge_color(col_ver, alpha = alpha)
#'
#'   ggcol <- colorRampPalette(c("azure1", "azure4"))
#'   logp <- ceiling(-log10(igraph::V(enrc_net)$p_value))[igraph::V(enrc_net)$nodes == "gene"]
#'  # maxp <- plyr::round_any(max(logp), 5, f = ceiling)
#'  # minp <- plyr::round_any(min(logp), -5, f = ceiling)
#'  # gene_col <- ggcol(maxp-minp)[logp-minp]
#'   igraph::V(enrc_net)$color <- c(col_ver[1:n_path], rep("gray", n_gene)) #gene_col
#'   igraph::E(enrc_net)$width <- floor(scales::rescale(igraph::E(enrc_net)$weight, to = c(1, 3)))
#'   igraph::E(enrc_net)$color <- col_edg[match(enrc_data$edges$from, igraph::V(enrc_net)$name)]
#'   pp <- ggraph::ggraph(enrc_net, layout = plot_layout) +
#'     ggraph::geom_edge_arc(aes(colour = color, width= 1.5),strength = arc, alpha = 0.5) +
#'     ggraph::geom_node_point(size = c(igraph::V(enrc_net)$size[1:n_path], rep(gene_size, n_gene)),
#'                             color = igraph::V(enrc_net)$color) +
#'     ggraph::geom_node_text(aes(label = c(rep("", n_path),attr(igraph::V(enrc_net), "name")[(n_path +1):(n_path + n_gene)])),
#'                            size = gene_size, repel = TRUE) +
#'     labs(title = title,tag = tag) +
#'     theme(plot.title = element_text(face = "bold", size = rel(title_size),hjust = 0.5),
#'           plot.tag = element_text(face = "bold", size = tag_size),
#'           panel.background = element_blank(),
#'           plot.background = element_rect(colour = NA),
#'           legend.position = "none"
#'     )
#'
#'   path_names <- split_merge(igraph::V(enrc_net)$name[enrc_data$vertices$nodes == "pathway"], 25)
#'   path_name <- factor(path_names, levels= path_names)
#'   path_size <- igraph::V(enrc_net)$size[enrc_data$vertices$nodes == "pathway"]
#'   path_col <- igraph::V(enrc_net)$color[enrc_data$vertices$nodes == "pathway"]
#'   df_leg <- data.frame(x=1:n_path,y=1:n_path, path_name= path_name)
#'   leg_plot <- ggplot(df_leg,aes(x=x,y=y,color=path_name, size = path_name))+
#'     geom_point()+
#'     scale_color_manual(values=path_col)+
#'     scale_size_manual(values =path_size/2)+
#'     theme(
#'       legend.title = element_blank(),
#'       # legend.position = "bottom",
#'       #legend.direction = "horizontal",
#'       legend.text = element_text(size= 13),
#'       legend.margin = margin(0,unit = "mm"),
#'       legend.box.margin= margin(0,unit = "mm"),
#'       legend.spacing.y = unit(-0.3, 'mm'),
#'       legend.spacing.x = unit(3, 'mm'),
#'       legend.key=element_blank())+
#'     guides(color = guide_legend(ncol = ifelse(n_path > 5, 2, 1), byrow = TRUE))
#'   leg_path <- common_legend(leg_plot)
#'
#'   # df_gene <- data.frame(x=logp,y=1:length(logp))
#'   # leg_gplot <- ggplot(df_gene,aes(x=x,y=y,color=x))+
#'   #   geom_point()+
#'   #   labs(color = "")+
#'   #   scale_colour_gradient(low = "azure1", high = "azure4")
#'   # leg_gene <-common_legend(leg_gplot)
#'
#'   path_logp <- ceiling(-log(igraph::V(enrc_net)$p_value[enrc_data$vertices$nodes == "pathway"]))
#'   df_psize <- data.frame(
#'     x= ceiling(c(min(path_logp), median(path_logp), max(path_logp))),
#'     y= 1:3,
#'     z= ceiling(c(min(path_size/2), median(path_size/2), max(path_size/2))))
#'   leg_ppath <- ggplot(df_psize,aes(x=x,y=y,size=z), color = "gray50")+
#'     geom_point(aes(size=x))+
#'     labs(size = expression("-log"[10]*" (p-value)"))+
#'     theme(
#'       legend.text = element_text(size= 10),
#'       legend.margin = margin(0,unit = "mm"),
#'       legend.box.margin= margin(0,unit = "mm"),
#'       legend.spacing.y = unit(-0.3, 'mm'),
#'       legend.spacing.x = unit(1, 'mm'),
#'       legend.key=element_blank())
#'   leg_psize <-common_legend(leg_ppath)
#'
#'   ppp <- grid.arrange(pp,
#'                       gridExtra::arrangeGrob(leg_path,leg_psize,
#'                                             # gridExtra::arrangeGrob(list(leg_psize, leg_gene),
#'                                             #                        widths = c(5,5), ncol = 2,
#'                                             #                        top = textGrob(expression("-log"[10]*" (p-value)"))),
#'                                              ncol=2,widths=c(10, 2)),
#'                       nrow=2,heights=c(10, 3.5))
#'
#'   return(ppp)
#' }
#'
#' # -End of enrichment plot
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' col25 <- c(
#'   "dodgerblue2", "#E31A1C", # red
#'   "green4",
#'   "#6A3D9A", # purple
#'   "#FF7F00", # orange
#'   "black", "gold1",
#'   "skyblue2", "#FB9A99", # lt pink
#'   "palegreen2",
#'   "#CAB2D6", # lt purple
#'   "#FDBF6F", # lt orange
#'   "gray70", "khaki2",
#'   "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
#'   "darkturquoise", "green1", "yellow4", "yellow3",
#'   "darkorange4", "brown"
#' )
#'
#'
#' plot_aopnet <- function(aop_data, plot_layout = "kk",arc = 0.5, edge_alpha = 0.5,
#'                         title = NULL, title_size = rel(1),
#'                         tag = NULL, tag_size =rel(1)) {
#'   aop_net <- igraph::graph_from_data_frame(
#'     d = aop_data$edges,
#'     vertices = aop_data$vertices,
#'     directed = TRUE
#'   )
#'   n_path <- sum(aop_data$vertices$nodes != "gene")
#'   n_gene <- sum(aop_data$vertices$nodes == "gene")
#'   path <- aop_data$vertices$nodes[aop_data$vertices$nodes != "gene"]
#'   path_class <- unique(aop_data$vertices$nodes[aop_data$vertices$nodes != "gene"])
#'   node_colors <- node_color(n_path)
#'   edge_colors <- edge_color(node_colors, alpha = edge_alpha)
#'   category_color <- colorRamps::primary.colors(n_path)
#'   igraph::V(aop_net)$fill <-  c(node_colors, rep("#D3D3D3", n_gene))
#'   igraph::V(aop_net)$color <-  c(category_color[match(path, path_class)], rep("#D3D3D3", n_gene))
#'   igraph::E(aop_net)$width <- round(scales::rescale(igraph::E(aop_net)$weight, to = c(0.1, 3)),2)
#'   # igraph::E(aop_net)$color <- "#BC8F8F"
#'   igraph::E(aop_net)$color <- edge_colors[match(aop_data$edges$from, igraph::V(aop_net)$name)]
#'
#'   pp <- ggraph::ggraph(aop_net, layout = plot_layout) +
#'     ggraph::geom_edge_arc(aes(color = color, width= width),strength = arc, alpha = 0.3) +
#'     ggraph::geom_node_point(size = c(igraph::V(aop_net)$size[1:n_path], rep(4, n_gene)),
#'                             color = igraph::V(aop_net)$color,
#'                             fill = igraph::V(aop_net)$fill, shape = 21, stroke = 2) +
#'     ggraph::geom_node_text(aes(label = c(rep("", n_path),attr(igraph::V(aop_net), "name")[(n_path +1):(n_path + n_gene)])), size = 5, repel = TRUE) +
#'     labs(title = "",tag = "") +
#'     theme(plot.title = element_text(face = "bold", size = rel(1),hjust = 0.5),
#'           plot.tag = element_text(face = "bold", size = 2),
#'           panel.background = element_blank(),
#'           plot.background = element_rect(colour = NA),
#'           legend.position = "none"
#'     )
#'
#'   path_name <- factor(igraph::V(aop_net)$name[aop_data$vertices$nodes != "gene"],
#'                       levels= igraph::V(aop_net)$name[aop_data$vertices$nodes != "gene"])
#'   path_fill <- igraph::V(aop_net)$fill[aop_data$vertices$nodes != "gene"]
#'   path_col <- igraph::V(aop_net)$color[aop_data$vertices$nodes != "gene"]
#'   df_leg <- data.frame(x = 1:n_path,y = 1:n_path, path_name = split_merge(path_name, 25))
#'   leg_plot <- ggplot(df_leg,aes(x=x,y=y,color=path_name, fill = path_name))+
#'     geom_point(shape =21, stroke = 2)+
#'     scale_color_manual(values=path_col)+
#'     scale_fill_manual(values = path_fill)+
#'     theme(
#'       legend.title = element_blank(),
#'       legend.text = element_text(size= 14),
#'       legend.margin = margin(0,unit = "mm"),
#'       legend.box.margin= margin(0,unit = "mm"),
#'       legend.spacing.y = unit(-0.3, 'mm'),
#'       legend.spacing.x = unit(3, 'mm'),
#'       legend.key=element_blank())+
#'     #guides(color = guide_legend(ncol = ifelse(n_path > 10, 2, 1), byrow = TRUE))
#'     guides(color = guide_legend(override.aes = list(size= 5), ncol = 1))
#'   leg_path <- common_legend(leg_plot)
#'
#'   temp_df <- data.frame(x=rnorm(100),y=rnorm(100),z=runif(100,0.1, 1),
#'                         g=sample(letters[1:2],10,replace = TRUE))
#'   leg_edgwt <- ggplot(temp_df) +
#'     geom_line(
#'       aes(x = x, y = y, linewidth = z, group = g),
#'       lineend = "round"
#'     ) +
#'     scale_linewidth("Edges: ", range = c(0.1, 3))+
#'     theme(
#'       legend.position = "bottom",
#'       legend.direction = "horizontal",
#'       legend.title = element_text(size= 14),
#'       legend.text = element_text(size= 14),
#'       legend.margin = margin(0,unit = "mm"),
#'       legend.box.margin= margin(0,unit = "mm"),
#'       legend.spacing.y = unit(-0.3, 'mm'),
#'       legend.spacing.x = unit(1, 'mm'),
#'       legend.key=element_blank())
#'   leg_wtsize <- common_legend(leg_edgwt)
#'
#'   ppp <- grid.arrange(gridExtra::arrangeGrob(grobs = list(pp, leg_wtsize),
#'                                              heights = c(9.5,0.5), nrow = 2),
#'                       leg_path,
#'                       ncol = 2, widths = c(7, 3))
#'
#'   return(ppp)
#' }
#'
#'
#'
#' plot_hubaop <- function(hub_netdata, plot_layout = "fr",gene_size = 5,
#'                         edge_alpha = 0.5, title = NULL, title_size = rel(1),
#'                         tag = NULL, tag_size =rel(1)) {
#'   hub_net <- igraph::graph_from_data_frame(
#'     d = hub_netdata$edges,
#'     vertices = hub_netdata$vertices,
#'     directed = TRUE
#'   )
#'   n_path <- sum(hub_netdata$vertices$nodes == "aop")
#'   n_gene <- sum(hub_netdata$vertices$nodes == "gene")
#'   col_ver <- c("#1F78B4", "#33A02C", "#E31A1C", "#FF9F00", "#FF00FF",
#'                "#8B4513", "#556B2F", "#6A3D9A", "yellow3", "paleturquoise4")
#'   col_edg <- c("deepskyblue", "lightgreen", "#FB9A99", "#FDBF6F", "#DDA0DD",
#'                "#CD853F", "#CAFF70", "#CAB2D6", "#FFFF99", "paleturquoise3")
#'   # col_ver <- node_color(n_path)
#'   # col_edg <- edge_color(col_ver, alpha = edge_alpha)
#'
#'   # ggcol <- colorRampPalette(c("azure1", "azure4"))
#'   # logp <- ceiling(-log10(igraph::V(enrc_net)$p_value))[igraph::V(enrc_net)$nodes == "gene"]
#'   # maxp <- plyr::round_any(max(logp), 5, f = ceiling)
#'   # minp <- plyr::round_any(min(logp), -5, f = ceiling)
#'   # gene_col <- ggcol(maxp-minp)[logp-minp]
#'   igraph::V(hub_net)$color <- c(col_ver[1:n_path],  rep("#D3D3D3", n_gene))
#'   igraph::E(hub_net)$width <- round(scales::rescale(igraph::E(hub_net)$weight, to = c(0.1, 3)),2)
#'   igraph::E(hub_net)$color <- col_edg[match(hub_netdata$edges$from, igraph::V(hub_net)$name)]
#'
#'   pp <- ggraph::ggraph(hub_net, layout = plot_layout) +
#'     ggraph::geom_edge_arc(aes(colour = color, width= width),strength = 0.5, alpha = 0.3) +
#'     ggraph::geom_node_point(size = c(igraph::V(hub_net)$size[1:n_path], rep(gene_size, n_gene)),
#'                             color = igraph::V(hub_net)$color) +
#'     ggraph::geom_node_text(aes(label = c(rep("", n_path),attr(igraph::V(hub_net), "name")[(n_path +1):(n_path + n_gene)])),
#'                            size = gene_size, repel = TRUE) +
#'     labs(title = title,tag = tag) +
#'     theme(plot.title = element_text(face = "bold", size = rel(title_size),hjust = 0.5),
#'           plot.tag = element_text(face = "bold", size = tag_size),
#'           panel.background = element_blank(),
#'           plot.background = element_rect(colour = NA),
#'           legend.position = "none"
#'     )
#'
#'   path_name <- factor(igraph::V(hub_net)$name[hub_netdata$vertices$nodes == "aop"],
#'                       levels= igraph::V(hub_net)$name[hub_netdata$vertices$nodes == "aop"])
#'   path_size <- igraph::V(hub_net)$size[hub_netdata$vertices$nodes == "aop"]
#'   path_col <- igraph::V(hub_net)$color[hub_netdata$vertices$nodes == "aop"]
#'   df_leg <- data.frame(x=1:n_path,y=1:n_path, path_name= path_name)
#'   leg_plot <- ggplot(df_leg,aes(x=x,y=y,color=path_name, size = path_name))+
#'     geom_point()+
#'     scale_color_manual(values=path_col)+
#'     scale_size_manual(values =path_size/2)+
#'     theme(
#'       legend.title = element_blank(),
#'       legend.position = "bottom",
#'       legend.direction = "horizontal",
#'       legend.text = element_text(size= 14),
#'       legend.margin = margin(0,unit = "mm"),
#'       legend.box.margin= margin(0,unit = "mm"),
#'       legend.spacing.y = unit(-0.3, 'mm'),
#'       legend.spacing.x = unit(3, 'mm'),
#'       legend.key=element_blank())+
#'     guides(color = guide_legend(ncol = ifelse(n_path > 3, 2, 1), byrow = TRUE))
#'   leg_path <- common_legend(leg_plot)
#'
#'   # temp_df <- data.frame(x=rnorm(10),y=rnorm(10),z=runif(10,0.1, 1),
#'   #                       g=sample(letters[1:2],10,replace = TRUE))
#'   # leg_edgwt <- ggplot(temp_df) +
#'   #   geom_line(
#'   #     aes(x = x, y = y, linewidth = z, group = g),
#'   #     lineend = "round"
#'   #   ) +
#'   #   scale_linewidth("Edges: ", range = c(0.1, 2))+
#'   #   theme(
#'   #     legend.position = "bottom",
#'   #     legend.direction = "horizontal",
#'   #     legend.text = element_text(size= 10),
#'   #     legend.margin = margin(0,unit = "mm"),
#'   #     legend.box.margin= margin(0,unit = "mm"),
#'   #     legend.spacing.y = unit(-0.3, 'mm'),
#'   #     legend.spacing.x = unit(1, 'mm'),
#'   #     legend.key=element_blank())
#'   # leg_wtsize <- common_legend(leg_edgwt)
#'
#'   ppp <- grid.arrange(pp,leg_path,
#'                       nrow=2, heights=c(8, 2))
#'   # gridExtra::arrangeGrob(grobs = list(, leg_wtsize), widths = c(6,4), ncol = 2)
#'   return(ppp)
#' }
#'
#' gene_data <- read_tsv("D:/ToxAssay/aop_data/compound_gene.tsv")
#' dise_voc <- read_csv("D:/ToxAssay/ToxAssay_genes/liver_diseases.csv")
#' dise_data <- read_tsv("D:/ToxAssay/aop_data/compound_disease.tsv")
#'
#' ctd <- get_ctd(compounds = comps_gr$Positive)
#'
#' mm <- get_transaction(v1 = ctd$compound_gene,
#'                       v2 = ctd$compound_disease,
#'                       compounds = NULL,
#'                       gene_id = res4$gene_name,
#'                       disease_id = dise_voc$DiseaseID)
#' oo <- get_aops(mm[[1]], mm[[2]], mm[[3]])
#' aop <- aop_network(oo, dise_voc, lift)
#'
#' plot_aopnet(aop, plot_layout = "nicely", node_color = col25,arc = 0.5)
#'
#' aop_res <- oo %>%
#'   filter(rhs %in% dise_voc$DiseaseID) %>%
#'   mutate(Disease = dise_voc$DiseaseName[match(rhs, dise_voc$DiseaseID)],
#'          Group = dise_voc$DiseaseGroup[match(rhs, dise_voc$DiseaseID)], .before = support) %>%
#'   arrange(desc(oddsRatio))
#'
#'
#'
#' network <- igraph::graph_from_data_frame(d = aop$edges, vertices = aop$vertices, directed = TRUE)
#' deg <- igraph::degree(network)            # Degree centrality
#' clo <- igraph::closeness(network)         # Closeness centrality
#' bet <- igraph::betweenness(network)       # Betweenness centrality
#' eig <- igraph::evcent(network)$vector     # Eigenvector centrality
#'
#'
#' hub_aop <- sort(deg[deg>=30], decreasing = TRUE)
#' hub_aop <- sort(eig[eig>0.7], decreasing = TRUE)
#'
#'
#' hub_edge <- aop$edges %>%
#'   dplyr::filter(from %in% names(hub_aop))
#'
#' hub_vertics <- aop$vertices %>%
#'   dplyr::filter(name %in% c(unique(hub_edge$from), unique(hub_edge$to))) %>%
#'   dplyr::mutate(nodes = ifelse(nodes == "gene", "gene", "aop"))
#'
#' hub_netdata <- list(vertices = hub_vertics, edges = hub_edge)
#'
#' plot_hubaop(hub_netdata, plot_layout = "graphopt") # graphopt
#'
#' unique(hub_edge$to)
#'
#'
#'
#'
#' oo %>%
#'   filter(rhs == "MESH:D056487") %>%
#'   arrange(desc(oddsRatio))
#'
#' ###------------------------------
#' df_ap2 <- oo %>%
#'   filter(lhs == "TXNRD1") %>%
#'   filter(rhs %in% dise_voc$DiseaseID) %>%
#'   mutate(Disease = dise_voc$DiseaseName[match(rhs, dise_voc$DiseaseID)], .before = support) %>%
#'   arrange(desc(oddsRatio))
#'
#'
#' df_ap3 <- oo %>%
#'   filter(rhs %in% dise_voc$DiseaseID) %>%
#'   mutate(Disease = dise_voc$DiseaseName[match(rhs, dise_voc$DiseaseID)], .before = support) %>%
#'   arrange(desc(lift)) %>%
#'   group_by(rhs, Disease) %>%
#'   summarise(total_lift = sum(count)) %>%
#'   arrange(desc(total_lift))
#' nest() %>%
#'   mutate(`Genes (OR)` = map_chr(data,~glue_collapse(
#'     glue("{block_fst(.x$lhs)} ({round(.x$oddsRatio,2)})"),", ", last = " and "))) %>%
#'   select(-data)
#'
#' write.csv(df_ap3, file = "G:/Thesis/Thesis/thesis_writing/Tables/dise_table41genes.csv")
#'
#'
#' length(unique(dese_net$rhs))
#'
#'
#' hub_dis <- dese_net %>%
#'   filter(Group %in% c("Cirrhosis", "Degeneration", "Encephalopathy", "Fatty Liver", "Hepatitis"))%>%
#'   filter(confidence >= 0.8)
#'
#' unique(hub_dis$lhs)
#'
#'
#' unique(dese_net$Group)
#'
#' write.csv(dese_net, file = "G:/Thesis/Thesis/thesis_writing/Tables/dise_net41genes.csv")
#'
#'
#'
#' #----------------------- start from here ---------------------
#' dese_net2 <- read.csv("E:/Thesis//thesis_writing/Tables/dise_net.csv")
#'
#'
#' unique(aop_data$edges$from)
#'
#'
#'
#' ggsave(filename = "G:/Thesis/Thesis/thesis_writing/Figure/fig_dese41genes.jpeg", plot = p, width = 350, height = 200, dpi = 500, units = "mm")
#'
#'
#' jpeg(filename = "G:/Thesis/Thesis/thesis_writing/Figure/fig_dese.jpeg", width = 340, height = 180, res = 500, units = "mm")
#'
#' ppp <- grid.arrange(pp,leg_path,
#'                     ncol=2, widths=c(10, 3.5))
#' dev.off()
#'
#'
#'
#'
#' #------------------------------
#'
#'
#'
#' ggplot(df_cirosis, aes(oddsRatio, reorder(genes, oddsRatio), fill = lift)) +
#'   geom_bar(stat = "identity")+
#'   labs(x="Number of genes", y= "Genes", fill = "lift")+
#'   #scale_fill_continuous(breaks = c(5,12,20))+
#'   theme_Publication(tag_size = rel(1.5),
#'                     leg_size = rel(0.8),
#'                     lab_size = rel(0.9)
#'                     #x_lab = element_text(size = rel(0.6)),
#'                     #y_lab = element_text(size = rel(0.6))
#'   )+
#'   theme(legend.direction = "horizontal",
#'         legend.position = "bottom")
#'
#' #------------------------
#'
#' lift_RR <- function(aop_data){
#'   support <- aop_data$support
#'   confidence <- aop_data$confidence
#'   lift <- aop_data$lift
#'   prob_E <- support / confidence 					            	#P(E)
#'   prob_O <- confidence / lift                            				#P(O)
#'   RR <- (1 - prob_E)*lift / (1 - prob_E*lift)
#'   OR <- (RR - confidence) / (1 - confidence)
#'   lift_O_negative <- (1 - confidence) / (1 - prob_O) 		#Lift(O~|E)
#'   ul_RR <- (1 - prob_E)*aop_data$CI.UL / (1 - prob_E*aop_data$CI.UL)
#'   ul_OR <- (ul_RR - confidence) / (1 - confidence)
#'   result <- data.frame(aop_data$lhs, aop_data$rhs, lift, RR, OR, ul_RR, ul_OR, aop_data$oddsRatio)
#'   return(result)
#' }
#' lift_RR(aop_data)
#' aop_data$CI.UL
#'
#'
#' #-------------------- SHAP value
#' shap_plot <- function (data_long,
#'                        x_bound = NULL,
#'                        hlt_feature = NULL,
#'                        scientific = FALSE,
#'                        legend_position = "right",
#'                        heading= ""
#' )
#' {
#'   if (scientific) {
#'     label_format = "%.1e"
#'   } else {
#'     label_format = "%.3f"
#'   }
#'   x_bound <- if (is.null(x_bound)) {
#'     max(abs(data_long$value)) * 1.1
#'   } else as.numeric(abs(x_bound))
#'
#'   if(is.null(hlt_feature)) {
#'     data_long <- data_long %>% mutate(
#'       name = fct_reorder(variable, mean_value))
#'   } else {
#'     data_long <- data_long %>%
#'       dplyr::mutate(genes = ifelse(variable %in% hlt_feature, "highlited", "others"))
#'     group_colors <- c(highlited = 'black', others = "gray")
#'     data_long <- data_long %>% mutate(
#'       name = glue("<p style='color:{group_colors[genes]}'>{variable}</p>"),
#'       name = fct_reorder(name, mean_value))
#'   }
#'   if (identical(legend_position, "right")) {
#'     bar_width <- 0.4
#'     bar_height <- 10
#'   } else {
#'     bar_width <- 10
#'     bar_height <- 0.4
#'   }
#'   plot_shap <- ggplot(data = data_long) + coord_flip(ylim = c(-x_bound, x_bound)) +
#'     geom_hline(yintercept = 0) +
#'     ggforce::geom_sina(aes(x = reorder(name, value), y = value, color = stdfvalue),
#'                        method = "counts", maxwidth = 0.7,  alpha = 0.7) +
#'     geom_text(data = unique(data_long[, c("name", "mean_value")]),
#'               aes(x = name, y = -Inf, label = sprintf(label_format, mean_value)),
#'               size = 4, alpha = 0.7, hjust = -0.2, fontface = "bold", check_overlap = TRUE) +
#'     scale_color_gradient(low = "yellow", high = "red",
#'                          breaks = c(0, 1), labels = c(" Low", "High "),
#'                          guide = guide_colorbar(barwidth = bar_width, barheight = bar_height)) +
#'     theme_bw() + theme(plot.title = element_text(size = rel(1.3), face = "bold", hjust = 0.5),
#'                        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
#'                        axis.text.y = ggtext::element_markdown(size = rel(1.3)),
#'                        # legend.justification = "top",
#'                        legend.position = legend_position,
#'                        legend.text = element_text(size = 8), axis.title.x = element_text(size = 10)) +
#'     scale_x_discrete(limits = levels(data_long$name),
#'                      labels = levels(data_long$name)) +
#'     labs(title = heading, y = "SHAP value (impact on model output)", x = "",
#'          color = "Expression    ")
#'   return(plot_shap)
#' }
#' #---------------------  END shap
#'
#'
#'
#' #----------  Start AOP plot     ----------------
#' plot_aop <- function(aop_res,
#'                      gene_order = NULL,
#'                      disease_order = NULL,
#'                      n_gclass = 3,
#'                      n_dclass = 3,
#'                      genename_block = FALSE,
#'                      adj_label = FALSE) {
#'   gene_sim <- gene_similarity(aop_res, disease_column = rhs, target_column = lhs)
#'   disease_sim <- disease_similarity(aop_res, gene_column = lhs, target_column = rhs)
#'   hclust_gene <- hclust(dist(gene_sim))
#'   cluster_gene <- cutree(hclust_gene, k = n_gclass)
#'   hclust_disease <- hclust(dist(disease_sim))
#'   cluster_disease <- cutree(hclust_disease, k = n_dclass)
#'   if(is.null(gene_order)) {
#'     gene_order <- names(sort(cluster_gene))
#'   }
#'   if(is.null(disease_order)) {
#'     disease_order <-  names(sort(cluster_disease))
#'   }
#'   aop_df <- aop_res %>%
#'     dplyr::select(c(lhs, rhs, lift)) %>%
#'     dplyr::mutate(disease_gr = cluster_disease[match(rhs, names(cluster_disease))],
#'                   gene_gr = n_dclass + cluster_gene[match(lhs, names(cluster_gene))])
#'   if (length(gene_order) > length(disease_order)) {
#'     gene_pos <- data.frame(label = gene_order, x = 1 : length(gene_order))
#'     add_incmt <- floor((length(gene_order) - length(disease_order))/2)
#'     disease_pos <- data.frame(label = disease_order, x = add_incmt + (1 : length(disease_order)))
#'   } else if (length(gene_order) < length(disease_order)) {
#'     add_incmt <- round((length(disease_order) - length(gene_order))/2)
#'     gene_pos <- data.frame(label = gene_order, x = add_incmt + (1 : length(gene_order)))
#'     disease_pos <- data.frame(label = disease_order, x = 1 : length(disease_order))
#'   } else {
#'     gene_pos <- data.frame(label = gene_order, x = 1 : length(gene_order))
#'     disease_pos <- data.frame(label = disease_order, x = 1 : length(disease_order))
#'   }
#'   aop_df <- merge(aop_df, gene_pos, by.x = "lhs", by.y = "label")
#'   aop_df <- merge(aop_df, disease_pos, by.x = "rhs", by.y = "label", suffixes = c("_gene", "_disease"))
#'   aop_df <- aop_df %>%
#'     dplyr::mutate(rhs = gsub("MESH:", "", rhs),
#'                   gene = "gene", disease = "disease",
#'                   relation = paste0("rel", 1:nrow(aop_df)),
#'                   disease_hjust = -0.2,
#'                   gene_hjust = 1.2)
#'   if (!genename_block) {
#'     aop_df <- aop_df %>%
#'       dplyr::mutate(lhs = block_fst(lhs))
#'   }
#'
#'   n_node <- sort(c(unique(aop_df$gene_gr), unique(aop_df$disease_gr)))
#'   node_colors <- node_color(length(n_node))
#'   #edge_colors <- edge_color(node_colors)
#'   data_long <- data.frame(
#'     nodes =  as.vector(t(aop_df[c("lhs", "rhs")])),
#'     y_idx = as.vector(t(aop_df[c("x_gene", "x_disease")])),
#'     x_var = factor(as.vector(t(aop_df[c("gene", "disease")])), levels = c("gene", "disease")),
#'     gr = rep(aop_df$relation, each = 2),
#'     #edge_col = as.vector(edge_colors[rep(aop_df$disease_gr, each = 2)]),
#'     edge_col = as.vector(node_colors[rep(aop_df$disease_gr, each = 2)]),
#'     node_col = as.vector(node_colors[as.vector(t(aop_df[c("gene_gr", "disease_gr")]))]),
#'     score = rep(aop_df$lift, each = 2),
#'     adj_hjust = 0.5
#'   )
#'   if (adj_label) {
#'     data_long <- data_long %>%
#'       mutate(adj_hjust = as.vector(t(aop_df[c("gene_hjust", "disease_hjust")])))
#'   }
#'   link_plot <- ggplot(data_long, aes(x = x_var, y = y_idx, group = gr)) +
#'     geom_line(aes(linewidth = score, color = edge_col), alpha = 0.3) +
#'     geom_point(aes(color = node_col) ,shape = 16,size = 4)+
#'     geom_text(aes(label = nodes, hjust = adj_hjust), vjust = 0.4,
#'               size = 3, color = "gray20", parse = TRUE)+
#'     scale_size_continuous(range = c(0.1, 3)) +
#'     scale_y_continuous(breaks = NULL) +  # Remove y-axis breaks for simplicity
#'     theme_void() +
#'     theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 2)) +
#'     labs(x = "", y = "", title = "") +
#'     theme(legend.position = "none")
#'   return(link_plot)
#' }
#' #-----------------------PPI corediagram plot-------------
#' library(circlize)
#'
#' net_enrichment <- function(net_gene = net_fgreedy$vertices,
#'                            path_gene = path_gene,
#'                            header_size = 0.6,
#'                            header_split = 12,
#'                            starting_angle = 0, ...) {
#'   join_data <- right_join(net_gene, path_gene, by="gene_symbol")
#'
#'   df <- join_data %>%
#'     dplyr::mutate(path = as.factor(path)) %>%
#'     dplyr::group_by(gene_class, path) %>%
#'     dplyr::summarise(value = n()) %>%
#'     dplyr::rename(from = gene_class, to = path)
#'   col_ver <- c("#1F78B4", "#33A02C", "#E31A1C","#6A3D9A", "#FF00FF",
#'                "#8B4513", "#556B2F", "yellow3","#FF9F00", "paleturquoise4")
#'   path_name <- unique(path_gene$path)
#'   cl_name <- levels(net_gene$gene_class)
#'   grid_col <- c(col_ver[1:length(cl_name)], rep("gray",length(path_name)))
#'   names(grid_col) <- c(cl_name, path_name)
#'   p1 <- function() {
#'     par(mar = c(0.2, 0.2, 0.2, 0.2), cex = 1)
#'     circlize::circos.clear()
#'     circlize::circos.par(start.degree = starting_angle,
#'                          track.margin = c(-0.1, 0.1),
#'                          gap.after = 1, #c(rep(2, length(cl_name)-1), 15, rep(2, length(path_name)-1), 15),
#'                          points.overflow.warning = FALSE
#'     )
#'     # plot the chord diagram
#'     circlize::chordDiagram(x = df,
#'                            directional = 1,
#'                            big.gap = 5,
#'                            order = sort(c(cl_name, path_name)),
#'                            scale = TRUE,
#'                            grid.col = grid_col,
#'                            annotationTrack = "grid",
#'                            transparency = 0.25,
#'                            annotationTrackHeight = c(0.05, 0.1),
#'                            direction.type = c("diffHeight", "arrows"),
#'                            link.arr.type = "big.arrow",
#'                            diffHeight  = -0.04,
#'                            link.sort = TRUE,
#'                            link.largest.ontop = TRUE
#'     )
#'     # add labels and axis
#'     circlize::circos.track(track.index = 1,
#'                            bg.border = NA,
#'                            panel.fun = function(x, y) {
#'
#'                              sector_name = circlize::get.cell.meta.data("sector.index")
#'                              xlim = circlize::get.cell.meta.data("xlim")
#'                              ylim = circlize::get.cell.meta.data("cell.ylim")
#'                              # header_split <-  15 #round(100* split_ratio/n_words(sector_name))
#'                              split_name <- rev(split_sentence(sector_name, header_split))
#'                              for (i in 1:length(split_name)) {
#'                                circlize::circos.text(x = mean(xlim),
#'                                                      y = ylim[2] + (ylim[2] - ylim[1]) * 0.1,
#'                                                      labels = split_name[[i]],
#'                                                      cex = header_size ,
#'                                                      font = 1,
#'                                                      facing = "bending.inside",
#'                                                      adj = c(0.5, -i*1.3),
#'                                                      niceFacing = FALSE
#'                                )
#'                              }
#'                            })
#'
#'     # circlize::circos.track(track.index = 1,
#'     #                        bg.border = NA,
#'     #                        panel.fun = function(x, y) {
#'     #
#'     #                          sector.name = circlize::get.cell.meta.data("sector.index")
#'     #                          split_name <- split_merge(sector.name, 20)
#'     #                          name1 <- split_name[[1]]
#'     #                          name2 <- split_name[[2]]
#'     #
#'     #                          #paste(unlist(str_split(circlize::get.cell.meta.data("sector.index"), " ", n=2)), collapse = "\n ")
#'     #
#'     #                          #sapply(circlize::get.cell.meta.data("sector.index"), function(x) paste(unlist(str_split(x, " ", n=2)), collapse = "\n "), USE.NAMES = FALSE)
#'     #
#'     #                          xlim = circlize::get.cell.meta.data("xlim")
#'     #                          ylim = circlize::get.cell.meta.data("ylim")
#'     #                          circlize::circos.text(x = mean(xlim),
#'     #                                                y = ylim[1] + 2.8,
#'     #                                                labels = name1,
#'     #                                                cex = lab_size ,
#'     #                                                font = 1,
#'     #                                                facing = "bending.inside",
#'     #                                                adj = c(0.5, 0),
#'     #                                                niceFacing = FALSE
#'     #                                                )
#'     #                          circlize::circos.text(x = mean(xlim),
#'     #                                                y = ylim[1] + 2,
#'     #                                                labels = name2,
#'     #                                                cex = lab_size ,
#'     #                                                font = 1,
#'     #                                                facing = "bending.inside",
#'     #                                                adj = c(0.5, 0),
#'     #                                                niceFacing = FALSE
#'     #                                                )
#'     #                          # circos.text(mean(xlim), ylim[1] + .1, sector.name,
#'     #                          #             facing = "bending", cex = 0.8)
#'     #                          #  circos.text(x = mean(xlimt), y = 3, labels = reg2, facing = "bending", cex = 0.8)
#'     #                          circlize::circos.axis(h = "top",
#'     #                                                labels.cex = 0.5,
#'     #                                                labels.niceFacing = FALSE,
#'     #                                                labels.pos.adjust = FALSE
#'     #                                                )
#'     #                        }
#'     #   )
#'     title_tag(...)
#'     circlize::circos.clear()
#'   }
#'   grob_obj <- cowplot::as_grob(p1)
#'   p <- cowplot::ggdraw(grob_obj)
#'   return(p)
#' }
#'
#'
#' #----------------------------  circlize Heatmap ------------------
#' subnet_heatmap <- function(..., probes, expr_data, attr_data, net_data,
#'                            lab_size = 1,axis_size = 0.8,starting_angle = 0,
#'                       outer_dose = "High", outer_time = "24 hr", inner_dose = "High", inner_time = "9 hr",
#'                       outer_title = "Outer circle", inner_title = "Inner circle",
#'                       labs_option = list( title = NULL, tag = NULL, tag_size = 1.5,
#'                                           note = NULL, title_size = 2, note_size = 0.5)) {
#'   df_link <- net_data$edges
#'   net_data <- net_data$vertices
#'   comps_gr <- test_group(...)
#'   temp_mat1 <- get_subset(
#'     comps_gr,
#'     ge_matrix = expr_data,
#'     metadata = attr_data,
#'     probes = probes,
#'     dose = outer_dose,
#'     time = outer_time
#'   )
#'   mat1 <- temp_mat1$expression
#'   temp_mat2 <- get_subset(
#'     comps_gr,
#'     ge_matrix = expr_data,
#'     metadata = attr_data,
#'     probes = probes,
#'     dose = inner_dose,
#'     time = inner_time
#'   )
#'   mat2 <- temp_mat2$expression
#'   mat1 <- mat1[match(net_data$probe_id , rownames(mat1)),]
#'   mat2 <- mat2[match(net_data$probe_id , rownames(mat1)),]
#'   split <- net_data$gene_class[match(rownames(mat1), net_data$probe_id)]
#'   gene_mat <- net_data$gene_symbol[match(rownames(mat1), net_data$probe_id)]
#'   #split <- split[!is.na(split)]# split[is.na(split)]= "Miscellaneous"
#'   rownames(mat1) <- gene_mat
#'   rownames(mat2) <- gene_mat
#'
#'   #######################
#'   col_ver <- c("#1F78B4", "#33A02C", "#E31A1C","#6A3D9A", "#FF00FF",
#'                "#8B4513", "#556B2F", "yellow3","#FF9F00", "paleturquoise4","gray")
#'   gene_col <- net_data$community[match(gene_mat, net_data$gene_symbol)]
#'   a <- net_data$community[match(df_link$from, net_data$gene_symbol)]
#'   b <- net_data$community[match(df_link$from, net_data$gene_symbol)]
#'   d <- ifelse(a == b, a, 11)
#'   idx_link = data.frame(
#'     from = match(df_link$from, gene_mat),
#'     to = match(df_link$to, gene_mat),
#'     cc = col_ver[d]
#'   )
#'   col_fun1 = circlize::colorRamp2(c(-3, 0, 3),
#'                                   c("blue", "white", "red"))
#'   col_fun2 = circlize::colorRamp2(c(-3, 0, 3),
#'                                   c("green", "white", "magenta"))
#'   ##############################################
#'   p1 <- function() {
#'     par(mar = c(0,0,1.5,0), cex=1)
#'     circlize::circos.clear()
#'     circlize::circos.par(start.degree = starting_angle,
#'                          gap.after = 1, #c(rep(2, length(cl_name)-1), 15, rep(2, length(path_name)-1), 15),
#'                          points.overflow.warning = FALSE
#'     )
#'     circlize::circos.heatmap(mat1, split = split, col = col_fun1,
#'                              rownames.side = "outside",
#'                              #track.height = 0.3,
#'                              rownames.col = col_ver[gene_col],
#'                              rownames.cex = axis_size,
#'                              bg.border = "gray", bg.lwd = 2, bg.lty = 1)
#'     circlize::circos.track(track.index = circlize::get.current.track.index(),
#'                            panel.fun = function(x, y) {
#'                              circlize::circos.text(CELL_META$xcenter,
#'                                                    CELL_META$cell.ylim[2] + circlize::convert_y(16, "mm"),
#'                                                    CELL_META$sector.index,
#'                                                    facing = "bending.inside", cex = lab_size, face = 2,
#'                                                    adj = c(0.5, 0), niceFacing = FALSE)
#'                            }, bg.border = NA)
#'     circlize::circos.heatmap(mat2, split = split, col = col_fun2,
#'                              #track.height = 0.3,
#'                              bg.border = "gray", bg.lwd = 2, bg.lty = 1)
#'     for(i in seq_len(nrow(idx_link))) {
#'       circlize::circos.heatmap.link(idx_link$from[i],
#'                                     idx_link$to[i],
#'                                     col = idx_link$cc[i]
#'       )
#'     }
#'
#'     do.call(title_tag, labs_option)
#'
#'     leg_out <- ComplexHeatmap::Legend(at = c(floor(min(mat1)), 0, ceiling(max(mat1))),
#'                                       col_fun = col_fun1,
#'                                       title_position = "leftcenter-rot", title = outer_title)
#'     leg_in <- ComplexHeatmap::Legend(at = c(floor(min(mat2)), 0, ceiling(max(mat2))),
#'                                      col_fun = col_fun2,
#'                                      title_position = "leftcenter-rot", title = inner_title)
#'     # ComplexHeatmap::draw(leg_out,
#'     #                      x = unit(1, "npc")- unit(1, "mm"),
#'     #                      y = unit(1, "npc")- unit(12, "mm"),
#'     #                      just = c("right", "bottom"))
#'     # ComplexHeatmap::draw(leg_out,
#'     #                      x = unit(1, "npc")- unit(1, "mm"),
#'     #                      y = unit(1, "npc")- unit(12, "mm"),
#'     #                      just = c("right", "bottom"))
#'     circlize::circos.clear()
#'   }
#'   grob_obj <- cowplot::as_grob(p1)
#'
#'   leg_out <- ComplexHeatmap::Legend(at = c(floor(min(mat1)), 0, ceiling(max(mat1))),
#'                                     col_fun = col_fun1,
#'                                     title_position = "leftcenter-rot", title = outer_title)
#'   leg_in <- ComplexHeatmap::Legend(at = c(floor(min(mat2)), 0, ceiling(max(mat2))),
#'                                    col_fun = col_fun2,
#'                                    title_position = "leftcenter-rot", title = inner_title)
#'   p2 <- cowplot::plot_grid(grob_obj,
#'                            cowplot::plot_grid(
#'                              grid::grid.grabExpr(ComplexHeatmap::draw(leg_out)),
#'                              grid::grid.grabExpr(ComplexHeatmap::draw(leg_in)),
#'                              nrow = 2),
#'                            ncol=2,rel_widths =c(9.3,0.7))
#'   p <- cowplot::ggdraw(p2)
#'   return(p)
#' }
#'
