library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)


classnodes <-
  create_node_df(
    n = 7,
    label = c("data.frame\\l",
              "mbnma.network\\l - plot()\\l - summary()\\l",
              "mbnma\\l - plot()\\l - summary()\\l",
              "mbnma.predict\\l - plot()\\l - summary()\\l",
              "nma\\l - plot()\\l",
              "nodesplit\\l - plot()\\l - summary()\\l",
              "mbnma.rank\\l - plot()\\l - summary()\\l"),
    #type = "lower",
    #style = "filled",
    #color = "aqua",
    color = "black",
    fontname="Consolas",
    shape = "rectangle",
    fillcolor = "Honeydew",
    fontcolor = "black",
    fontalign="left",
    height=0.7,
    width=2
    #data = c(3.5, 2.6, 9.4, 2.7)
    )

classedges <-
  create_edge_df(
    from = c("1", "2", "2", "2", "3", "3"),
    to =   c("2", "3", "5", "6", "4", "7"),
    color="black",
    #label=c("mbnma.network()", "mbnma.run()", "nma.run()", "mbnma.nodesplit()", "predict()", "rank()"),
    rel = "a",
    fontname="Consolas"
    )

g <- create_graph(nodes_df = classnodes,
                  #edges_df = classedges,
                  attr_theme = "tb")

funnodes <-
  create_node_df(
    n=9,
    label=c("mbnma.network()", "mbnma.run()", "nma.run()", "mbnma.nodesplit()", "predict()", "rank()",
            "fitplot()", "devplot()", "cumrank()"),
    shape="rectangle",
    fontname="Consolas",
    fillcolor = "white",
    fontcolor = "black",
    color="black",
    width=1.5,
    style="dashed"
  )
g <- add_node_df(g, funnodes)



funedges <-
  create_edge_df(
    from = c("1", "8", "2", "9", "2", "10", "2", "11", "3", "12", "3", "13", "4", "3", "3", "7"),
    to =   c("8", "2", "9", "3", "10", "5", "11", "6", "12", "4", "13", "7", "13", "14", "15", "16"),
    color="black",
    #label=c("mbnma.network()", "mbnma.run()", "nma.run()", "mbnma.nodesplit()", "predict()", "rank()"),
    rel = "a",
    fontname="Consolas",
    arrowhead=c(rep(c("None", "normal"),6), rep("none",4))
  )
g <- add_edge_df(g, funedges)

g <- g %>% #set_node_attrs(node_attr = width, values=2) %>%
  add_global_graph_attrs( attr = "splines",
                          value = "ortho",
                          attr_type = "graph")


render_graph(g)

render_graph(g) %>% export_svg %>% charToRaw %>%
  rsvg_png("~/MBNMA/MBNMA R Package/Dose/MBNMAdose/man/figures/functionstructure.png")




# #####################################################################################
#
#         EDGE LAYOUT
#
# #####################################################################################
#
#
# classnodes <-
#   create_node_df(
#     n = 8,
#     label = c("data.frame\\l",
#               "mbnma.network\\l - plot()\\l - summary()\\l",
#               "mbnma\\l - plot()\\l - summary()\\l",
#               "mbnma.predict\\l - plot()\\l - summary()\\l",
#               "nma\\l - plot()\\l",
#               "nodesplit\\l - plot()\\l - summary()\\l",
#               "mbnma.rank\\l - plot()\\l - summary()\\l",
#               "Plots a figure\\l"),
#     #type = "lower",
#     #style = "filled",
#     #color = "aqua",
#     fontname="Consolas",
#     shape = "none",
#     fillcolor = c("LemonChiffon", rep("Honeydew",6), "LemonChiffon"),
#     fontcolor = "black",
#     fontalign="left"
#     #data = c(3.5, 2.6, 9.4, 2.7)
#   )
#
# classedges <-
#   create_edge_df(
#     from = c("1", "2", "2", "2", "3", "3", "4", "5", "5", "7"),
#     to =   c("2", "3", "5", "6", "4", "7", "7", "8", "8", "8"),
#     color="black",
#     label=c("mbnma.network()", "mbnma.run()", "nma.run()", "mbnma.nodesplit()", "predict()", "rank()", "rank()",
#             "fitplot()", "devplot()", "cumrank()"),
#     rel = "a",
#     fontname="Consolas"
#     #splines=ortho
#
#   )
#
# g <- create_graph(nodes_df = classnodes,
#                   edges_df = classedges,
#                   attr_theme = "tb")
#
# g <- g %>% set_node_attrs(node_attr = width, values=2) %>%
#   add_global_graph_attrs( attr = "splines",
#                           value = "ortho",
#                           attr_type = "graph")
#
# render_graph(g)
#
