#library(DiagrammeR)
#library(DiagrammeRsvg)
#library(rsvg)


# Generate nodes for MBNMA classes
classnodes <-
  DiagrammeR::create_node_df(
    n = 8,
    label = c("data.frame\\l",
              "mbnma.network\\l - plot()\\l - summary()\\l",
              "mbnma\\l - plot()\\l - summary()\\l",
              "mbnma.predict\\l - plot()\\l - summary()\\l",
              "nma\\l - plot()\\l",
              "nodesplit\\l - plot()\\l - summary()\\l",
              "mbnma.rank\\l - plot()\\l - summary()\\l",
              "relative.array\\l"),
    color = "black",
    fontname="Consolas",
    shape = "rectangle",
    fillcolor = "Honeydew",
    fontcolor = "black",
    fontalign="left",
    height=0.7,
    width=2
  )

# Generate edges between classes
# classedges <-
#   create_edge_df(
#     from = c("1", "2", "2", "2", "3", "3"),
#     to =   c("2", "3", "5", "6", "4", "7"),
#     color="black",
#     rel = "a",
#     fontname="Consolas"
#     )

g <- DiagrammeR::create_graph(nodes_df = classnodes, attr_theme = "tb")

# Generate function nodes
funnodes <-
  DiagrammeR::create_node_df(
    n=11,
    label=c("mbnma.network()", "mbnma.run()", "nma.run()", "mbnma.nodesplit()", "predict()", "rank()",
            "fitplot()", "devplot()", "cumrank()", "devdev()", "get.relative()"),
    shape="rectangle",
    fontname="Consolas",
    fillcolor = "white",
    fontcolor = "black",
    color="black",
    width=1.5,
    style="dashed"
  )
g <- DiagrammeR::add_node_df(g, funnodes)


# Generate edges between classes and functions
funedges <-
  DiagrammeR::create_edge_df(
    from = c("1", "9", "2", "10", "2", "11", "2", "12", "3", "13", "3", "14", "4", "3", "3", "7", "3", "5", "19", "3", "8"),
    to =   c("9", "2", "10", "3", "11", "5", "12", "6", "13", "4", "14", "7", "14", "15", "16", "17", "18", "18", "8", "19", "14"),
    color="black",
    rel = "a",
    fontname="Consolas",
    arrowhead=c(rep(c("None", "normal"),6), rep("none",4), rep("none",2), "normal", "none", "normal")
  )
g <- DiagrammeR::add_edge_df(g, funedges)


# Set graph attributes
g <- g %>%
  DiagrammeR::add_global_graph_attrs( attr = "splines",
                          value = "ortho",
                          attr_type = "graph")

# Render graph
DiagrammeR::render_graph(g)

# Save graph
DiagrammeR::render_graph(g) %>% DiagrammeRsvg::export_svg %>% base::charToRaw %>%
  rsvg::rsvg_png("~/MBNMA/MBNMA R Package/Dose/MBNMAdose/man/figures/functionstructure.png")




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
