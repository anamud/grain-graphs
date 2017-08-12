#
# Setup
#

# Clean slate
rm(list=ls())

# Include support functions
grain_graphs_root_dir <- Sys.getenv("GRAIN_GRAPHS_ROOT")
source(paste(grain_graphs_root_dir,"/prototype/common.R", sep=""))

# Default configs
config_dir <- paste(grain_graphs_root_dir, "/prototype/configs", sep="")

# Read arguments
Rstudio_mode <- F
if (Rstudio_mode) {
  cl_args <- list(graph="grain-graph.graphml",
                 data="aggregated-grain-graph-data.csv",
                 groupname="fj1",
                 enumcriticalpath=F,
                 grainpropertyconfig="grain-properties.cfg",
                 edgepropertyconfig="edge-properties.cfg",
                 verbose=T,
                 timing=F,
                 out="subgraph")
} else {
  option_list <- list(
    make_option(c("--graph"), default="grain-graph.graphml", help="Non-aggregated grain graph [default \"%default\"].", metavar="FILE"), # Cannot add -g as option since it is parsed by Rscript as "gui" option
    make_option(c("--data"), default="aggregated-grain-graph-data.csv", help="Aggregated grain graph data [default \"%default\"].", metavar="FILE"),
    make_option(c("--groupname"), default="fj1", help="Group at root of subgraph [default \"%default\"].", metavar="INT"),
    make_option(c("--enumcriticalpath"), action="store_true", default=FALSE, help="Enumerate nodes on subgraph critical path."),
    make_option(c("--grainpropertyconfig"), default=paste(config_dir, "/grain-properties.cfg", sep=""), help = "Grain property configuration file [default \"%default\"].", metavar="FILE"),
    make_option(c("--edgepropertyconfig"), default=paste(config_dir, "/edge-properties.cfg", sep=""), help = "Edge property configuration file [default \"%default\"].", metavar="FILE"),
    make_option(c("--verbose"), action="store_true", default=TRUE, help="Print output [default]."),
    make_option(c("--timing"), action="store_true", default=FALSE, help="Print timing information."),
    make_option(c("--out"), default="subgraph", help="Subgraph output file prefix [default \"%default\"]", metavar="STRING"),
    make_option(c("--quiet"), action="store_false", dest="verbose", help="Print little output."))

  cl_args <- parse_args(OptionParser(option_list = option_list), args = commandArgs(TRUE))

  if (!exists("groupname", where=cl_args)) {
    my_print("Error: Group name argument missing. Check help (-h).")
    quit("no", 1)
  }
}

if (cl_args$verbose) my_print("Reading inputs ...")
# Read graph
g_graphml <- read.graph(cl_args$graph, format="graphml")

# Read data
g_data <- read.csv(cl_args$data, header=FALSE, sep=';', comment.char='#', na.strings="NA")
colnames(g_data) <- c("group", "member", "grouptype")
g_data[g_data == "NA"] <- NA
is.na(g_data) <- is.na(g_data)

if (cl_args$verbose) my_print("Initializing ...")
if (cl_args$timing) tic(type="elapsed")

# Property query functions
get_value <- function(prop_cfg, type_type, property_property)
{
    value1 <- subset(prop_cfg, type == type_type & property == property_property, select = value1)
    if (nrow(value1) != 1) {
        my_print(paste("Error: Ambiguous values for property", property_property, "for type", type_type, "!"))
        quit("no", 1)
    }
    value1 <- as.character(unlist(value1))
    value2 <- as.character(unlist(subset(prop_cfg, type == type_type & property == property_property, select = value2)))
    value1_first <- substr(value1, 1, 1)
    value1_last <- substr(value1, nchar(value1), nchar(value1))
    if (value1_first == '[' && value1_last == ']') {
        ret_val <- c(substr(value1, 2, nchar(value1) - 1), value2)
    } else {
        ret_val <- c(value1, NA)
    }

    return(ret_val)
}

# Read property configuration files
#grain_prop_cfg <- read.csv(cl_args$grainpropertyconfig, header=TRUE, comment.char='#', na.strings="NA")
edge_prop_cfg <- read.csv(cl_args$edgepropertyconfig, header=TRUE, comment.char='#', na.strings="NA")

# Prepare information output file
g_subgraph_info_out_file <- paste(gsub(". $", "", cl_args$out), paste("-",cl_args$groupname, ".info", sep=""), sep="")
sink(g_subgraph_info_out_file)
sink()

# Function to return members of group
get_members <- function(group_name)
{
    #my_print(paste("Processing group", group_name, "..."))

    temp <- unname(unlist(g_data %>% filter(group == group_name) %>% select(member)))
    members <- unlist(strsplit(temp, ','))

    stopifnot(length(members) > 0)

    ungrouped_members <- c()

    for (m in members) {
        if (grepl("^s", m) | grepl("^fj", m) | grepl("^l", m)) {
            temp <- get_members(m)
            ungrouped_members <- append(ungrouped_members, temp)
        } else {
            ungrouped_members <- append(ungrouped_members, m)
        }
    }

    ungrouped_members <- unique(ungrouped_members)

    return (ungrouped_members)
}

if (cl_args$timing) toc("Initializing")

if (cl_args$verbose) my_print(paste("Getting subgraph enclosed by group", cl_args$groupname, "..."))
if (cl_args$timing) tic(type="elapsed")

# Get members of group provided
# ... group exists in data
# ... else abort
members <- get_members(cl_args$groupname)

g_subgraph <- induced_subgraph(g_graphml, members, impl="auto")

if (cl_args$timing) toc("Getting subgraph")

if (cl_args$verbose) my_print("Computing metrics ...")
if (cl_args$timing) tic(type="elapsed")

if (length(V(g_subgraph)) != 1) {
    # Compute metrics
    if (!cl_args$enumcriticalpath) {
        # Compute critical path length using Bellman Ford algorithm with negative weights
        start_index = which(degree(g_subgraph, mode = "in") == 0)
        stopifnot(length(start_index) == 1)
        end_index = which(degree(g_subgraph, mode = "out") == 0)
        stopifnot(length(end_index) == 1)
        shortest_path <- shortest.paths(g_subgraph, v=start_index, to=end_index, mode="out")
        critical_path <- -as.numeric(shortest_path)
    } else {
        num_vertices <- length(V(g_subgraph))
        if (cl_args$verbose) {
            progress_bar <- txtProgressBar(min = 0, max = num_vertices, style = 3)
            ctr <- 0
        }
        # Topological sort
        top_sort_graph <- topological.sort(g_subgraph)
        # Set root path attributes
        V(g_subgraph)[top_sort_graph[1]]$root_dist <- 0
        V(g_subgraph)[top_sort_graph[1]]$depth <- 0
        V(g_subgraph)[top_sort_graph[1]]$root_path <- top_sort_graph[1]
        # Get data frame of g_subgraph object
        g_subgraph_vertices <- get.data.frame(g_subgraph, what="vertices")
        # Get longest paths from root
        for(node in top_sort_graph[-1])
        {
            # Get distance from node's predecessors
            incident_edges <- incident(g_subgraph, node, mode="in")
            incident_edge_weights <- -E(g_subgraph)[incident_edges]$weight
            # Get distance from root to node's predecessors
            adjacent_nodes <- neighbors(g_subgraph, node, mode="in")
            adjacent_nodes_root_dist <- g_subgraph_vertices$root_dist[adjacent_nodes]
            # Add distances (assuming one-one corr.)
            root_dist <- incident_edge_weights + adjacent_nodes_root_dist
            # Set node's distance from root to max of added distances
            max_root_dist <- max(root_dist)
            g_subgraph_vertices$root_dist[node] <- max_root_dist
            # Set node's path from root to path of max of added distances
            nodes_on_root_path <- as.vector(adjacent_nodes)[match(max_root_dist,root_dist)]
            root_path <- list(c(unlist(g_subgraph_vertices$root_path[nodes_on_root_path]),node))
            g_subgraph_vertices$root_path[node] <- root_path
            # Set node's depth as one greater than the largest depth its predecessors
            g_subgraph_vertices$depth[node] <- max(g_subgraph_vertices$depth[adjacent_nodes]) + 1
            if (cl_args$verbose) {
                ctr <- ctr + 1
                setTxtProgressBar(progress_bar, ctr)
            }
        }
        ## Critical path is the one with the largest root distance
        critical_path <- max(g_subgraph_vertices$root_dist)
        # Enumerate nodes on critical path
        critical_nodes <- unlist(g_subgraph_vertices$root_path[match(critical_path,g_subgraph_vertices$root_dist)])
        g_subgraph_vertices$on_crit_path <- 0
        g_subgraph_vertices$on_crit_path[critical_nodes] <- 1
        # Mark critical nodes and edges on grain graph
        g_subgraph <- set.vertex.attribute(g_subgraph, name="subgraph_on_crit_path", index=V(g_subgraph), value=g_subgraph_vertices$on_crit_path)
        g_subgraph <- set.vertex.attribute(g_subgraph, name="subgraph_root_dist", index=V(g_subgraph), value=g_subgraph_vertices$root_dist)
        g_subgraph <- set.vertex.attribute(g_subgraph, name="subgraph_depth", index=V(g_subgraph), value=g_subgraph_vertices$depth)
        critical_edges <- E(g_subgraph)[V(g_subgraph)[subgraph_on_crit_path==1] %--% V(g_subgraph)[subgraph_on_crit_path==1]]
        g_subgraph <- set.edge.attribute(g_subgraph, name="subgraph_on_crit_path", index=critical_edges, value=1)
        g_subgraph <- set.edge.attribute(g_subgraph, name="color", index=critical_edges, value="#0000FF")
        #g_subgraph <- set.vertex.attribute(g_subgraph, name="border-color", index=critical_nodes, value="#0000FF")
        # Cleanup
        g_subgraph <- remove.vertex.attribute(g_subgraph,"root_path")
        if (cl_args$verbose) {
            ctr <- ctr + 1
            setTxtProgressBar(progress_bar, ctr)
            close(progress_bar)
        }
    }

    #
    # Compute basic information about grain graph
    #

    # Get subgraph data restricted to tasks
    g_subgraph_vertices <- get.data.frame(g_subgraph, what="vertices")
    g_subgraph_data <- g_subgraph_vertices
    g_subgraph_data[g_subgraph_data == "NA"] <- NA
    is.na(g_subgraph_data) <- is.na(g_subgraph_data)
    g_subgraph_data <- subset(g_subgraph_data, !is.na(task))

    # Set edge weights
    common_edge_weight <- get_value(edge_prop_cfg, "common", "weight")
    if (!is.na(common_edge_weight[2])) {
        if (!(common_edge_weight[1] %in% colnames(g_subgraph_data))) {
            my_print(paste("Error: Mapped variable", common_edge_weight[1], "not found in profiling data!"))
            quit("no", 1)
        }
    }

    # Compute simple metrics
    fork_nodes_index <- startsWith(V(g_subgraph)$name, 'f')
    sink(g_subgraph_info_out_file, append=T)
    my_print("# Structure")
    my_print(paste("Number of nodes =", length(V(g_subgraph))))
    my_print(paste("Number of edges =", length(E(g_subgraph))))
    my_print(paste("Number of tasks =", length(g_subgraph_data$task)))
    my_print(paste("Number of forks =", length(unique(as.character(get.vertex.attribute(g_subgraph, name="name", index=fork_nodes_index))))))
    my_print("# Out-degree distribution of forks")
    degree.distribution(g_subgraph, v=fork_nodes_index, mode="out")
    if (!is.na(common_edge_weight[2])) {
        my_print(paste("# Cilk theory parallelism (metric = ", common_edge_weight[1], ")", sep=""))
        my_print(paste("Span (critical path) =", critical_path))
        work <- sum(as.numeric(g_subgraph_data[,common_edge_weight[1]]))
        my_print(paste("Work =", work))
        my_print(paste("Parallelism (Work/Span) =", work/critical_path))
    } else {
        my_print(paste("# Cilk theory parallelism (metric = constant", common_edge_weight[1], ")", sep=""))
        my_print(paste("Span (critical path) =", critical_path))
        work <- nrow(g_subgraph_data)*common_edge_weight[1]
        my_print(paste("Work =", work))
        my_print(paste("Parallelism (Work/Span) =", work/critical_path))
    }
    if (cl_args$enumcriticalpath) {
        my_print(paste("Number of critical tasks =", length(g_subgraph_vertices$name[g_subgraph_vertices$subgraph_on_crit_path == 1 & g_subgraph_vertices$type == "task"])))
    }
    sink()
    my_print(paste("Wrote file:", g_subgraph_info_out_file))
} else {
    #
    # Compute basic information about grain graph
    #
    sink(g_subgraph_info_out_file, append=T)
    my_print("# Structure")
    my_print(paste("Number of nodes =", length(V(g_subgraph))))
    sink()
    my_print(paste("Wrote file:", g_subgraph_info_out_file))
}

if (cl_args$timing) toc("Computing metrics")

if (cl_args$verbose) my_print("Plotting subgraph ...")
if (cl_args$timing) tic(type="elapsed")

# Plot subgraph
temp_out_file <- paste(gsub(". $", "", cl_args$out), paste("-",cl_args$groupname, ".graphml", sep=""), sep="")
res <- write.graph(g_subgraph, file=temp_out_file, format="graphml")
my_print(paste("Wrote file:", temp_out_file))

# Plot subgraph without attributes
temp_out_file <- paste(gsub(". $", "", cl_args$out), paste("-",cl_args$groupname, "-noattrib.graphml", sep=""), sep="")
g_subgraph_noattrib <- g_subgraph
attribs <- vertex_attr_names(g_subgraph_noattrib)
for (attrib in attribs) {
    g_subgraph_noattrib <- delete_vertex_attr(g_subgraph_noattrib, attrib)
}
attribs <- edge_attr_names(g_subgraph_noattrib)
for (attrib in attribs) {
    g_subgraph_noattrib <- delete_edge_attr(g_subgraph_noattrib, attrib)
}
res <- write.graph(g_subgraph_noattrib, file=temp_out_file, format="graphml")
my_print(paste("Wrote file:", temp_out_file))

if (cl_args$timing) toc("Plotting subgraph")

# Warn
w <- warnings()
if (class(w) != "NULL")
  my_print(w)
