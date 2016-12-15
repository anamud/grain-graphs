# Clean slate
rm(list=ls())

# Include
mir_root <- Sys.getenv("GRAIN_GRAPHS_ROOT")
source(paste(mir_root,"/prototype/common.R",sep=""))

# Parse args
Rstudio_mode <- F
if (Rstudio_mode) {
    cl_args <- list(data="task-stats.processed",
                   grainpropertyconfig="grain-properties.cfg",
                   edgepropertyconfig="edge-properties.cfg",
                   out="grain-graph",
                   enumcriticalpath=F,
                   forloop=F,
                   layout=F,
                   verbose=T,
                   timing=F)
} else {
    option_list <- list(make_option(c("-d","--data"), help = "Task profiling data.", metavar="FILE"),
                        make_option(c("--grainpropertyconfig"), default="grain-properties.cfg", help = "Grain property configuration file [default \"%default\"].", metavar="FILE"),
                        make_option(c("--edgepropertyconfig"), default="edge-properties.cfg", help = "Edge property configuration file [default \"%default\"].", metavar="FILE"),
                        make_option(c("-o","--out"), default="grain-graph", help = "Output file suffix [default \"%default\"].", metavar="STRING"),
                        make_option(c("--enumcriticalpath"), action="store_true", default=FALSE, help="Enumerate nodes critical path."),
                        make_option(c("--forloop"), action="store_true", default=FALSE, help="Task profiling data obtained from a for-loop program."),
                        make_option(c("--layout"), action="store_true", default=FALSE, help="Layout using Sugiyama style and plot to PDF."),
                        make_option(c("--verbose"), action="store_true", default=TRUE, help="Print output [default]."),
                        make_option(c("--quiet"), action="store_false", dest="verbose", help="Print little output."),
                        make_option(c("--timing"), action="store_true", default=FALSE, help="Print processing time."))

    cl_args <- parse_args(OptionParser(option_list = option_list), args = commandArgs(TRUE))

    if (!exists("data", where=cl_args)) {
        my_print("Error: Invalid arguments. Check help (-h)")
        quit("no", 1)
    }
}

if (cl_args$verbose) my_print("Initializing ...")

# Set grain and edge property configuration
grain_prop_cfg <- read.csv(cl_args$grainpropertyconfig, header=TRUE)
edge_prop_cfg <- read.csv(cl_args$edgepropertyconfig, header=TRUE)

# Grain sizes
fork_size <- as.numeric(unlist(subset(grain_prop_cfg, type == "fork" & property == "size", select = value1)))
join_size <- as.numeric(unlist(subset(grain_prop_cfg, type == "join" & property == "size", select = value1)))
start_size <- as.numeric(unlist(subset(grain_prop_cfg, type == "start" & property == "size", select = value1)))
end_size <- as.numeric(unlist(subset(grain_prop_cfg, type == "end" & property == "size", select = value1)))
task_size <- as.numeric(unlist(subset(grain_prop_cfg, type == "fragment" & property == "size", select = value1)))
task_size_mult <- as.numeric(unlist(subset(grain_prop_cfg, type == "fragment" & property == "mult", select = value1)))
task_size_bins <- as.numeric(unlist(subset(grain_prop_cfg, type == "fragment" & property == "bins", select = value1)))

# Grain shapes
fork_shape <- as.character(unlist(subset(grain_prop_cfg, type == "fork" & property == "shape", select = value1)))
join_shape <- as.character(unlist(subset(grain_prop_cfg, type == "join" & property == "shape", select = value1)))
start_shape <- as.character(unlist(subset(grain_prop_cfg, type == "start" & property == "shape", select = value1)))
end_shape <- as.character(unlist(subset(grain_prop_cfg, type == "end" & property == "shape", select = value1)))
task_shape <- as.character(unlist(subset(grain_prop_cfg, type == "fragment" & property == "shape", select = value1)))

# Grain colors
fork_color <- as.character(unlist(subset(grain_prop_cfg, type == "fork" & property == "color", select = value1)))
join_color <- as.character(unlist(subset(grain_prop_cfg, type == "join" & property == "color", select = value1)))
start_color <- as.character(unlist(subset(grain_prop_cfg, type == "start" & property == "color", select = value1)))
end_color <- as.character(unlist(subset(grain_prop_cfg, type == "end" & property == "color", select = value1)))
task_color <- as.character(unlist(subset(grain_prop_cfg, type == "fragment" & property == "color", select = value1)))
color_fun <- colorRampPalette(as.character(unlist(subset(grain_prop_cfg, type == "fragment" & property == "color_gradient", select = c(value1,value2)))))
task_color_bins <- as.numeric(unlist(subset(grain_prop_cfg, type == "fragment" & property == "bins", select = value1)))
task_color_pal <- color_fun(task_color_bins)

# Edge properties
create_edge_color <- as.character(unlist(subset(edge_prop_cfg, type == "create" & property == "color", select = value1)))
sync_edge_color <- as.character(unlist(subset(edge_prop_cfg, type == "sync" & property == "color", select = value1)))
scope_edge_color <- as.character(unlist(subset(edge_prop_cfg, type == "scope" & property == "color", select = value1)))
cont_edge_color <- as.character(unlist(subset(edge_prop_cfg, type == "continuation" & property == "color", select = value1)))
temp <- as.character(unlist(subset(edge_prop_cfg, type == "common" & property == "weight", select = value1)))
path_weight <- substr(temp, 2, nchar(temp) - 1)

# Read data
prof_data <- read.csv(cl_args$data, header=TRUE)

# Information output
grain_graph_info_out_file <- paste(gsub(". $", "", cl_args$out), ".info", sep="")
sink(grain_graph_info_out_file)
sink()

# Remove background task
prof_data <- prof_data[!is.na(prof_data$parent),]

# Path weight avaiability
if (!(path_weight %in% colnames(prof_data))) {
    my_print(paste("Path weight variable", path_weight, "not found in profiling data!"))
    path_weight <- NA
}

if (cl_args$forloop) {
    # Remove idle task without children
    prof_data <- prof_data[!(prof_data$tag == "idle_task" & prof_data$num_children == 0),]
}

if (cl_args$verbose) my_print("Creating grain graph ...")

# Create node lists
if (cl_args$timing) tic(type="elapsed")

# Create join nodes list
join_nodes <- mapply(function(x, y, z) {paste('j', x, y, sep='.')}, x=prof_data$parent, y=prof_data$joins_at)
join_nodes_unique <- unique(unlist(join_nodes, use.names=FALSE))

# Create parent nodes list
parent_nodes_unique <- unique(prof_data$parent)

# Create fork nodes list
fork_nodes <- mapply(function(x, y, z) {paste('f', x, y, sep='.')}, x=prof_data$parent, y=prof_data$joins_at)
fork_nodes_unique <- unique(unlist(fork_nodes, use.names=FALSE))

if (cl_args$timing) toc("Node list creation")

# Create grain graph
if (cl_args$timing) tic(type="elapsed")

grain_graph <- graph.empty(directed=TRUE) + vertices('E',
                                            unique(c(join_nodes_unique,
                                                     fork_nodes_unique,
                                                     parent_nodes_unique,
                                                     prof_data$task)))

if (cl_args$timing) toc("grain graph creation")

# Connect parent fork to task
if (cl_args$verbose) my_print("Connecting nodes ...")
if (cl_args$timing) tic(type="elapsed")

grain_graph[from=fork_nodes, to=prof_data$task, attr='type'] <- 'create'
grain_graph[from=fork_nodes, to=prof_data$task, attr='color'] <- create_edge_color

if (cl_args$timing) toc("Connect parent fork to task")

# Connect parent task to first fork
if (cl_args$timing) tic(type="elapsed")

first_forks_index <- which(grepl("f.[0-9]+.0$", fork_nodes_unique))
parent_first_forks <- as.vector(sapply(fork_nodes_unique[first_forks_index], function(x) {gsub('f.(.*)\\.+.*','\\1', x)}))
first_forks <- fork_nodes_unique[first_forks_index]

grain_graph[to=first_forks, from=parent_first_forks, attr='type'] <- 'scope'
grain_graph[to=first_forks, from=parent_first_forks, attr='color'] <- scope_edge_color

if (!is.na(path_weight)) {
    temp <- as.numeric(prof_data[match(parent_first_forks, prof_data$task),path_weight])

    grain_graph[to=first_forks, from=parent_first_forks, attr=path_weight] <- temp
    grain_graph[to=first_forks, from=parent_first_forks, attr='weight'] <- -temp
}

if (cl_args$timing) toc("Connect parent to first fork")

# Connect leaf task to join
if (cl_args$timing) tic(type="elapsed")

leaf_tasks <- prof_data$task[prof_data$leaf == T]
leaf_join_nodes <- join_nodes[match(leaf_tasks, prof_data$task)]

grain_graph[from=leaf_tasks, to=leaf_join_nodes, attr='type'] <- 'sync'
grain_graph[from=leaf_tasks, to=leaf_join_nodes, attr='color'] <- sync_edge_color

if (!is.na(path_weight)) {
    temp <- as.numeric(prof_data[match(leaf_tasks, prof_data$task),path_weight])
    grain_graph[from=leaf_tasks, to=leaf_join_nodes, attr=path_weight] <- temp
    grain_graph[from=leaf_tasks, to=leaf_join_nodes, attr='weight'] <- -temp
}

if (cl_args$timing) toc("Connect leaf task to join")

# Connect join to next fork
if (cl_args$timing) tic(type="elapsed")

#Rprof("profile-jointonext.out")
find_next_fork <- function(node)
{
    #my_print(paste('Processing node',node, sep=" "))

    # Get node info
    node_split <- unlist(strsplit(node, "\\."))
    parent <- as.numeric(node_split[2])
    join_count <- as.numeric(node_split[3])

    # Find next fork
    next_fork <- paste('f', as.character(parent), as.character(join_count+1), sep=".")
    if (is.na(match(next_fork, fork_nodes_unique)) == F) {
        next_fork <- next_fork # Connect to next fork
    } else {

        # Next fork is part of grandfather
        parent_index <- match(parent, prof_data$task)
        gfather <- prof_data[parent_index,]$parent
        gfather_join <- paste('j', as.character(gfather), as.character(prof_data[parent_index,]$joins_at), sep=".")

        if (is.na(match(gfather_join, join_nodes_unique)) == F) {
            next_fork <- gfather_join # Connect to grandfather's join
        } else {
            next_fork <- 'E' # Connect to end node
        }
    }
    next_fork
}

next_forks <- as.vector(sapply(join_nodes_unique, find_next_fork))

grain_graph[from=join_nodes_unique, to=next_forks, attr='type'] <- 'continue'
grain_graph[from=join_nodes_unique, to=next_forks, attr='color'] <- cont_edge_color

#Rprof(NULL)
if (cl_args$timing) toc("Connect join to next fork")

# Set attributes
if (cl_args$verbose) my_print("Setting attributes ...")
if (cl_args$timing) tic(type="elapsed")

# Common vertex attributes
V(grain_graph)$label <- V(grain_graph)$name

# Set task vertex attributes
task_index <- match(as.character(prof_data$task), V(grain_graph)$name)

# Set annotations
for (annot in colnames(prof_data)) {
    values <- as.character(prof_data[,annot])
    grain_graph <- set.vertex.attribute(grain_graph, name=annot, index=task_index, value=values)
}

# Set size constants
grain_graph <- set.vertex.attribute(grain_graph, name='size', index=task_index, value=task_size)
grain_graph <- set.vertex.attribute(grain_graph, name='width', index=task_index, value=task_size)
grain_graph <- set.vertex.attribute(grain_graph, name='height', index=task_index, value=task_size)
grain_graph <- set.vertex.attribute(grain_graph, name='shape', index=task_index, value=task_shape)

# Scale size to attributes
size_scaled <- c("ins_count", "work_cycles", "overhead_cycles", "exec_cycles")
for(attrib in size_scaled) {
    if (attrib %in% colnames(prof_data)) {

        # Set size
        attrib_unique <- unique(prof_data[,attrib])
        if (length(attrib_unique) == 1) {
            p_task_size <- task_size_mult
        } else if(length(attrib_unique) == 2 & identical(c(0,NA), as.numeric(attrib_unique[order(attrib_unique)]))) {
            p_task_size <- task_size_mult
        } else {
            p_task_size <- task_size_mult * as.numeric(cut(prof_data[,attrib], task_size_bins))
        }
        annot_name <- paste(attrib, "_to_size", sep="")
        grain_graph <- set.vertex.attribute(grain_graph, name=annot_name, index=task_index, value=p_task_size)

        # Set height
        attrib_val <- prof_data[,attrib]
        attrib_val_norm <- 1 + ((attrib_val - min(attrib_val)) / (max(attrib_val) - min(attrib_val)))
        annot_name <- paste(attrib, "_to_height", sep="")
        grain_graph <- set.vertex.attribute(grain_graph, name=annot_name, index=task_index, value=attrib_val_norm*task_size)
    }
}

# Set color constants
grain_graph <- set.vertex.attribute(grain_graph, name='color', index=task_index, value=task_color)

# Scale color to attributes
# "-" in attribute name implies higher is better
attrib_color_scaled <- c("mem_fp", "-compute_int", "PAPI_RES_STL_sum", "-mem_hier_util", "work_deviation", "overhead_deviation", "-parallel_benefit", "-min_shape_contrib", "-max_shape_contrib","-median_shape_contrib", "sibling_work_balance")
if (!cl_args$forloop) {
    attrib_color_scaled <- c(attrib_color_scaled, c("sibling_scatter"))
} else {
    attrib_color_scaled <- c(attrib_color_scaled, c("chunk_work_balance", "chunk_work_cpu_balance"))
}
for(attrib in attrib_color_scaled) {
    invert_colors <- F
    if (substring(attrib, 1, 1) == "-") {
        invert_colors <- T
        attrib <- substring(attrib, 2, nchar(attrib))
    }

    if (attrib %in% colnames(prof_data)) {

        # Set color in proportion to attrib
        attrib_unique <- unique(prof_data[,attrib])
        if (invert_colors) {
            if (length(attrib_unique) == 1) {
                p_task_color <- task_color_pal[task_color_bins]
            } else if(length(attrib_unique) == 2 & identical(c(0,NA), as.numeric(attrib_unique[order(attrib_unique)]))) {
                p_task_color <- task_color_pal[task_color_bins]
            } else {
                p_task_color <- rev(task_color_pal)[as.numeric(cut(prof_data[,attrib], task_color_bins))]
            }
        } else {
            if (length(attrib_unique) == 1) {
                p_task_color <- task_color_pal[1]
            } else if(length(attrib_unique) == 2 & identical(c(0,NA), as.numeric(attrib_unique[order(attrib_unique)]))) {
                p_task_color <- task_color_pal[1]
            } else {
                p_task_color <- task_color_pal[as.numeric(cut(prof_data[,attrib], task_color_bins))]
            }
        }
        annot_name <- paste(attrib, "_to_color", sep="")
        grain_graph <- set.vertex.attribute(grain_graph, name=annot_name, index=task_index, value=p_task_color)

        # Write colors for reference
        temp_out_file <- paste(gsub(". $", "", cl_args$out), annot_name, sep=".")
        if (length(attrib_unique) == 1) {
            write.csv(data.frame(value=attrib_unique, color=p_task_color), temp_out_file, row.names=F)
        } else if(length(attrib_unique) == 2 & identical(c(0,NA), as.numeric(attrib_unique[order(attrib_unique)]))) {
            write.csv(data.frame(value=attrib_unique, color=p_task_color), temp_out_file, row.names=F)
        } else {
            v <- unique(cut(prof_data[,attrib], task_color_bins))
            write.csv(data.frame(value=v, color=task_color_pal[as.numeric(v)]), temp_out_file, row.names=F)
        }
        my_print(paste("Wrote file:", temp_out_file))
    }
}

# Set attributes to distinct color
attrib_color_distinct <- c("cpu_id", "outl_func", "tag", "outline_function")
for(attrib in attrib_color_distinct) {
    if (attrib %in% colnames(prof_data)) {

        # Map distinct color to attrib
        attrib_val <- as.character(prof_data[,attrib])
        unique_attrib_val <- unique(attrib_val)
        #attrib_color <- rainbow(length(unique_attrib_val), start=0.4, end=0.8)
        attrib_color <- rev(topo.colors(length(unique_attrib_val)))
        annot_name <- paste(attrib, "_to_color", sep="")
        grain_graph <- set.vertex.attribute(grain_graph, name=annot_name, index=task_index, value=attrib_color[match(attrib_val, unique_attrib_val)])

        # Write colors for reference
        temp_out_file <- paste(gsub(". $", "", cl_args$out), annot_name, sep=".")
        write.csv(data.frame(value=unique_attrib_val, color=attrib_color), temp_out_file, row.names=F)
        my_print(paste("Wrote file:", temp_out_file))
    }
}

# Set label and color of 'task 0'
start_index <- V(grain_graph)$name == '0'
grain_graph <- set.vertex.attribute(grain_graph, name='color', index=start_index, value=start_color)
grain_graph <- set.vertex.attribute(grain_graph, name='label', index=start_index, value='S')
grain_graph <- set.vertex.attribute(grain_graph, name='size', index=start_index, value=start_size)
grain_graph <- set.vertex.attribute(grain_graph, name='shape', index=start_index, value=start_shape)

# Set label and color of 'task E'
end_index <- V(grain_graph)$name == "E"
grain_graph <- set.vertex.attribute(grain_graph, name='color', index=end_index, value=end_color)
grain_graph <- set.vertex.attribute(grain_graph, name='label', index=end_index, value='E')
grain_graph <- set.vertex.attribute(grain_graph, name='size', index=end_index, value=end_size)
grain_graph <- set.vertex.attribute(grain_graph, name='shape', index=end_index, value=end_shape)

# Set fork vertex attributes
fork_nodes_index <- startsWith(V(grain_graph)$name, 'f')
grain_graph <- set.vertex.attribute(grain_graph, name='size', index=fork_nodes_index, value=fork_size)
grain_graph <- set.vertex.attribute(grain_graph, name='color', index=fork_nodes_index, value=fork_color)
grain_graph <- set.vertex.attribute(grain_graph, name='label', index=fork_nodes_index, value='^')
grain_graph <- set.vertex.attribute(grain_graph, name='shape', index=fork_nodes_index, value=fork_shape)

# Set join vertex attributes
join_nodes_index <- startsWith(V(grain_graph)$name, 'j')
grain_graph <- set.vertex.attribute(grain_graph, name='size', index=join_nodes_index, value=join_size)
grain_graph <- set.vertex.attribute(grain_graph, name='color', index=join_nodes_index, value=join_color)
grain_graph <- set.vertex.attribute(grain_graph, name='label', index=join_nodes_index, value='*')
grain_graph <- set.vertex.attribute(grain_graph, name='shape', index=join_nodes_index, value=join_shape)

# Set edge attributes
if (!is.na(path_weight)) {
    grain_graph <- set.edge.attribute(grain_graph, name="weight", index=which(is.na(E(grain_graph)$weight)), value=0)
    grain_graph <- set.edge.attribute(grain_graph, name=path_weight, index=which(is.na(get.edge.attribute(grain_graph, name=path_weight, index=E(grain_graph)))), value=0)
}

if (cl_args$timing) toc("Attribute setting")

# Check if grain_graph has bad structure
if (cl_args$verbose) my_print("Checking for bad structure ...")
if (cl_args$timing) tic(type="elapsed")

if (is.element(0, degree(grain_graph, fork_nodes_index, mode = c("in")))) {
    my_print("Warning! One or more fork nodes have zero degree since one or more tasks in the program performed empty synchronization.")
    my_print("Aborting on error!")
    quit("no", 1)
}
if (is.element(0, degree(grain_graph, fork_nodes_index, mode = c("out")))) {
    my_print("Warning! One or more fork nodes have zero degree since one or more tasks in the program performed empty synchronization.")
    my_print("Aborting on error!")
    quit("no", 1)
}
if (is.element(0, degree(grain_graph, join_nodes_index, mode = c("in")))) {
    my_print("Warning! One or more join nodes have zero degree since one or more tasks in the program performed empty synchronization.")
    my_print("Aborting on error!")
    quit("no", 1)
}
if (is.element(0, degree(grain_graph, join_nodes_index, mode = c("out")))) {
    my_print("Warning! One or more join nodes have zero degree since one or more tasks in the program performed empty synchronization.")
    my_print("Aborting on error!")
    quit("no", 1)
}

if (cl_args$timing) toc("Checking for bad structure")

# Calculate critical path
if (!is.na(path_weight)) {
    if (cl_args$verbose) my_print("Calculating critical path ...")
    if (cl_args$timing) tic(type="elapsed")
    # Simplify - DO NOT USE. Fucks up the critical path analysis.
    #grain_graph <- simplify(grain_graph, edge.attr.comb=toString)

    # Compute critical path
    #Rprof("profile-critpathcalc.out")
    if (!cl_args$enumcriticalpath) {
        # Compute critical path length using Bellman Ford algorithm with negative weights
        shortest_path <- shortest.paths(grain_graph, v=start_index, to=end_index, mode="out")
        critical_path <- -as.numeric(shortest_path)
    } else {
        num_vertices <- length(V(grain_graph))
        if (cl_args$verbose) {
            progress_bar <- txtProgressBar(min = 0, max = num_vertices, style = 3)
            ctr <- 0
        }
        # Topological sort
        top_sort_graph <- topological.sort(grain_graph)
        # Set root path attributes
        V(grain_graph)[top_sort_graph[1]]$root_dist <- 0
        V(grain_graph)[top_sort_graph[1]]$depth <- 0
        V(grain_graph)[top_sort_graph[1]]$root_path <- top_sort_graph[1]
        # Get data frame of grain_graph object
        graph_vertices <- get.data.frame(grain_graph, what="vertices")
        # Get longest paths from root
        for(node in top_sort_graph[-1])
        {
            # Get distance from node's predecessors
            incident_edges <- incident(grain_graph, node, mode="in")
            incident_edge_weights <- -E(grain_graph)[incident_edges]$weight
            # Get distance from root to node's predecessors
            adjacent_nodes <- neighbors(grain_graph, node, mode="in")
            adjacent_nodes_root_dist <- graph_vertices$root_dist[adjacent_nodes]
            # Add distances (assuming one-one corr.)
            root_dist <- incident_edge_weights + adjacent_nodes_root_dist
            # Set node's distance from root to max of added distances
            max_root_dist <- max(root_dist)
            graph_vertices$root_dist[node] <- max_root_dist
            # Set node's path from root to path of max of added distances
            nodes_on_root_path <- as.vector(adjacent_nodes)[match(max_root_distances,root_dist)]
            root_path <- list(c(unlist(graph_vertices$root_path[nodes_on_root_path]),node))
            graph_vertices$root_path[node] <- root_path
            # Set node's depth as one greater than the largest depth its predecessors
            graph_vertices$depth[node] <- max(graph_vertices$depth[adjacent_nodes]) + 1
            if (cl_args$verbose) {
                ctr <- ctr + 1
                setTxtProgressBar(progress_bar, ctr)
            }
        }
        ## Critical path is the one with the largest root distance
        critical_path <- max(graph_vertices$root_dist)
        # Enumerate nodes on critical path
        critical_nodes <- unlist(graph_vertices$root_path[match(critical_path,graph_vertices$root_dist)])
        graph_vertices$on_crit_path <- 0
        graph_vertices$on_crit_path[critical_nodes] <- 1
        # Mark critical nodes and edges on grain graph
        grain_graph <- set.vertex.attribute(grain_graph, name="on_crit_path", index=V(grain_graph), value=graph_vertices$on_crit_path)
        grain_graph <- set.vertex.attribute(grain_graph, name="root_dist", index=V(grain_graph), value=graph_vertices$root_dist)
        grain_graph <- set.vertex.attribute(grain_graph, name="depth", index=V(grain_graph), value=graph_vertices$depth)
        critical_edges <- E(grain_graph)[V(grain_graph)[on_crit_path==1] %--% V(grain_graph)[on_crit_path==1]]
        grain_graph <- set.edge.attribute(grain_graph, name="on_crit_path", index=critical_edges, value=1)
        if (cl_args$verbose) {
            ctr <- ctr + 1
            setTxtProgressBar(progress_bar, ctr)
            close(progress_bar)
        }
    }
    #Rprof(NULL)

    if (cl_args$timing) toc("Critical path calculation")
}

# Write basic grain graph info
sink(grain_graph_info_out_file, append=T)
my_print("# Grain graph structure")
my_print(paste("Number of nodes =", length(V(grain_graph))))
my_print(paste("Number of edges =", length(E(grain_graph))))
my_print(paste("Number of tasks =", length(prof_data$task)))
my_print(paste("Number of forks =", length(fork_nodes_unique)))
my_print("# Out-degree distribution of forks")
degree.distribution(grain_graph, v=fork_nodes_index, mode="out")
my_print(paste("# Cilk theory parallelism (metric = ", path_weight, ")", sep=""))
my_print(paste("Span (critical path) =", critical_path))
work <- sum(as.numeric(prof_data[,path_weight]))
my_print(paste("Work =", work))
my_print(paste("Parallelism (Work/Span) =", work/critical_path))
if (cl_args$enumcriticalpath)
    my_print(paste("Number of critical tasks =", length(graph_vertices$task[graph_vertices$on_crit_path == 1])))
my_print()
sink()

# Write grain graph to file
if (cl_args$verbose) my_print("Writing grain graph to various file formats ...")

## Layout in Sugiyama style and write to PDF
if (cl_args$layout) {
    if (cl_args$timing) tic(type="elapsed")
    temp_out_file <- paste(gsub(". $", "", cl_args$out), ".pdf", sep="")
    lyt <- layout_with_sugiyama(grain_graph, attributes="all")
    pdf(temp_out_file)
    res <- plot(grain_graph, layout=lyt$layout)
    res <- dev.off()
    my_print(paste("Wrote file:", temp_out_file))
    if (cl_args$timing) toc("Write Sugiyama layout PDF")
}

## Write dot file
#if (cl_args$timing) tic(type="elapsed")
#temp_out_file <- paste(gsub(". $", "", cl_args$out), ".dot", sep="")
#res <- write.graph(grain_graph, file=temp_out_file, format="dot")
#my_print(paste("Wrote file:", temp_out_file))
#if (cl_args$timing) toc("Write dot")

# Write graphml file
if (cl_args$timing) tic(type="elapsed")
temp_out_file <- paste(gsub(". $", "", cl_args$out), ".graphml", sep="")
res <- write.graph(grain_graph, file=temp_out_file, format="graphml")
my_print(paste("Wrote file:", temp_out_file))
if (cl_args$timing) toc("Write graphml")

# Write graphml file without attributes
if (cl_args$timing) tic(type="elapsed")
temp_out_file <- paste(gsub(". $", "", cl_args$out), "-noattrib.graphml", sep="")
grain_graph_noattrib <- grain_graph
for (attrib in vertex_attr_names(grain_graph)) {
    grain_graph_noattrib <- delete_vertex_attr(grain_graph_noattrib, attrib)
}
for (attrib in edge_attr_names(grain_graph)) {
    grain_graph_noattrib <- delete_edge_attr(grain_graph_noattrib, attrib)
}
res <- write.graph(grain_graph_noattrib, file=temp_out_file, format="graphml")
my_print(paste("Wrote file:", temp_out_file))
if (cl_args$timing) toc("Write graphml without attributes")

## Write adjacency matrix file
#if (cl_args$timing) tic(type="elapsed")
#temp_out_file <- paste(gsub(". $", "", cl_args$out), ".adjmat", sep="")
#sink(temp_out_file)
#get.adjacency(grain_graph,names=T)
#sink()
#my_print(paste("Wrote file:", temp_out_file))
#if (cl_args$timing) toc("Write adjacency matrix")

## Write edgelist file
#if (cl_args$timing) tic(type="elapsed")
#temp_out_file <- paste(gsub(". $", "", cl_args$out), ".edgelist", sep="")
#sink(temp_out_file)
#get.edgelist(grain_graph, names=T)
#sink()
#my_print(paste("Wrote file:", temp_out_file))
#if (cl_args$timing) toc("Write edgelist")

# Write node attributes
if (cl_args$timing) tic(type="elapsed")
temp_out_file <- paste(gsub(". $", "", cl_args$out), ".nodeattr", sep="")
write.table(get.data.frame(grain_graph, what="vertices"), sep=",", file=temp_out_file)
my_print(paste("Wrote file:", temp_out_file))
if (cl_args$timing) toc("Write node attributes")

my_print(paste("Wrote file:", grain_graph_info_out_file))

# Warn
wa <- warnings()
if (class(wa) != "NULL")
    print(wa)
