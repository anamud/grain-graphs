# Clean slate
rm(list=ls())

# Include
mir_root <- Sys.getenv("GRAIN_GRAPHS_ROOT")
source(paste(mir_root,"/prototype/common.R",sep=""))

# Parse args
Rstudio_mode <- F
if (Rstudio_mode) {
    cl_args <- list(data="task-stats.processed",
                   out="task-graph",
                   enumcriticalpath=F,
                   showproblems=T,
                   grainproblemconfig="grain-problems.cfg",
                   grainpropertyconfig="grain-properties.cfg",
                   edgepropertyconfig="edge-properties.cfg",
                   verbose=T,
                   timing=F,
                   layout=F)
} else {
    option_list <- list(make_option(c("-d","--data"), help = "Task stats.", metavar="FILE"),
                        make_option(c("-o","--out"), default="task-graph", help = "Output file suffix [default \"%default\"].", metavar="STRING"),
                        make_option(c("--enumcriticalpath"), action="store_true", default=FALSE, help="Enumerate critical path."),
                        make_option(c("--forloop"), action="store_true", default=FALSE, help="Task stats obtained from a for-loop program."),
                        make_option(c("--showproblems"), action="store_true", default=FALSE, help="Analyze graph for problems."),
                        make_option(c("--grainproblemconfig"), default="grain-problems.cfg", help = "Grain problem configuration file [default \"%default\"].", metavar="FILE"),
                        make_option(c("--grainpropertyconfig"), default="grain-properties.cfg", help = "Grain property configuration file [default \"%default\"].", metavar="FILE"),
                        make_option(c("--edgepropertyconfig"), default="edge-properties.cfg", help = "Edge property configuration file [default \"%default\"].", metavar="FILE"),
                        make_option(c("--verbose"), action="store_true", default=TRUE, help="Print output [default]."),
                        make_option(c("--quiet"), action="store_false", dest="verbose", help="Print little output."),
                        make_option(c("--timing"), action="store_true", default=FALSE, help="Print processing time."),
                        make_option(c("--layout"), action="store_true", default=FALSE, help="Layout using Sugiyama style and plot to PDF."))

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
tg_data <- read.csv(cl_args$data, header=TRUE)

# Information output
tg_info_out_file <- paste(gsub(". $", "", cl_args$out), ".info", sep="")
sink(tg_info_out_file)
sink()

# Remove background task
tg_data <- tg_data[!is.na(tg_data$parent),]

# Path weight avaiability
if (path_weight %in% colnames(tg_data))
    path_weight <- NA

if (cl_args$forloop) {
    # Remove idle task without children
    tg_data <- tg_data[!(tg_data$tag == "idle_task" & tg_data$num_children == 0),]
}

if (cl_args$verbose) my_print("Creating graph ...")

# Create node lists
if (cl_args$timing) tic(type="elapsed")

# Create join nodes list
join_nodes <- mapply(function(x, y, z) {paste('j', x, y, sep='.')}, x=tg_data$parent, y=tg_data$joins_at)
join_nodes_unique <- unique(unlist(join_nodes, use.names=FALSE))

# Create parent nodes list
parent_nodes_unique <- unique(tg_data$parent)

# Create fork nodes list
fork_nodes <- mapply(function(x, y, z) {paste('f', x, y, sep='.')}, x=tg_data$parent, y=tg_data$joins_at)
fork_nodes_unique <- unique(unlist(fork_nodes, use.names=FALSE))

if (cl_args$timing) toc("Node list creation")

# Create graph
if (cl_args$timing) tic(type="elapsed")

tg <- graph.empty(directed=TRUE) + vertices('E',
                                            unique(c(join_nodes_unique,
                                                     fork_nodes_unique,
                                                     parent_nodes_unique,
                                                     tg_data$task)))

if (cl_args$timing) toc("Graph creation")

# Connect parent fork to task
if (cl_args$verbose) my_print("Connecting nodes ...")
if (cl_args$timing) tic(type="elapsed")

tg[from=fork_nodes, to=tg_data$task, attr='kind'] <- 'create'
tg[from=fork_nodes, to=tg_data$task, attr='color'] <- create_edge_color

if (cl_args$timing) toc("Connect parent fork to task")

# Connect parent task to first fork
if (cl_args$timing) tic(type="elapsed")

first_forks_index <- which(grepl("f.[0-9]+.0$", fork_nodes_unique))
parent_first_forks <- as.vector(sapply(fork_nodes_unique[first_forks_index], function(x) {gsub('f.(.*)\\.+.*','\\1', x)}))
first_forks <- fork_nodes_unique[first_forks_index]

tg[to=first_forks, from=parent_first_forks, attr='kind'] <- 'scope'
tg[to=first_forks, from=parent_first_forks, attr='color'] <- scope_edge_color

if (!is.na(path_weight)) {
    temp <- as.numeric(tg_data[match(parent_first_forks, tg_data$task),path_weight])

    tg[to=first_forks, from=parent_first_forks, attr=path_weight] <- temp
    tg[to=first_forks, from=parent_first_forks, attr='weight'] <- -temp
}

if (cl_args$timing) toc("Connect parent to first fork")

# Connect leaf task to join
if (cl_args$timing) tic(type="elapsed")

leaf_tasks <- tg_data$task[tg_data$leaf == T]
leaf_join_nodes <- join_nodes[match(leaf_tasks, tg_data$task)]

tg[from=leaf_tasks, to=leaf_join_nodes, attr='kind'] <- 'sync'
tg[from=leaf_tasks, to=leaf_join_nodes, attr='color'] <- sync_edge_color

if (!is.na(path_weight)) {
    temp <- as.numeric(tg_data[match(leaf_tasks, tg_data$task),path_weight])
    tg[from=leaf_tasks, to=leaf_join_nodes, attr=path_weight] <- temp
    tg[from=leaf_tasks, to=leaf_join_nodes, attr='weight'] <- -temp
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
        parent_index <- match(parent, tg_data$task)
        gfather <- tg_data[parent_index,]$parent
        gfather_join <- paste('j', as.character(gfather), as.character(tg_data[parent_index,]$joins_at), sep=".")

        if (is.na(match(gfather_join, join_nodes_unique)) == F) {
            next_fork <- gfather_join # Connect to grandfather's join
        } else {
            next_fork <- 'E' # Connect to end node
        }
    }
    next_fork
}

next_forks <- as.vector(sapply(join_nodes_unique, find_next_fork))

tg[from=join_nodes_unique, to=next_forks, attr='kind'] <- 'continue'
tg[from=join_nodes_unique, to=next_forks, attr='color'] <- cont_edge_color

#Rprof(NULL)
if (cl_args$timing) toc("Connect join to next fork")

# Set attributes
if (cl_args$verbose) my_print("Setting attributes ...")
if (cl_args$timing) tic(type="elapsed")

# Common vertex attributes
V(tg)$label <- V(tg)$name

# Set task vertex attributes
task_index <- match(as.character(tg_data$task), V(tg)$name)

# Set annotations
for (annot in colnames(tg_data)) {
    values <- as.character(tg_data[,annot])
    tg <- set.vertex.attribute(tg, name=annot, index=task_index, value=values)
}

# Set size constants
tg <- set.vertex.attribute(tg, name='size', index=task_index, value=task_size)
tg <- set.vertex.attribute(tg, name='width', index=task_index, value=task_size)
tg <- set.vertex.attribute(tg, name='height', index=task_index, value=task_size)
tg <- set.vertex.attribute(tg, name='shape', index=task_index, value=task_shape)

# Scale size to attributes
size_scaled <- c("ins_count", "work_cycles", "overhead_cycles", "exec_cycles")
for(attrib in size_scaled) {
    if (attrib %in% colnames(tg_data)) {

        # Set size
        attrib_unique <- unique(tg_data[,attrib])
        if (length(attrib_unique) == 1) {
            p_task_size <- task_size_mult
        } else if(length(attrib_unique) == 2 & identical(c(0,NA), as.numeric(attrib_unique[order(attrib_unique)]))) {
            p_task_size <- task_size_mult
        } else {
            p_task_size <- task_size_mult * as.numeric(cut(tg_data[,attrib], task_size_bins))
        }
        annot_name <- paste(attrib, "_to_size", sep="")
        tg <- set.vertex.attribute(tg, name=annot_name, index=task_index, value=p_task_size)

        # Set height
        attrib_val <- tg_data[,attrib]
        attrib_val_norm <- 1 + ((attrib_val - min(attrib_val)) / (max(attrib_val) - min(attrib_val)))
        annot_name <- paste(attrib, "_to_height", sep="")
        tg <- set.vertex.attribute(tg, name=annot_name, index=task_index, value=attrib_val_norm*task_size)
    }
}

# Set color constants
tg <- set.vertex.attribute(tg, name='color', index=task_index, value=task_color)

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

    if (attrib %in% colnames(tg_data)) {

        # Set color in proportion to attrib
        attrib_unique <- unique(tg_data[,attrib])
        if (invert_colors) {
            if (length(attrib_unique) == 1) {
                p_task_color <- task_color_pal[task_color_bins]
            } else if(length(attrib_unique) == 2 & identical(c(0,NA), as.numeric(attrib_unique[order(attrib_unique)]))) {
                p_task_color <- task_color_pal[task_color_bins]
            } else {
                p_task_color <- rev(task_color_pal)[as.numeric(cut(tg_data[,attrib], task_color_bins))]
            }
        } else {
            if (length(attrib_unique) == 1) {
                p_task_color <- task_color_pal[1]
            } else if(length(attrib_unique) == 2 & identical(c(0,NA), as.numeric(attrib_unique[order(attrib_unique)]))) {
                p_task_color <- task_color_pal[1]
            } else {
                p_task_color <- task_color_pal[as.numeric(cut(tg_data[,attrib], task_color_bins))]
            }
        }
        annot_name <- paste(attrib, "_to_color", sep="")
        tg <- set.vertex.attribute(tg, name=annot_name, index=task_index, value=p_task_color)

        # Write colors for reference
        tg_out_file <- paste(gsub(". $", "", cl_args$out), annot_name, sep=".")
        if (length(attrib_unique) == 1) {
            write.csv(data.frame(value=attrib_unique, color=p_task_color), tg_out_file, row.names=F)
        } else if(length(attrib_unique) == 2 & identical(c(0,NA), as.numeric(attrib_unique[order(attrib_unique)]))) {
            write.csv(data.frame(value=attrib_unique, color=p_task_color), tg_out_file, row.names=F)
        } else {
            v <- unique(cut(tg_data[,attrib], task_color_bins))
            write.csv(data.frame(value=v, color=task_color_pal[as.numeric(v)]), tg_out_file, row.names=F)
        }
        my_print(paste("Wrote file:", tg_out_file))
    }
}

# Set attributes to distinct color
attrib_color_distinct <- c("cpu_id", "outl_func", "tag", "outline_function")
for(attrib in attrib_color_distinct) {
    if (attrib %in% colnames(tg_data)) {

        # Map distinct color to attrib
        attrib_val <- as.character(tg_data[,attrib])
        unique_attrib_val <- unique(attrib_val)
        #attrib_color <- rainbow(length(unique_attrib_val), start=0.4, end=0.8)
        attrib_color <- rev(topo.colors(length(unique_attrib_val)))
        annot_name <- paste(attrib, "_to_color", sep="")
        tg <- set.vertex.attribute(tg, name=annot_name, index=task_index, value=attrib_color[match(attrib_val, unique_attrib_val)])

        # Write colors for reference
        tg_out_file <- paste(gsub(". $", "", cl_args$out), annot_name, sep=".")
        write.csv(data.frame(value=unique_attrib_val, color=attrib_color), tg_out_file, row.names=F)
        my_print(paste("Wrote file:", tg_out_file))
    }
}

# Set label and color of 'task 0'
start_index <- V(tg)$name == '0'
tg <- set.vertex.attribute(tg, name='color', index=start_index, value=start_color)
tg <- set.vertex.attribute(tg, name='label', index=start_index, value='S')
tg <- set.vertex.attribute(tg, name='size', index=start_index, value=start_size)
tg <- set.vertex.attribute(tg, name='shape', index=start_index, value=start_shape)

# Set label and color of 'task E'
end_index <- V(tg)$name == "E"
tg <- set.vertex.attribute(tg, name='color', index=end_index, value=end_color)
tg <- set.vertex.attribute(tg, name='label', index=end_index, value='E')
tg <- set.vertex.attribute(tg, name='size', index=end_index, value=end_size)
tg <- set.vertex.attribute(tg, name='shape', index=end_index, value=end_shape)

# Set fork vertex attributes
fork_nodes_index <- startsWith(V(tg)$name, 'f')
tg <- set.vertex.attribute(tg, name='size', index=fork_nodes_index, value=fork_size)
tg <- set.vertex.attribute(tg, name='color', index=fork_nodes_index, value=fork_color)
tg <- set.vertex.attribute(tg, name='label', index=fork_nodes_index, value='^')
tg <- set.vertex.attribute(tg, name='shape', index=fork_nodes_index, value=fork_shape)

# Set join vertex attributes
join_nodes_index <- startsWith(V(tg)$name, 'j')
tg <- set.vertex.attribute(tg, name='size', index=join_nodes_index, value=join_size)
tg <- set.vertex.attribute(tg, name='color', index=join_nodes_index, value=join_color)
tg <- set.vertex.attribute(tg, name='label', index=join_nodes_index, value='*')
tg <- set.vertex.attribute(tg, name='shape', index=join_nodes_index, value=join_shape)

# Set edge attributes
if (!is.na(path_weight)) {
    tg <- set.edge.attribute(tg, name="weight", index=which(is.na(E(tg)$weight)), value=0)
    tg <- set.edge.attribute(tg, name=path_weight, index=which(is.na(get.edge.attribute(tg, name=path_weight, index=E(tg)))), value=0)
}

if (cl_args$timing) toc("Attribute setting")

# Check if graph has bad structure
if (cl_args$verbose) my_print("Checking for bad structure ...")
if (cl_args$timing) tic(type="elapsed")

if (is.element(0, degree(tg, fork_nodes_index, mode = c("in")))) {
    my_print("Warning! One or more fork nodes have zero degree since one or more tasks in the program performed empty synchronization.")
    my_print("Aborting on error!")
    quit("no", 1)
}
if (is.element(0, degree(tg, fork_nodes_index, mode = c("out")))) {
    my_print("Warning! One or more fork nodes have zero degree since one or more tasks in the program performed empty synchronization.")
    my_print("Aborting on error!")
    quit("no", 1)
}
if (is.element(0, degree(tg, join_nodes_index, mode = c("in")))) {
    my_print("Warning! One or more join nodes have zero degree since one or more tasks in the program performed empty synchronization.")
    my_print("Aborting on error!")
    quit("no", 1)
}
if (is.element(0, degree(tg, join_nodes_index, mode = c("out")))) {
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
    #tg <- simplify(tg, edge.attr.comb=toString)

    # Get critical path
    #Rprof("profile-critpathcalc.out")
    if (!cl_args$enumcriticalpath) {
        # Compute critical path length
        sp <- shortest.paths(tg, v=start_index, to=end_index, mode="out")
        lpl <- -as.numeric(sp)
    } else {
        # TODO: Make variable names in this block meaningfull.
        lntg <- length(V(tg))
        if (cl_args$verbose) {
            pb <- txtProgressBar(min = 0, max = lntg, style = 3)
            ctr <- 0
        }
        # Topological sort
        tsg <- topological.sort(tg)
        # Set root path attributes
        V(tg)[tsg[1]]$rdist <- 0
        V(tg)[tsg[1]]$depth <- 0
        V(tg)[tsg[1]]$rpath <- tsg[1]
        # Get data frame of graph object
        vgdf <- get.data.frame(tg, what="vertices")
        # Get longest paths from root
        for(node in tsg[-1])
        {
            # Get distance from node's predecessors
            ni <- incident(tg, node, mode="in")
            w <- -E(tg)[ni]$weight
            # Get distance from root to node's predecessors
            nn <- neighbors(tg, node, mode="in")
            d <- vgdf$rdist[nn]
            # Add distances (assuming one-one corr.)
            wd <- w+d
            # Set node's distance from root to max of added distances
            mwd <- max(wd)
            vgdf$rdist[node] <- mwd
            # Set node's path from root to path of max of added distances
            mwdn <- as.vector(nn)[match(mwd,wd)]
            nrp <- list(c(unlist(vgdf$rpath[mwdn]), node))
            vgdf$rpath[node] <- nrp
            # Set node's depth as one greater than the largest depth its predecessors
            vgdf$depth[node] <- max(vgdf$depth[nn]) + 1
            if (cl_args$verbose) {
                ctr <- ctr + 1
                setTxtProgressBar(pb, ctr)
            }
        }
        ## Longest path is the largest root distance
        lpl <- max(vgdf$rdist)
        # Enumerate longest path
        lpm <- unlist(vgdf$rpath[match(lpl,vgdf$rdist)])
        vgdf$on_crit_path <- 0
        vgdf$on_crit_path[lpm] <- 1
        # Set back on graph
        tg <- set.vertex.attribute(tg, name="on_crit_path", index=V(tg), value=vgdf$on_crit_path)
        tg <- set.vertex.attribute(tg, name="rdist", index=V(tg), value=vgdf$rdist)
        tg <- set.vertex.attribute(tg, name="depth", index=V(tg), value=vgdf$depth)
        critical_edges <- E(tg)[V(tg)[on_crit_path==1] %--% V(tg)[on_crit_path==1]]
        tg <- set.edge.attribute(tg, name="on_crit_path", index=critical_edges, value=1)
        if (cl_args$verbose) {
            ctr <- ctr + 1
            setTxtProgressBar(pb, ctr)
            close(pb)
        }
    }
    #Rprof(NULL)

    # Calculate and write info
    sink(tg_info_out_file, append=T)
    my_print(paste("# Cilk theory parallelism (unit =", path_weight, "):", sep=""))
    my_print(paste("Span (critical path) =", lpl))
    work <- sum(as.numeric(tg_data[,path_weight]))
    my_print(paste("Work =", work))
    my_print(paste("Parallelism (Work/Span) =", work/lpl))
    my_print()
    sink()

    # Shape calculation
    if (cl_args$enumcriticalpath) {
        # Clear rpath since dot/table writing complains
        tg <- remove.vertex.attribute(tg,"rpath")

        # Calc shape
        tg_df <- get.data.frame(tg, what="vertices")
        tg_df <- tg_df[!is.na(as.numeric(tg_df$label)),]
        #tg_shape_interval_width <- work/(length(unique(tg_data$cpu_id))*mean(tg_data[,path_weight]))
        tg_shape_interval_width <- median(as.numeric(tg_df[,path_weight], na.rm=T))
        stopifnot(tg_shape_interval_width > 0)
        tg_shape_breaks <- seq(0, max(tg_df$rdist) + 1 + tg_shape_interval_width, by=tg_shape_interval_width)
        tg_shape <- hist(tg_df$rdist, breaks=tg_shape_breaks, plot=F)

        # Write out shape
        tg_out_file <- paste(gsub(". $", "", cl_args$out), "-shape.pdf", sep="")
        pdf(tg_out_file)
        plot(tg_shape, freq=T, xlab=paste("Elapsed ", path_weight), ylab="Tasks", main="Instantaneous task parallelism", col="white")
        abline(h = length(unique(tg_data$cpu_id)), col = "blue", lty=2)
        abline(h = work/lpl , col = "red", lty=1)
        legend("top", legend = c("Number of cores", "Exposed task parallelism"), fill = c("blue", "red"))
        dev.off()
        my_print(paste("Wrote file:", tg_out_file))
    }
    if (cl_args$timing) toc("Critical path calculation")
}

# Write basic graph info
sink(tg_info_out_file, append=T)
my_print("# Task graph structure:")
my_print(paste("Number of nodes =", length(V(tg))))
my_print(paste("Number of edges =", length(E(tg))))
my_print(paste("Number of tasks =", length(tg_data$task)))
if (cl_args$enumcriticalpath)
    my_print(paste("Number of critical tasks =", length(tg_df$task[tg_df$on_crit_path == 1])))
my_print(paste("Number of forks =", length(fork_nodes_unique)))
my_print("Out-degree distribution of forks:")
degree.distribution(tg, v=fork_nodes_index, mode="out")
my_print()
sink()

# Write graph to file
if (cl_args$verbose) my_print("Writing graph files ...")

## Layout in Sugiyama style and write to PDF
if (cl_args$layout) {
    if (cl_args$timing) tic(type="elapsed")
    tg_out_file <- paste(gsub(". $", "", cl_args$out), ".pdf", sep="")
    lyt <- layout_with_sugiyama(tg, attributes="all")
    pdf(tg_out_file)
    res <- plot(tg, layout=lyt$layout)
    res <- dev.off()
    my_print(paste("Wrote file:", tg_out_file))
    if (cl_args$timing) toc("Write Sugiyama layout PDF")
}

## Write dot file
#if (cl_args$timing) tic(type="elapsed")
#tg_out_file <- paste(gsub(". $", "", cl_args$out), ".dot", sep="")
#res <- write.graph(tg, file=tg_out_file, format="dot")
#my_print(paste("Wrote file:", tg_out_file))
#if (cl_args$timing) toc("Write dot")

# Write gml file
if (cl_args$timing) tic(type="elapsed")
tg_out_file <- paste(gsub(". $", "", cl_args$out), ".graphml", sep="")
res <- write.graph(tg, file=tg_out_file, format="graphml")
my_print(paste("Wrote file:", tg_out_file))
if (cl_args$timing) toc("Write graphml")

# Write graphml file with no attributes
if (cl_args$timing) tic(type="elapsed")
tg_out_file <- paste(gsub(". $", "", cl_args$out), "-noattrib.graphml", sep="")
tg_noattrib <- tg
for (attrib in vertex_attr_names(tg)) {
    tg_noattrib <- delete_vertex_attr(tg_noattrib, attrib)
}
for (attrib in edge_attr_names(tg)) {
    tg_noattrib <- delete_edge_attr(tg_noattrib, attrib)
}
res <- write.graph(tg_noattrib, file=tg_out_file, format="graphml")
my_print(paste("Wrote file:", tg_out_file))
if (cl_args$timing) toc("Write graphml without attributes")

## Write adjacency matrix file
#if (cl_args$timing) tic(type="elapsed")
#tg_out_file <- paste(gsub(". $", "", cl_args$out), ".adjmat", sep="")
#sink(tg_out_file)
#get.adjacency(tg,names=T)
#sink()
#my_print(paste("Wrote file:", tg_out_file))
#if (cl_args$timing) toc("Write adjacency matrix")

## Write edgelist file
#if (cl_args$timing) tic(type="elapsed")
#tg_out_file <- paste(gsub(". $", "", cl_args$out), ".edgelist", sep="")
#sink(tg_out_file)
#get.edgelist(tg, names=T)
#sink()
#my_print(paste("Wrote file:", tg_out_file))
#if (cl_args$timing) toc("Write edgelist")

# Write node attributes
if (cl_args$timing) tic(type="elapsed")
tg_out_file <- paste(gsub(". $", "", cl_args$out), ".nodeattr", sep="")
write.table(get.data.frame(tg, what="vertices"), sep=",", file=tg_out_file)
my_print(paste("Wrote file:", tg_out_file))
if (cl_args$timing) toc("Write node attributes")

# Show problems
if (cl_args$showproblems) {
    if (cl_args$verbose) my_print("Analyzing graph for problems ...")
    if (cl_args$timing) tic(type="elapsed")

    # Base task graph with transparent elements
    base_tg <- tg
    base_tg_vertex_color <- add.alpha(get.vertex.attribute(base_tg, name='color'), alpha=0.2)
    base_tg <- set.vertex.attribute(base_tg, name='color', value=base_tg_vertex_color)
    base_tg <- set.vertex.attribute(base_tg, name='problematic', value=0)
    base_tg_edge_color <- add.alpha(get.edge.attribute(base_tg, name='color'), alpha=0.2)
    base_tg <- set.edge.attribute(base_tg, name='color', value=base_tg_edge_color)
    # Set base tg edge colors to gray
    #base_tg <- set.edge.attribute(base_tg, name='color', value="#c0c0c0")

    # problem text output
    tg_problem_out_file <- paste(gsub(". $", "", cl_args$out), "-problem.info", sep="")
    sink(tg_problem_out_file)
    my_print("# Problems:")
    my_print()
    sink()

    # Memory hierarchy utilization problem
    if ("mem_hier_util" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        mem_hier_util_thresh <- 0.5
        prob_task <- subset(tg_data, mem_hier_util > mem_hier_util_thresh, select=task)

        sink(tg_problem_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have mem_hier_util >", mem_hier_util_thresh))
        sink()

        if (cl_args$enumcriticalpath) {
            prob_task_critical <- subset(tg_df, mem_hier_util > mem_hier_util_thresh & on_crit_path == 1, select=task)

            sink(tg_problem_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_problem_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='mem_hier_util_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", cl_args$out), "-problem-memory-hierarchy-utilization.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Memory footprint problem
    if ("mem_fp" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        mem_fp_thresh <- 512000
        prob_task <- subset(tg_data, mem_fp > mem_fp_thresh, select=task)

        sink(tg_problem_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have mem_fp >", mem_fp_thresh))
        sink()

        if (cl_args$enumcriticalpath) {
            prob_task_critical <- subset(tg_df, mem_fp > mem_fp_thresh & on_crit_path == 1, select=task)

            sink(tg_problem_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_problem_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='mem_fp_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", cl_args$out), "-problem-memory-footprint.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Compute intensity problem
    if ("compute_int" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        compute_int_thresh <- 2
        prob_task <- subset(tg_data, compute_int < compute_int_thresh, select=task)

        sink(tg_problem_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have compute_int <", compute_int_thresh))
        sink()

        if (cl_args$enumcriticalpath) {
            prob_task_critical <- subset(tg_df, compute_int < compute_int_thresh & on_crit_path == 1, select=task)

            sink(tg_problem_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_problem_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='compute_int_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", cl_args$out), "-problem-compute-intensity.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Work deviation problem
    if ("work_deviation" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        work_deviation_thresh <- 2
        prob_task <- subset(tg_data, work_deviation > work_deviation_thresh, select=task)

        sink(tg_problem_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have work_deviation >", work_deviation_thresh))
        sink()

        if (cl_args$enumcriticalpath) {
            prob_task_critical <- subset(tg_df, work_deviation > work_deviation_thresh & on_crit_path == 1, select=task)

            sink(tg_problem_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_problem_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='work_deviation_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", cl_args$out), "-problem-work-deviation.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Parallel benefit problem
    if ("parallel_benefit" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        parallel_benefit_thresh <- 1
        prob_task <- subset(tg_data, parallel_benefit < parallel_benefit_thresh, select=task)

        sink(tg_problem_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have parallel_benefit <", parallel_benefit_thresh))
        sink()

        if (cl_args$enumcriticalpath) {
            prob_task_critical <- subset(tg_df, parallel_benefit < parallel_benefit_thresh & on_crit_path == 1, select=task)

            sink(tg_problem_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_problem_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='parallel_benefit_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", cl_args$out), "-problem-parallel-benefit.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Parallelism problem
    if (cl_args$enumcriticalpath) {# {{{
        prob_tg <- base_tg
        parallelism_thresh <- length(unique(tg_data$cpu_id))
        ranges <- which(tg_shape$counts > 0 && tg_shape$counts < parallelism_thresh)

        sink(tg_problem_out_file, append=T)
        my_print(paste(length(ranges), "shape bins out of", length(tg_shape$counts), "have parallelism <", parallelism_thresh))
        sink()

        for (r in ranges) {
            prob_task <- subset(tg_df, rdist < tg_shape$breaks[r+1] & rdist > tg_shape$breaks[r], select=label)
            prob_task_index <- match(as.character(prob_task$label), V(prob_tg)$name)
            #prob_task_color <- get.vertex.attribute(prob_tg, name='cpu_id_to_color', index=prob_task_index)
            prob_task_color <- "#FF0000"

            prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
            prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)
            # TODO: Highlight range.
        }

        sink(tg_problem_out_file, append=T)
        my_print()
        sink()

        tg_out_file <- paste(gsub(". $", "", cl_args$out), "-problem-parallelism.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Shape contribution problem (median)
    if ("median_shape_contrib" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        shape_contrib_thresh <- length(unique(tg_data$cpu_id))
        prob_task <- subset(tg_data, median_shape_contrib < shape_contrib_thresh, select=task)

        sink(tg_problem_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have median_shape_contrib <", shape_contrib_thresh))
        sink()

        if (cl_args$enumcriticalpath) {
            prob_task_critical <- subset(tg_df, median_shape_contrib < shape_contrib_thresh & on_crit_path == 1, select=task)

            sink(tg_problem_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_problem_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='median_shape_contrib_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", cl_args$out), "-problem-median-shape-contrib.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Shape contribution problem (min)
    if ("min_shape_contrib" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        shape_contrib_thresh <- length(unique(tg_data$cpu_id))
        prob_task <- subset(tg_data, min_shape_contrib < shape_contrib_thresh, select=task)

        sink(tg_problem_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have min_shape_contrib <", shape_contrib_thresh))
        sink()

        if (cl_args$enumcriticalpath) {
            prob_task_critical <- subset(tg_df, min_shape_contrib < shape_contrib_thresh & on_crit_path == 1, select=task)

            sink(tg_problem_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_problem_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)

        # Set color in proportion to shape_attrib of all tasks
        #prob_task_color <- get.vertex.attribute(prob_tg, name='min_shape_contrib_to_color', index=prob_task_index)

        # Set color in proportion to shape_attrib of problem tasks
        prob_task_shape_contrib <- as.numeric(get.vertex.attribute(prob_tg, name='min_shape_contrib', index=prob_task_index))
        prob_task_shape_contrib_unique <- unique(prob_task_shape_contrib)
        invert_colors <- 1

        breaks <- pretty(prob_task_shape_contrib, n=10)
        num_bins <- length(breaks)
        color_pal <- color_fun(num_bins)
        if (invert_colors) {
            if (length(prob_task_shape_contrib_unique) == 1) {
                prob_task_color <- color_pal[num_bins]
            } else if (length(prob_task_shape_contrib_unique) == 2 & identical(c(0,NA), as.numeric(prob_task_shape_contrib_unique[order(prob_task_shape_contrib_unique)]))) {
                prob_task_color <- color_pal[num_bins]
            } else {
                #prob_task_color <- rev(color_pal)[as.numeric(cut(prob_task_shape_contrib, num_bins))]
                prob_task_color <- rev(color_pal)[as.numeric(cut(prob_task_shape_contrib, breaks=breaks))]
            }
        } else {
            if (length(prob_task_shape_contrib_unique) == 1) {
                prob_task_color <- color_pal[1]
            } else if (length(prob_task_shape_contrib_unique) == 2 & identical(c(0,NA), as.numeric(prob_task_shape_contrib_unique[order(prob_task_shape_contrib_unique)]))) {
                prob_task_color <- color_pal[1]
            } else {
                #prob_task_color <- color_pal[as.numeric(cut(prob_task_shape_contrib, num_bins))]
                prob_task_color <- color_pal[as.numeric(cut(prob_task_shape_contrib, breaks=breaks))]
            }
        }
        # Write colors for reference
        tg_out_file <- paste(gsub(". $", "", cl_args$out), "prob_task_min_shape_contrib_to_color", sep=".")
        if (length(prob_task_shape_contrib_unique) == 1) {
            write.csv(data.frame(value=prob_task_shape_contrib_unique, color=prob_task_color), tg_out_file, row.names=F)
        } else if (length(prob_task_shape_contrib_unique) == 2 & identical(c(0,NA), as.numeric(prob_task_shape_contrib_unique[order(prob_task_shape_contrib_unique)]))) {
            write.csv(data.frame(value=prob_task_shape_contrib_unique, color=prob_task_color), tg_out_file, row.names=F)
        } else {
            #v <- unique(cut(prob_task_shape_contrib, num_bins))
            v <- unique(cut(prob_task_shape_contrib, breaks=breaks))
            if (invert_colors) {
                write.csv(data.frame(value=v, color=rev(color_pal)[as.numeric(v)]), tg_out_file, row.names=F)
            } else {
                write.csv(data.frame(value=v, color=color_pal[as.numeric(v)]), tg_out_file, row.names=F)
            }
        }
        my_print(paste("Wrote file:", tg_out_file))

        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", cl_args$out), "-problem-min-shape-contrib.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Shape contribution problem (max)
    if ("max_shape_contrib" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        shape_contrib_thresh <- length(unique(tg_data$cpu_id))
        prob_task <- subset(tg_data, max_shape_contrib < shape_contrib_thresh, select=task)

        sink(tg_problem_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have max_shape_contrib <", shape_contrib_thresh))
        sink()

        if (cl_args$enumcriticalpath) {
            prob_task_critical <- subset(tg_df, max_shape_contrib < shape_contrib_thresh & on_crit_path == 1, select=task)

            sink(tg_problem_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_problem_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='max_shape_contrib_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", cl_args$out), "-problem-max-shape-contrib.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    if(!cl_args$forloop) {
        # Sibling work balance problem
        if ("sibling_work_balance" %in% colnames(tg_data)) {# {{{
            prob_tg <- base_tg
            sibling_work_balance_thresh <- 2
            prob_task <- subset(tg_data, sibling_work_balance > sibling_work_balance_thresh, select=task)

            sink(tg_problem_out_file, append=T)
            my_print(paste(length(prob_task$task), "tasks have sibling_work_balance >", sibling_work_balance_thresh))
            sink()

            if (cl_args$enumcriticalpath) {
                prob_task_critical <- subset(tg_df, sibling_work_balance > sibling_work_balance_thresh & on_crit_path == 1, select=task)

                sink(tg_problem_out_file, append=T)
                my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
                sink()
            }

            sink(tg_problem_out_file, append=T)
            my_print()
            sink()

            prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
            prob_task_color <- get.vertex.attribute(prob_tg, name='sibling_work_balance_to_color', index=prob_task_index)
            prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
            prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

            tg_out_file <- paste(gsub(". $", "", cl_args$out), "-problem-sibling-work-balance.graphml", sep="")
            res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
            my_print(paste("Wrote file:", tg_out_file))
        }# }}}

        # Sibling scatter problem
        if ("sibling_scatter" %in% colnames(tg_data)) {# {{{
            prob_tg <- base_tg
            sibling_scatter_thresh <- (length(unique(tg_data$cpu_id))/4)
            prob_task <- subset(tg_data, sibling_scatter > sibling_scatter_thresh, select=task)

            sink(tg_problem_out_file, append=T)
            my_print(paste(length(prob_task$task), "tasks have sibling_scatter >", sibling_scatter_thresh))
            sink()

            if (cl_args$enumcriticalpath) {
                prob_task_critical <- subset(tg_df, sibling_scatter > sibling_scatter_thresh & on_crit_path == 1, select=task)

                sink(tg_problem_out_file, append=T)
                my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
                sink()
            }

            sink(tg_problem_out_file, append=T)
            my_print()
            sink()

            prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
            #prob_task_color <- get.vertex.attribute(prob_tg, name='sibling_scatter_to_color', index=prob_task_index)
            prob_task_color <- get.vertex.attribute(prob_tg, name='cpu_id_to_color', index=prob_task_index)
            prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
            prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

            tg_out_file <- paste(gsub(". $", "", cl_args$out), "-problem-sibling-scatter.graphml", sep="")
            res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
            my_print(paste("Wrote file:", tg_out_file))
        }# }}}
    }

    if(cl_args$forloop) {
        # Chunk work balance problem
        if ("chunk_work_balance" %in% colnames(tg_data)) {# {{{
            prob_tg <- base_tg
            chunk_work_balance_thresh <- 2
            prob_task <- subset(tg_data, !is.na(chunk_work_balance) & chunk_work_balance > chunk_work_balance_thresh, select=task)

            sink(tg_analysis_out_file, append=T)
            my_print(paste(length(prob_task$task), "tasks have chunk_work_balance >", chunk_work_balance_thresh))
            sink()

            if (!parsed$cplengthonly) {
                prob_task_critical <- subset(tg_df, !is.na(chunk_work_balance) & chunk_work_balance > chunk_work_balance_thresh & on_crit_path == 1, select=task)

                sink(tg_analysis_out_file, append=T)
                my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
                sink()
            }

            sink(tg_analysis_out_file, append=T)
            my_print()
            sink()

            prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
            prob_task_color <- get.vertex.attribute(prob_tg, name='chunk_work_balance_to_color', index=prob_task_index)
            prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
            prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

            tg_out_file <- paste(gsub(". $", "", parsed$out), "-problem-chunk-work-balance.graphml", sep="")
            res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
            my_print(paste("Wrote file:", tg_out_file))
        }# }}}

        # Chunk work CPU balance problem
        if ("chunk_work_cpu_balance" %in% colnames(tg_data)) {# {{{
            prob_tg <- base_tg
            chunk_work_cpu_balance_thresh <- 2
            prob_task <- subset(tg_data, !is.na(chunk_work_cpu_balance) & chunk_work_cpu_balance > chunk_work_cpu_balance_thresh, select=task)

            sink(tg_analysis_out_file, append=T)
            my_print(paste(length(prob_task$task), "tasks have chunk_work_cpu_balance >", chunk_work_cpu_balance_thresh))
            sink()

            if (!parsed$cplengthonly) {
                prob_task_critical <- subset(tg_df, !is.na(chunk_work_cpu_balance) & chunk_work_cpu_balance > chunk_work_cpu_balance_thresh & on_crit_path == 1, select=task)

                sink(tg_analysis_out_file, append=T)
                my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
                sink()
            }

            sink(tg_analysis_out_file, append=T)
            my_print()
            sink()

            prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
            prob_task_color <- get.vertex.attribute(prob_tg, name='chunk_work_cpu_balance_to_color', index=prob_task_index)
            prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
            prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

            tg_out_file <- paste(gsub(". $", "", parsed$out), "-problem-chunk-work-cpu-balance.graphml", sep="")
            res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
            my_print(paste("Wrote file:", tg_out_file))
        }# }}}
    }

    my_print(paste("Wrote file:", tg_problem_out_file))
    if (cl_args$timing) toc("Analyzing graph for problems")
}

my_print(paste("Wrote file:", tg_info_out_file))

# Warn
wa <- warnings()
if (class(wa) != "NULL")
    print(wa)
