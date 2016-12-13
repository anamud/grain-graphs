# Clean slate
rm(list=ls())

# Include
mir_root <- Sys.getenv("GRAIN_GRAPHS_ROOT")
source(paste(mir_root,"prototype/common.R",sep=""))

# Graph element sizes
join_size <- 10
fork_size <- join_size
fork_size_mult <- 10
fork_size_bins <- 10
start_size <- 15
end_size <- start_size
task_size <- 30
task_size_mult <- 10
task_size_bins <- 10

# Graph element shapes
task_shape <- "rectangle"
fork_shape <- "circle"
join_shape <- fork_shape
start_shape <- fork_shape
end_shape <- start_shape

# Parse args
Rstudio_mode <- F
if (Rstudio_mode) {
    parsed <- list(data="task-stats.processed",
                   palette="color",
                   out="task-graph",
                   tree=F,
                   cplengthonly=F,
                   analyze=T,
                   config="task-graph-analysis.cfg",
                   verbose=T,
                   timing=F,
                   layout=F)
} else {
    option_list <- list(
                        make_option(c("-d","--data"), help = "Task stats.", metavar="FILE"),
                        make_option(c("-p","--palette"), default="color", help = "Color palette for graph elements [default \"%default\"]."),
                        make_option(c("-o","--out"), default="task-graph", help = "Output file suffix [default \"%default\"].", metavar="STRING"),
                        make_option(c("-t", "--tree"), action="store_true", default=FALSE, help="Plot task graph as tree."),
                        make_option(c("--cplengthonly"), action="store_true", default=FALSE, help="Calculate critical path length only. Skip critical path enumeration."),
                        make_option(c("--analyze"), action="store_true", default=FALSE, help="Analyze task graph for problems."),
                        make_option(c("--config"), default="task-graph-analysis.cfg", help = "Analysis configuration file [default \"%default\"].", metavar="FILE"),
                        make_option(c("--verbose"), action="store_true", default=TRUE, help="Print output [default]."),
                        make_option(c("--quiet"), action="store_false", dest="verbose", help="Print little output."),
                        make_option(c("--timing"), action="store_true", default=FALSE, help="Print processing time."),
                        make_option(c("--layout"), action="store_true", default=FALSE, help="Layout using Sugiyama style and plot to PDF."))

    parsed <- parse_args(OptionParser(option_list = option_list), args = commandArgs(TRUE))

    if (!exists("data", where=parsed)) {
        my_print("Error: Invalid arguments. Check help (-h)")
        quit("no", 1)
    }
}

if (parsed$verbose) my_print("Initializing ...")

# Read data
tg_data <- read.csv(parsed$data, header=TRUE)

# Information output
tg_info_out_file <- paste(gsub(". $", "", parsed$out), ".info", sep="")
sink(tg_info_out_file)
sink()

# Remove background task
tg_data <- tg_data[!is.na(tg_data$parent),]

# Critical path calculation weight
if ("ins_count" %in% colnames(tg_data)) {
    path_weight <- "ins_count"
} else if ("work_cycles" %in% colnames(tg_data)) {
    path_weight <- "work_cycles"
} else if ("exec_cycles" %in% colnames(tg_data)) {
    path_weight <- "exec_cycles"
} else {
    path_weight <- NA
}

# Set colors
join_color <- "#FF7F50"  # coral
fork_color <- "#2E8B57"  # seagreen
task_color <- "#4682B4" #steelblue
other_color <- "#DEB887" # burlywood
create_edge_color <- fork_color
sync_edge_color <- join_color
scope_edge_color <- "#000000"
cont_edge_color <- "#000000"
color_fun <- colorRampPalette(c("yellow", "orange"))

if (parsed$palette == "gray") {
    join_color <- "#D3D3D3"  # light gray
    fork_color <- "#D3D3D3"  # light gray
    task_color <- "#6B6B6B"  # gray42
    other_color <- "#D3D3D3" # light gray
    create_edge_color <- "#000000"
    sync_edge_color <- "#000000"
    scope_edge_color <-  "#000000"
    cont_edge_color <- "#000000"
    color_fun <- colorRampPalette(c("gray10", "gray90"))
} else if (parsed$palette != "color") {
    my_print("Unsupported color format. Supported formats: color, gray. Defaulting to color.")
}

# Task color binning
task_color_bins <- 10
task_color_pal <- color_fun(task_color_bins)
fork_color_bins <- 10
fork_color_pal <- color_fun(fork_color_bins)

if (parsed$verbose) my_print("Creating graph ...")

# Create node lists
if (parsed$timing) tic(type="elapsed")

if (!parsed$tree) {

    # Create join nodes list
    join_nodes <- mapply(function(x, y, z) {paste('j', x, y, sep='.')}, x=tg_data$parent, y=tg_data$joins_at)
    join_nodes_unique <- unique(unlist(join_nodes, use.names=FALSE))
}

# Create parent nodes list
parent_nodes_unique <- unique(tg_data$parent)

# Create fork nodes list
fork_nodes <- mapply(function(x, y, z) {paste('f', x, y, sep='.')}, x=tg_data$parent, y=tg_data$joins_at)
fork_nodes_unique <- unique(unlist(fork_nodes, use.names=FALSE))

if (parsed$timing) toc("Node list creation")

# Create graph
if (parsed$timing) tic(type="elapsed")

if (!parsed$tree) {
    tg <- graph.empty(directed=TRUE) + vertices('E',
                                                unique(c(join_nodes_unique,
                                                         fork_nodes_unique,
                                                         parent_nodes_unique,
                                                         tg_data$task)))
} else {
    tg <- graph.empty(directed=TRUE) + vertices('E',
                                                unique(c(fork_nodes_unique,
                                                         parent_nodes_unique,
                                                         tg_data$task)))
}

if (parsed$timing) toc("Graph creation")

# Connect parent fork to task
if (parsed$verbose) my_print("Connecting nodes ...")
if (parsed$timing) tic(type="elapsed")

tg[from=fork_nodes, to=tg_data$task, attr='kind'] <- 'create'
tg[from=fork_nodes, to=tg_data$task, attr='color'] <- create_edge_color

if (parsed$timing) toc("Connect parent fork to task")

# Connect parent task to first fork
if (parsed$timing) tic(type="elapsed")

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

if (parsed$timing) toc("Connect parent to first fork")

if (!parsed$tree) {

    # Connect leaf task to join
    if (parsed$timing) tic(type="elapsed")

    leaf_tasks <- tg_data$task[tg_data$leaf == T]
    leaf_join_nodes <- join_nodes[match(leaf_tasks, tg_data$task)]

    tg[from=leaf_tasks, to=leaf_join_nodes, attr='kind'] <- 'sync'
    tg[from=leaf_tasks, to=leaf_join_nodes, attr='color'] <- sync_edge_color

    if (!is.na(path_weight)) {
        temp <- as.numeric(tg_data[match(leaf_tasks, tg_data$task),path_weight])
        tg[from=leaf_tasks, to=leaf_join_nodes, attr=path_weight] <- temp
        tg[from=leaf_tasks, to=leaf_join_nodes, attr='weight'] <- -temp
    }

    if (parsed$timing) toc("Connect leaf task to join")

    # Connect join to next fork
    if (parsed$timing) tic(type="elapsed")

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
    if (parsed$timing) toc("Connect join to next fork")
} else {

    # Connect fork to next fork
    if (parsed$timing) tic(type="elapsed")

    #Rprof("profile-forktonext.out")
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
            next_fork <- node # Connect to myself (self-loop)
        }
        next_fork
    }

    next_forks <- as.vector(sapply(fork_nodes_unique, find_next_fork))

    tg[from=fork_nodes_unique, to=next_forks, attr='kind'] <- 'continue'
    tg[from=fork_nodes_unique, to=next_forks, attr='color'] <- cont_edge_color

    #Rprof(NULL)

    if (parsed$timing) toc("Connect fork to next fork")

    # Connext E to last fork of task 0
    if (parsed$timing) tic(type="elapsed")

    get_join_count <- function(node)
    {
        node_split <- unlist(strsplit(node, "\\."))
        parent <- as.numeric(node_split[2])
        join_count <- as.numeric(node_split[3])
        join_count
    }

    fork_nodes_of_zero <- fork_nodes_unique[which(grepl("f.0.[0-9]+$", fork_nodes_unique))]
    largest_join_count_of_zero <- max(as.vector(sapply(fork_nodes_of_zero, get_join_count)))

    tg[from=paste("f.0.",largest_join_count_of_zero,sep=""), to='E', attr='kind'] <- 'continue'
    tg[from=paste("f.0.",largest_join_count_of_zero,sep=""), to='E', attr='color'] <- cont_edge_color

    if (parsed$timing) toc("Connect last fork of 0 to node E")
}

# Set attributes
if (parsed$verbose) my_print("Setting attributes ...")
if (parsed$timing) tic(type="elapsed")

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
attrib_color_scaled <- c("mem_fp", "-compute_int", "PAPI_RES_STL_sum", "-mem_hier_util", "work_deviation", "overhead_deviation", "-parallel_benefit", "-min_shape_contrib", "-max_shape_contrib","-median_shape_contrib", "sibling_work_balance", "sibling_scatter")
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
        tg_out_file <- paste(gsub(". $", "", parsed$out), annot_name, sep=".")
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
        tg_out_file <- paste(gsub(". $", "", parsed$out), annot_name, sep=".")
        write.csv(data.frame(value=unique_attrib_val, color=attrib_color), tg_out_file, row.names=F)
        my_print(paste("Wrote file:", tg_out_file))
    }
}

# Set label and color of 'task 0'
start_index <- V(tg)$name == '0'
tg <- set.vertex.attribute(tg, name='color', index=start_index, value=other_color)
tg <- set.vertex.attribute(tg, name='label', index=start_index, value='S')
tg <- set.vertex.attribute(tg, name='size', index=start_index, value=start_size)
tg <- set.vertex.attribute(tg, name='shape', index=start_index, value=start_shape)

# Set label and color of 'task E'
end_index <- V(tg)$name == "E"
tg <- set.vertex.attribute(tg, name='color', index=end_index, value=other_color)
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
if (!parsed$tree) {
    join_nodes_index <- startsWith(V(tg)$name, 'j')
    tg <- set.vertex.attribute(tg, name='size', index=join_nodes_index, value=join_size)
    tg <- set.vertex.attribute(tg, name='color', index=join_nodes_index, value=join_color)
    tg <- set.vertex.attribute(tg, name='label', index=join_nodes_index, value='*')
    tg <- set.vertex.attribute(tg, name='shape', index=join_nodes_index, value=join_shape)
}

# Set edge attributes
if (!is.na(path_weight)) {
    tg <- set.edge.attribute(tg, name="weight", index=which(is.na(E(tg)$weight)), value=0)
    tg <- set.edge.attribute(tg, name=path_weight, index=which(is.na(get.edge.attribute(tg, name=path_weight, index=E(tg)))), value=0)
}

if (parsed$timing) toc("Attribute setting")

# Check if graph has bad structure
if (parsed$verbose) my_print("Checking for bad structure ...")
if (parsed$timing) tic(type="elapsed")

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
if (!parsed$tree) {
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
}

if (parsed$timing) toc("Checking for bad structure")

# Calculate critical path
if (!is.na(path_weight) && !parsed$tree) {
    if (parsed$verbose) my_print("Calculating critical path ...")
    if (parsed$timing) tic(type="elapsed")
    # Simplify - DO NOT USE. Fucks up the critical path analysis.
    #tg <- simplify(tg, edge.attr.comb=toString)

    # Get critical path
    #Rprof("profile-critpathcalc.out")
    if (parsed$cplengthonly) {
        # Get critical path length
        sp <- shortest.paths(tg, v=start_index, to=end_index, mode="out")
        lpl <- -as.numeric(sp)
    } else {
        # TODO: Make variable names in this block meaningfull.
        lntg <- length(V(tg))
        if (parsed$verbose) {
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
            if (parsed$verbose) {
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
        if (parsed$verbose) {
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

    if (!parsed$cplengthonly) {
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
        tg_out_file <- paste(gsub(". $", "", parsed$out), "-shape.pdf", sep="")
        pdf(tg_out_file)
        plot(tg_shape, freq=T, xlab=paste("Elapsed ", path_weight), ylab="Tasks", main="Instantaneous task parallelism", col="white")
        abline(h = length(unique(tg_data$cpu_id)), col = "blue", lty=2)
        abline(h = work/lpl , col = "red", lty=1)
        legend("top", legend = c("Number of cores", "Exposed task parallelism"), fill = c("blue", "red"))
        dev.off()
        my_print(paste("Wrote file:", tg_out_file))
    }
    if (parsed$timing) toc("Critical path calculation")
} else {
    if (parsed$verbose) my_print("Simplifying graph ...")
    if (parsed$timing) tic(type="elapsed")
    tg <- simplify(tg, remove.multiple=T, remove.loops=T)
    if (parsed$timing) toc("Simplify")
}

# Write basic graph info
sink(tg_info_out_file, append=T)
my_print("# Task graph structure:")
my_print(paste("Number of nodes =", length(V(tg))))
my_print(paste("Number of edges =", length(E(tg))))
my_print(paste("Number of tasks =", length(tg_data$task)))
if (!parsed$cplengthonly)
    my_print(paste("Number of critical tasks =", length(tg_df$task[tg_df$on_crit_path == 1])))
my_print(paste("Number of forks =", length(fork_nodes_unique)))
my_print("Out-degree distribution of forks:")
degree.distribution(tg, v=fork_nodes_index, mode="out")
my_print()
sink()

# Write graph to file
if (parsed$verbose) my_print("Writing graph files ...")

## Layout in Sugiyama style and write to PDF
if (parsed$layout) {
    if (parsed$timing) tic(type="elapsed")
    tg_out_file <- paste(gsub(". $", "", parsed$out), ".pdf", sep="")
    lyt <- layout_with_sugiyama(tg, attributes="all")
    pdf(tg_out_file)
    res <- plot(tg, layout=lyt$layout)
    res <- dev.off()
    my_print(paste("Wrote file:", tg_out_file))
    if (parsed$timing) toc("Write Sugiyama layout PDF")
}

## Write dot file
#if (parsed$timing) tic(type="elapsed")
#tg_out_file <- paste(gsub(". $", "", parsed$out), ".dot", sep="")
#res <- write.graph(tg, file=tg_out_file, format="dot")
#my_print(paste("Wrote file:", tg_out_file))
#if (parsed$timing) toc("Write dot")

# Write gml file
if (parsed$timing) tic(type="elapsed")
tg_out_file <- paste(gsub(". $", "", parsed$out), ".graphml", sep="")
res <- write.graph(tg, file=tg_out_file, format="graphml")
my_print(paste("Wrote file:", tg_out_file))
if (parsed$timing) toc("Write graphml")

# Write graphml file with no attributes
if (parsed$timing) tic(type="elapsed")
tg_out_file <- paste(gsub(". $", "", parsed$out), "-noattrib.graphml", sep="")
tg_noattrib <- tg
for (attrib in vertex_attr_names(tg)) {
    tg_noattrib <- delete_vertex_attr(tg_noattrib, attrib)
}
for (attrib in edge_attr_names(tg)) {
    tg_noattrib <- delete_edge_attr(tg_noattrib, attrib)
}
res <- write.graph(tg_noattrib, file=tg_out_file, format="graphml")
my_print(paste("Wrote file:", tg_out_file))
if (parsed$timing) toc("Write graphml without attributes")

## Write adjacency matrix file
#if (parsed$timing) tic(type="elapsed")
#tg_out_file <- paste(gsub(". $", "", parsed$out), ".adjmat", sep="")
#sink(tg_out_file)
#get.adjacency(tg,names=T)
#sink()
#my_print(paste("Wrote file:", tg_out_file))
#if (parsed$timing) toc("Write adjacency matrix")

## Write edgelist file
#if (parsed$timing) tic(type="elapsed")
#tg_out_file <- paste(gsub(". $", "", parsed$out), ".edgelist", sep="")
#sink(tg_out_file)
#get.edgelist(tg, names=T)
#sink()
#my_print(paste("Wrote file:", tg_out_file))
#if (parsed$timing) toc("Write edgelist")

# Write node attributes
if (parsed$timing) tic(type="elapsed")
tg_out_file <- paste(gsub(". $", "", parsed$out), ".nodeattr", sep="")
write.table(get.data.frame(tg, what="vertices"), sep=",", file=tg_out_file)
my_print(paste("Wrote file:", tg_out_file))
if (parsed$timing) toc("Write node attributes")

# Show problems
if (parsed$analyze) {
    if (parsed$verbose) my_print("Analyzing graph for problems ...")
    if (parsed$timing) tic(type="elapsed")

    # Base task graph with transparent elements
    base_tg <- tg
    base_tg_vertex_color <- add.alpha(get.vertex.attribute(base_tg, name='color'), alpha=0.2)
    base_tg <- set.vertex.attribute(base_tg, name='color', value=base_tg_vertex_color)
    base_tg <- set.vertex.attribute(base_tg, name='problematic', value=0)
    if (!parsed$tree) {
        base_tg_edge_color <- add.alpha(get.edge.attribute(base_tg, name='color'), alpha=0.2)
        base_tg <- set.edge.attribute(base_tg, name='color', value=base_tg_edge_color)
        # Set base tg edge colors to gray
        #base_tg <- set.edge.attribute(base_tg, name='color', value="#c0c0c0")
    }

    # Analysis text output
    tg_analysis_out_file <- paste(gsub(". $", "", parsed$out), "-analysis.info", sep="")
    sink(tg_analysis_out_file)
    my_print("# Analysis:")
    my_print()
    sink()

    # Memory hierarchy utilization problem
    if ("mem_hier_util" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        mem_hier_util_thresh <- 0.5
        prob_task <- subset(tg_data, mem_hier_util > mem_hier_util_thresh, select=task)

        sink(tg_analysis_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have mem_hier_util >", mem_hier_util_thresh))
        sink()

        if (!parsed$cplengthonly) {
            prob_task_critical <- subset(tg_df, mem_hier_util > mem_hier_util_thresh & on_crit_path == 1, select=task)

            sink(tg_analysis_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_analysis_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='mem_hier_util_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", parsed$out), "-problem-memory-hierarchy-utilization.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Memory footprint problem
    if ("mem_fp" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        mem_fp_thresh <- 512000
        prob_task <- subset(tg_data, mem_fp > mem_fp_thresh, select=task)

        sink(tg_analysis_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have mem_fp >", mem_fp_thresh))
        sink()

        if (!parsed$cplengthonly) {
            prob_task_critical <- subset(tg_df, mem_fp > mem_fp_thresh & on_crit_path == 1, select=task)

            sink(tg_analysis_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_analysis_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='mem_fp_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", parsed$out), "-problem-memory-footprint.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Compute intensity problem
    if ("compute_int" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        compute_int_thresh <- 2
        prob_task <- subset(tg_data, compute_int < compute_int_thresh, select=task)

        sink(tg_analysis_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have compute_int <", compute_int_thresh))
        sink()

        if (!parsed$cplengthonly) {
            prob_task_critical <- subset(tg_df, compute_int < compute_int_thresh & on_crit_path == 1, select=task)

            sink(tg_analysis_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_analysis_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='compute_int_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", parsed$out), "-problem-compute-intensity.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Work deviation problem
    if ("work_deviation" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        work_deviation_thresh <- 2
        prob_task <- subset(tg_data, work_deviation > work_deviation_thresh, select=task)

        sink(tg_analysis_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have work_deviation >", work_deviation_thresh))
        sink()

        if (!parsed$cplengthonly) {
            prob_task_critical <- subset(tg_df, work_deviation > work_deviation_thresh & on_crit_path == 1, select=task)

            sink(tg_analysis_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_analysis_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='work_deviation_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", parsed$out), "-problem-work-deviation.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Parallel benefit problem
    if ("parallel_benefit" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        parallel_benefit_thresh <- 1
        prob_task <- subset(tg_data, parallel_benefit < parallel_benefit_thresh, select=task)

        sink(tg_analysis_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have parallel_benefit <", parallel_benefit_thresh))
        sink()

        if (!parsed$cplengthonly) {
            prob_task_critical <- subset(tg_df, parallel_benefit < parallel_benefit_thresh & on_crit_path == 1, select=task)

            sink(tg_analysis_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_analysis_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='parallel_benefit_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", parsed$out), "-problem-parallel-benefit.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Parallelism problem
    if (!parsed$cplengthonly && !parsed$tree) {# {{{
        prob_tg <- base_tg
        parallelism_thresh <- length(unique(tg_data$cpu_id))
        ranges <- which(tg_shape$counts > 0 && tg_shape$counts < parallelism_thresh)

        sink(tg_analysis_out_file, append=T)
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

        sink(tg_analysis_out_file, append=T)
        my_print()
        sink()

        tg_out_file <- paste(gsub(". $", "", parsed$out), "-problem-parallelism.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Shape contribution problem (median)
    if ("median_shape_contrib" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        shape_contrib_thresh <- length(unique(tg_data$cpu_id))
        prob_task <- subset(tg_data, median_shape_contrib < shape_contrib_thresh, select=task)

        sink(tg_analysis_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have median_shape_contrib <", shape_contrib_thresh))
        sink()

        if (!parsed$cplengthonly) {
            prob_task_critical <- subset(tg_df, median_shape_contrib < shape_contrib_thresh & on_crit_path == 1, select=task)

            sink(tg_analysis_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_analysis_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='median_shape_contrib_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", parsed$out), "-problem-median-shape-contrib.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Shape contribution problem (min)
    if ("min_shape_contrib" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        shape_contrib_thresh <- length(unique(tg_data$cpu_id))
        prob_task <- subset(tg_data, min_shape_contrib < shape_contrib_thresh, select=task)

        sink(tg_analysis_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have min_shape_contrib <", shape_contrib_thresh))
        sink()

        if (!parsed$cplengthonly) {
            prob_task_critical <- subset(tg_df, min_shape_contrib < shape_contrib_thresh & on_crit_path == 1, select=task)

            sink(tg_analysis_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_analysis_out_file, append=T)
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
        tg_out_file <- paste(gsub(". $", "", parsed$out), "prob_task_min_shape_contrib_to_color", sep=".")
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

        tg_out_file <- paste(gsub(". $", "", parsed$out), "-problem-min-shape-contrib.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Shape contribution problem (max)
    if ("max_shape_contrib" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        shape_contrib_thresh <- length(unique(tg_data$cpu_id))
        prob_task <- subset(tg_data, max_shape_contrib < shape_contrib_thresh, select=task)

        sink(tg_analysis_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have max_shape_contrib <", shape_contrib_thresh))
        sink()

        if (!parsed$cplengthonly) {
            prob_task_critical <- subset(tg_df, max_shape_contrib < shape_contrib_thresh & on_crit_path == 1, select=task)

            sink(tg_analysis_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_analysis_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='max_shape_contrib_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", parsed$out), "-problem-max-shape-contrib.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Sibling work balance problem
    if ("sibling_work_balance" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        sibling_work_balance_thresh <- 2
        prob_task <- subset(tg_data, sibling_work_balance > sibling_work_balance_thresh, select=task)

        sink(tg_analysis_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have sibling_work_balance >", sibling_work_balance_thresh))
        sink()

        if (!parsed$cplengthonly) {
            prob_task_critical <- subset(tg_df, sibling_work_balance > sibling_work_balance_thresh & on_crit_path == 1, select=task)

            sink(tg_analysis_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_analysis_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        prob_task_color <- get.vertex.attribute(prob_tg, name='sibling_work_balance_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", parsed$out), "-problem-sibling-work-balance.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    # Sibling scatter problem
    if ("sibling_scatter" %in% colnames(tg_data)) {# {{{
        prob_tg <- base_tg
        sibling_scatter_thresh <- (length(unique(tg_data$cpu_id))/4)
        prob_task <- subset(tg_data, sibling_scatter > sibling_scatter_thresh, select=task)

        sink(tg_analysis_out_file, append=T)
        my_print(paste(length(prob_task$task), "tasks have sibling_scatter >", sibling_scatter_thresh))
        sink()

        if (!parsed$cplengthonly) {
            prob_task_critical <- subset(tg_df, sibling_scatter > sibling_scatter_thresh & on_crit_path == 1, select=task)

            sink(tg_analysis_out_file, append=T)
            my_print(paste("    ", length(prob_task_critical$task), " of which are on the critical path."))
            sink()
        }

        sink(tg_analysis_out_file, append=T)
        my_print()
        sink()

        prob_task_index <- match(as.character(prob_task$task), V(prob_tg)$name)
        #prob_task_color <- get.vertex.attribute(prob_tg, name='sibling_scatter_to_color', index=prob_task_index)
        prob_task_color <- get.vertex.attribute(prob_tg, name='cpu_id_to_color', index=prob_task_index)
        prob_tg <- set.vertex.attribute(prob_tg, name='color', index=prob_task_index, value=prob_task_color)
        prob_tg <- set.vertex.attribute(prob_tg, name='problematic', index=prob_task_index, value=1)

        tg_out_file <- paste(gsub(". $", "", parsed$out), "-problem-sibling-scatter.graphml", sep="")
        res <- write.graph(prob_tg, file=tg_out_file, format="graphml")
        my_print(paste("Wrote file:", tg_out_file))
    }# }}}

    my_print(paste("Wrote file:", tg_analysis_out_file))
    if (parsed$timing) toc("Analyzing graph for problems")
}

my_print(paste("Wrote file:", tg_info_out_file))

# Warn
wa <- warnings()
if (class(wa) != "NULL")
    print(wa)
