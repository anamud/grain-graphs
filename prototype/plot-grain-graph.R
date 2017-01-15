#
# Setup
#

# Clean slate
rm(list=ls())

# Include support functions
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
                        make_option(c("-o","--out"), default="grain-graph", help = "Output file prefix [default \"%default\"].", metavar="STRING"),
                        make_option(c("--enumcriticalpath"), action="store_true", default=FALSE, help="Enumerate nodes on critical path."),
                        make_option(c("--forloop"), action="store_true", default=FALSE, help="Task profiling data obtained from a for-loop program."),
                        make_option(c("--unreduced"), action="store_true", default=FALSE, help="Plot unreduced graph (with fragments)"),
                        make_option(c("--overlap"), default="any", help = "Overlap type for instantaneous parallelism calculation. Choices: any, within. [default \"%default\"].", metavar="STRING"),
                        make_option(c("--layout"), action="store_true", default=FALSE, help="Layout using Sugiyama style and plot to PDF."),
                        make_option(c("--verbose"), action="store_true", default=TRUE, help="Print output [default]."),
                        make_option(c("--quiet"), action="store_false", dest="verbose", help="Print little output."),
                        make_option(c("--timing"), action="store_true", default=FALSE, help="Print processing time."))

    cl_args <- parse_args(OptionParser(option_list = option_list), args = commandArgs(TRUE))

    if (!exists("data", where=cl_args)) {
        my_print("Error: Invalid arguments. Check help (-h)")
        quit("no", 1)
    }

    if (cl_args$unreduced & cl_args$forloop) {
        my_print("Error: Unreduced graph for for-loop programs is not supported yet!")
        quit("no", 1)
    }
}

#
# Initialize varibles
#
if (cl_args$verbose) my_print("Initializing ...")
if (cl_args$timing) tic(type="elapsed")

# Prepare information output file
grain_graph_info_out_file <- paste(gsub(". $", "", cl_args$out), ".info", sep="")
sink(grain_graph_info_out_file)
sink()

# Read profiling data
prof_data <- read.csv(cl_args$data, header=TRUE, comment.char='#', na.strings="NA")

# Remove background task
prof_data <- prof_data[!is.na(prof_data$parent),]

if (cl_args$forloop) {
    # Remove idle task without children
    prof_data <- prof_data[!(prof_data$tag == "idle_task" & prof_data$num_children == 0),]
}

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

# Mapping utility functions
apply_task_size_mapping <- function(data,type,out_file=NA)
{
    data_unique <- unique(data)
    if (type == "direct") {
        ret_val <- data
        if(!is.na(out_file)) {
            write.csv(data.frame(value=data_unique, size=data_unique), out_file, row.names=F)
            my_print(paste("Wrote file:", out_file))
        }
    } else if (type == "linear") {
        linear_size_map_base_multiplier <- 30
        data_unique_norm <- 1 + ((data_unique - min(data_unique)) / (max(data_unique) - min(data_unique)))
        temp <- data_norm * linear_size_map_base_multiplier
        ret_val <-  temp[match(data, data_unique)]
        if(!is.na(out_file)) {
            write.csv(data.frame(value=data_unique, size=temp), out_file, row.names=F)
            my_print(paste("Wrote file:", out_file))
        }
    } else if (type == "linear-step") {
        linear_size_map_base_multiplier <- 30
        linear_size_map_num_steps <- 10
        if ((length(data_unique) == 1) || (length(data_unique) == 2 & identical(c(0,NA), as.numeric(data_unique[order(data_unique)])))) {
            ret_val <- linear_size_map_base_multiplier
            if(!is.na(out_file)) {
                write.csv(data.frame(value=data_unique, size=ret_val), out_file, row.names=F)
                my_print(paste("Wrote file:", out_file))
            }
        } else {
            temp <- cut(data, linear_size_map_num_steps)
            ret_val <- linear_size_map_base_multiplier * as.numeric(temp)
            if(!is.na(out_file)) {
                data_binned_unique <- unique(temp)
                write.csv(data.frame(value=data_binned_unique, size=linear_size_map_base_multiplier * as.numeric(data_binned_unique), out_file, row.names=F))
                my_print(paste("Wrote file:", out_file))
            }
        }
    } else {
        my_print(paste("Error: Unsupported task size mapping type", type))
        quit("no", 1)
    }

    return(ret_val)
}

apply_task_color_mapping <- function(data,type,out_file=NA)
{
    data_unique <- unique(data)
    if (type == "direct") {
        ret_val <- data
        if(!is.na(out_file)) {
            write.csv(data.frame(value=data_unique, color=data_unique), out_file, row.names=F)
            my_print(paste("Wrote file:", out_file))
        }
    } else if (type == "linear") {
        linear_color_map_palette <- rev(topo.colors(length(data_unique)))
        ret_val <-  linear_color_map_palette[match(data, data_unique)]
        if(!is.na(out_file)) {
            write.csv(data.frame(value=data_unique, color=linear_color_map_palette), out_file, row.names=F)
            my_print(paste("Wrote file:", out_file))
        }
    } else if (type == "linear-step") {
        color_func <- colorRampPalette(c("yellow", "orange"))
        #color_func <- colorRampPalette(c("gray10", "gray90"))
        linear_color_map_num_steps <- 10
        linear_color_map_palette <- color_func(linear_color_map_num_steps)
        if ((length(data_unique) == 1) || (length(data_unique) == 2 & identical(c(0,NA), as.numeric(data_unique[order(data_unique)])))) {
            ret_val <- linear_color_map_palette[1]
            if(!is.na(out_file)) {
                write.csv(data.frame(value=data_unique, color=ret_val), out_file, row.names=F)
                my_print(paste("Wrote file:", out_file))
            }
        } else {
            temp <- cut(data, linear_color_map_num_steps)
            ret_val <- linear_color_map_palette[as.numeric(temp)]
            if(!is.na(out_file)) {
                data_binned_unique <- unique(temp)
                write.csv(data.frame(value=data_binned_unique, color=linear_color_map_palette[as.numeric(data_binned_unique)]), out_file, row.names=F)
                my_print(paste("Wrote file:", out_file))
            }
        }
    } else {
        my_print(paste("Error: Unsupported task color mapping type", type))
        quit("no", 1)
    }

    return(ret_val)
}

apply_edge_weight_mapping <- function(data,type,out_file=NA)
{
    data_unique <- unique(data)
    if (type == "direct") {
        ret_val <- data
        if(!is.na(out_file)) {
            write.csv(data.frame(value=data_unique, weight=data_unique), out_file, row.names=F)
            my_print(paste("Wrote file:", out_file))
        }
    } else if (type == "linear") {
        my_print(paste("Error: Unsupported edge weight mapping type", type))
        quit("no", 1)
    } else if (type == "linear-step") {
        my_print(paste("Error: Unsupported edge weight mapping type", type))
        quit("no", 1)
    } else {
        my_print(paste("Error: Unsupported edge weight mapping type", type))
        quit("no", 1)
    }

    return(ret_val)
}


# Read property configuration files
grain_prop_cfg <- read.csv(cl_args$grainpropertyconfig, header=TRUE, comment.char='#', na.strings="NA")
edge_prop_cfg <- read.csv(cl_args$edgepropertyconfig, header=TRUE, comment.char='#', na.strings="NA")

# Set grain sizes
fork_width <- as.numeric(get_value(grain_prop_cfg, "fork", "width")[1])
fork_height <- as.numeric(get_value(grain_prop_cfg, "fork", "width")[1])
join_width <- as.numeric(get_value(grain_prop_cfg, "join", "width")[1])
join_height <- as.numeric(get_value(grain_prop_cfg, "join", "height")[1])
start_width <- as.numeric(get_value(grain_prop_cfg, "start", "width")[1])
start_height <- as.numeric(get_value(grain_prop_cfg, "start", "height")[1])
end_width <- as.numeric(get_value(grain_prop_cfg, "end", "width")[1])
end_height <- as.numeric(get_value(grain_prop_cfg, "end", "height")[1])

task_width <- get_value(grain_prop_cfg, "task", "width")
if (!is.na(task_width[2])) {
    if (!(task_width[1] %in% colnames(prof_data))) {
        my_print(paste("Error: Mapped variable", task_width[1], "not found in profiling data!"))
        quit("no", 1)
    }
}
task_height <- get_value(grain_prop_cfg, "task", "height")
if (!is.na(task_height[2])) {
    if (!(task_height[1] %in% colnames(prof_data))) {
        my_print(paste("Error: Mapped variable", task_height[1], "not found in profiling data!"))
        quit("no", 1)
    }
}

# Set grain shapes
fork_shape <- get_value(grain_prop_cfg, "fork", "shape")[1]
join_shape <- get_value(grain_prop_cfg, "join", "shape")[1]
start_shape <- get_value(grain_prop_cfg, "start", "shape")[1]
end_shape <- get_value(grain_prop_cfg, "end", "shape")[1]
task_shape <- get_value(grain_prop_cfg, "end", "shape")[1]

# Set grain colors
fork_color <- get_value(grain_prop_cfg, "fork", "color")[1]
join_color <- get_value(grain_prop_cfg, "join", "color")[1]
start_color <- get_value(grain_prop_cfg, "start", "color")[1]
end_color <- get_value(grain_prop_cfg, "end", "color")[1]

task_color <- get_value(grain_prop_cfg, "task", "color")
if (!is.na(task_color[2])) {
    if (!(task_color[1] %in% colnames(prof_data))) {
        my_print(paste("Error: Mapped variable", task_color[1], "not found in profiling data!"))
        quit("no", 1)
    }
}

# Set edge colors
create_edge_color <- get_value(edge_prop_cfg, "create", "color")[1]
sync_edge_color <- get_value(edge_prop_cfg, "sync", "color")[1]
scope_edge_color <- get_value(edge_prop_cfg, "scope", "color")[1]
cont_edge_color <- get_value(edge_prop_cfg, "continuation", "color")[1]

# Set edge weights
common_edge_weight <- get_value(edge_prop_cfg, "common", "weight")
if (!is.na(common_edge_weight[2])) {
    if (!(common_edge_weight[1] %in% colnames(prof_data))) {
        my_print(paste("Error: Mapped variable", common_edge_weight[1], "not found in profiling data!"))
        quit("no", 1)
    }
}

if (cl_args$timing) toc("Initializing")

#
# Create base graph structure
#
if (cl_args$verbose) my_print("Creating base grain graph structure ...")
if (cl_args$timing) tic(type="elapsed")

if (cl_args$unreduced) {
    # Compute join frequeny
    if (cl_args$timing) tic(type="elapsed")

    join_freq <- prof_data %>% arrange(parent, joins_at) %>% group_by(parent, joins_at) %>% summarise(count = n())

    if (cl_args$timing) toc("Computing join frequency")

    # Compute fragment chains for tasks
    if (cl_args$timing) tic(type="elapsed")

    fragmentize <- function (task, num_children, parent, child_number, joins_at)
    {
        num_fragments <- 1
        # Connect to parent fork
        new_fragment <- paste(task, num_fragments, sep='.')
        parent_fork <- paste('f', parent, child_number, sep='.')
        fragment_chain <- c(parent_fork, new_fragment)
        last_fragment <- new_fragment
        # Connect fragments
        if (num_children > 0)
        {
            joins <- join_freq[join_freq$parent == task,]$count
            joins_index <- 1
            num_forks <- 1
            num_joins <- 0
            while(1)
            {
                stopifnot(joins[joins_index] > 0)
                for(i in 1:joins[joins_index])
                {
                    num_fragments <- num_fragments + 1
                    new_fragment <- paste(task, num_fragments, sep='.')
                    fork <- paste('f', task, num_forks, sep=".")
                    fragment_chain <- c(fragment_chain, last_fragment, fork, fork, new_fragment)
                    num_forks <- num_forks + 1
                    last_fragment <- new_fragment
                    if (i == joins[joins_index])
                    {
                        num_fragments <- num_fragments + 1
                        new_fragment <- paste(task, num_fragments, sep='.')
                        join <- paste('j', task, num_joins, sep=".")
                        fragment_chain <- c(fragment_chain, last_fragment, join, join, new_fragment)
                        num_joins <- num_joins + 1
                        last_fragment <- new_fragment
                    }
                }
                joins_index <- joins_index + 1
                if (joins_index > length(joins))
                    break
            }
        }
        # Connect to parent join
        parent_join <- paste('j', parent, joins_at, sep='.')
        fragment_chain <- c(fragment_chain, last_fragment, parent_join)
        return(fragment_chain)
    }

    fragment_chains <- unlist(mapply(fragmentize,
                              task=prof_data$task,
                              num_children=prof_data$num_children,
                              parent=prof_data$parent,
                              child_number=prof_data$child_number,
                              joins_at=prof_data$joins_at))
    fragment_chains <- matrix(fragment_chains, nc=2, byrow=TRUE)

    if (cl_args$timing) toc("Computing fragment chains")

    # Construct graph
    if (cl_args$timing) tic(type="elapsed")

    grain_graph <- graph.edgelist(fragment_chains, directed = TRUE)

    if (cl_args$timing) toc("Constructing grain graph using fragment chains")
} else {
    # Create join nodes list
    join_nodes <- mapply(function(x, y, z) {paste('j', x, y, sep='.')}, x=prof_data$parent, y=prof_data$joins_at)
    join_nodes_unique <- unique(unlist(join_nodes, use.names=FALSE))

    # Create parent nodes list
    parent_nodes_unique <- unique(prof_data$parent)

    # Create fork nodes list
    fork_nodes <- mapply(function(x, y, z) {paste('f', x, y, sep='.')}, x=prof_data$parent, y=prof_data$joins_at)
    fork_nodes_unique <- unique(unlist(fork_nodes, use.names=FALSE))

    if (cl_args$timing) toc("Node list creation")

    # Place nodes in graph
    if (cl_args$timing) tic(type="elapsed")

    grain_graph <- graph.empty(directed=TRUE) + vertices('E',
                                                unique(c(join_nodes_unique,
                                                         fork_nodes_unique,
                                                         parent_nodes_unique,
                                                         prof_data$task)))

    if (cl_args$timing) toc("Adding nodes to grain graph")

    #
    # Connect nodes
    #
    if (cl_args$verbose) my_print("Connecting nodes ...")

    # Connect parent fork to task
    if (cl_args$timing) tic(type="elapsed")

    grain_graph[from=fork_nodes, to=prof_data$task, attr="type"] <- "create"
    grain_graph[from=fork_nodes, to=prof_data$task, attr="color"] <- create_edge_color

    if (cl_args$timing) toc("Connecting parent fork to task")

    # Connect parent task to first fork
    if (cl_args$timing) tic(type="elapsed")

    first_forks_index <- which(grepl("f.[0-9]+.0$", fork_nodes_unique))
    parent_first_forks <- as.vector(sapply(fork_nodes_unique[first_forks_index], function(x) {gsub("f.(.*)\\.+.*","\\1", x)}))
    first_forks <- fork_nodes_unique[first_forks_index]

    grain_graph[to=first_forks, from=parent_first_forks, attr="type"] <- "scope"
    grain_graph[to=first_forks, from=parent_first_forks, attr="color"] <- scope_edge_color

    if (!is.na(common_edge_weight[2])) {
        temp <- apply_edge_weight_mapping(as.numeric(prof_data[match(parent_first_forks, prof_data$task),common_edge_weight[1]]), common_edge_weight[2])
        grain_graph[to=first_forks, from=parent_first_forks, attr="weight"] <- -temp
    } else {
        grain_graph[to=first_forks, from=parent_first_forks, attr="weight"] <- -common_edge_weight[1]
    }

    if (cl_args$timing) toc("Connecting parent to first fork")

    # Connect leaf task to join
    if (cl_args$timing) tic(type="elapsed")

    leaf_tasks <- prof_data$task[prof_data$leaf == T]
    leaf_join_nodes <- join_nodes[match(leaf_tasks, prof_data$task)]

    grain_graph[from=leaf_tasks, to=leaf_join_nodes, attr="type"] <- "sync"
    grain_graph[from=leaf_tasks, to=leaf_join_nodes, attr="color"] <- sync_edge_color

    if (!is.na(common_edge_weight[2])) {
        temp <- apply_edge_weight_mapping(as.numeric(prof_data[match(leaf_tasks, prof_data$task),common_edge_weight[1]]), common_edge_weight[2])
        grain_graph[from=leaf_tasks, to=leaf_join_nodes, attr="weight"] <- -temp
    } else {
        grain_graph[from=leaf_tasks, to=leaf_join_nodes, attr="weight"] <- -common_edge_weight[1]
    }

    if (cl_args$timing) toc("Connecting leaf task to join")

    # Connect join to next fork
    if (cl_args$timing) tic(type="elapsed")

    #Rprof("profile-jointonext.out")
    find_next_fork <- function(node)
    {
        #my_print(paste("Processing node",node, sep=" "))

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

    grain_graph[from=join_nodes_unique, to=next_forks, attr="type"] <- "continue"
    grain_graph[from=join_nodes_unique, to=next_forks, attr="color"] <- cont_edge_color

    #Rprof(NULL)
    if (cl_args$timing) toc("Connecting join to next fork")
}

#
# Check if grain graph has bad structure
#
if (cl_args$verbose) my_print("Checking for bad structure ...")
if (cl_args$timing) tic(type="elapsed")

bad_structure <- 0

if (cl_args$unreduced) {
    bad_structure <- 0
} else {
    fork_nodes_index <- startsWith(V(grain_graph)$name, 'f')
    join_nodes_index <- startsWith(V(grain_graph)$name, 'j')
    if ((is.element(0, degree(grain_graph, fork_nodes_index, mode = c("in")))) ||
        (is.element(0, degree(grain_graph, fork_nodes_index, mode = c("out")))) ||
        (is.element(0, degree(grain_graph, join_nodes_index, mode = c("in")))) ||
        (is.element(0, degree(grain_graph, join_nodes_index, mode = c("out"))))) {
        my_print("Warning! One or more join nodes have zero degree since one or more tasks in the program performed empty synchronization.")
        bad_structure <- 1
    }
}

if (bad_structure == 1) {
    my_print("Graph has bad structure. Aborting on error!")
    quit("no", 1)
}

if (cl_args$timing) toc("Checking for bad structure")

#
# Set attributes
#
if (cl_args$verbose) my_print("Setting attributes ...")
if (cl_args$timing) tic(type="elapsed")

# Common attributes
V(grain_graph)$label <- V(grain_graph)$name

if (cl_args$unreduced) {
    # Set fragment grain attributes
    vertex_names <- get.vertex.attribute(grain_graph, name="name", index=V(grain_graph))
    fragment_index <- which(grepl("^[0-9]+.[0-9]+$", vertex_names))
    #start_index <- which(grepl("^f.0.1$", vertex_names))
    #end_index <- which(grepl("^j.0.0$", vertex_names))
    #fork_index <- which(grepl("^f.[0-9]+.[0-9]+$", vertex_names))
    #join_index <- which(grepl("^j.[0-9]+.[0-9]+$", vertex_names))

    grain_graph <- set.vertex.attribute(grain_graph, name="type", index=fragment_index, value="fragment")
    grain_graph <- set.vertex.attribute(grain_graph, name="shape", index=fragment_index, value=task_shape)

    fragment_tasks <- as.integer(get.vertex.attribute(grain_graph, name="name", index=fragment_index))
    for (annot in colnames(prof_data)) {
        if(annot %in% c("name","label","type","shape","exec_cycles")) {
           next
        }
        temp <- prof_data[match(fragment_tasks, prof_data$task),annot]
        if (is.numeric(temp)) {
            values <- as.numeric(temp)
        } else if (is.logical(temp)) {
            values <- as.logical(temp)
        } else {
            values <- as.character(temp)
        }
        grain_graph <- set.vertex.attribute(grain_graph, name=annot, index=fragment_index, value=values)
    }

    # Compute fragment execution cycles
    compute_fragment_duration <- function(task, wait, exec_cycles, choice)
    {
        # Each task has breaks at these instants: (execution start, child creation, child wait, execution end)
        # A fragments executes upto the next break.
        wait_instants <- as.numeric(unlist(strsplit(substring(wait, 2, nchar(wait)-1), ";", fixed = TRUE)))
        create_instants <- as.numeric(prof_data$create_instant[prof_data$parent == task])
        instants <- c(wait_instants, create_instants, 1, 1 + exec_cycles)
        # Sort to line up breaks.
        instants <- sort(instants)
        # Remove placeholder breaks.
        instants <- instants[instants != 0]
        # Assert if last break is not execution end
        stopifnot(instants[length(instants)] == 1 + exec_cycles)
        durations <- diff(instants)
        # Assert if there are no durations!
        stopifnot(length(durations) > 0)
        # Get fragment identifiers
        fragments <- paste(as.character(task),paste(".",as.character(seq(1:length(durations))),sep=""),sep="")
        # This is an ugly hack because mapply and do.call(rbind) are not working together.
        if (choice == 1)
            return(fragments)
        else
            return(durations)
    }

    fragment_duration <- data.table(fragment=unlist(mapply(compute_fragment_duration,
                                            task=prof_data$task,
                                            wait=prof_data$wait_instants,
                                            exec_cycles=prof_data$exec_cycles,
                                            choice=1)),
                     duration=unlist(mapply(compute_fragment_duration,
                                            task=prof_data$task,
                                            wait=prof_data$wait_instants,
                                            exec_cycles=prof_data$exec_cycles,
                                            choice=2)))

    fragment_exec_cycles <- fragment_duration$duration[match(fragment_duration$fragment, vertex_names[fragment_index])]
    grain_graph <- set.vertex.attribute(grain_graph, name="exec_cycles", index=V(grain_graph), value=0)
    grain_graph <- set.vertex.attribute(grain_graph, name="exec_cycles", index=fragment_index, value=fragment_exec_cycles)

    # Set fragment size to constant or based on execution cycles
    if (!is.na(task_width[2])) {
        if (task_width[1] != "exec_cycles") {
            my_print(paste("Error: Cannot map fragment width to", task_width[1], ". Available mapping options are: exec_cyles."))
            quit("no", 1)
        } else {
            temp <- apply_task_size_mapping(fragment_exec_cycles, task_width[2])
            grain_graph <- set.vertex.attribute(grain_graph, name="height", index=fragment_index, value=temp)
        }
    } else {
        grain_graph <- set.vertex.attribute(grain_graph, name="width", index=fragment_index, value=task_width[1])
    }
    if (!is.na(task_height[2])) {
        if (task_height[1] != "exec_cycles") {
            my_print(paste("Error: Cannot map fragment height to", task_height[1], ". Available mapping options are: exec_cyles."))
            quit("no", 1)
        } else {
            temp <- apply_task_size_mapping(fragment_exec_cycles, task_height[2])
            grain_graph <- set.vertex.attribute(grain_graph, name="height", index=fragment_index, value=temp)
        }
    } else {
        grain_graph <- set.vertex.attribute(grain_graph, name="height", index=fragment_index, value=task_height[1])
    }

    # Set fragment color to constant or based on execution cycles, CPU, outline function, tag, or source line
    if (!is.na(task_color[2])) {
        sensible_mappings <- c("exec_cycles","cpu_id","outline_function","tag","source_line")
        if (!(task_color[1] %in% sensible_mappings)) {
            temp <- paste(sensible_mappings, collapse = ", ")
            my_print(paste("Error: Cannot map fragment color to", task_color[1], ". Available mapping options are:", temp))
            quit("no", 1)
        } else {
            if (task_color[1] == "exec_cycles") {
                temp <- apply_task_color_mapping(fragment_exec_cycles, task_color[2], paste("task-", task_color[1], "-", task_color[2], ".colormap", sep=""))
            } else {
                temp <- apply_task_color_mapping(as.character(get.vertex.attribute(grain_graph, name=task_color[1], index=fragment_index)), task_color[2], paste("task-", task_color[1], "-", task_color[2], ".colormap", sep=""))
            }
            grain_graph <- set.vertex.attribute(grain_graph, name="color", index=fragment_index, value=temp)
        }
    } else {
        grain_graph <- set.vertex.attribute(grain_graph, name="color", index=fragment_index, value=task_color[1])
    }

    # Set edge weight to constant or based on execution cycles
    if (!is.na(common_edge_weight[2])) {
        if (common_edge_weight[1] != "exec_cycles") {
            my_print(paste("Error: Cannot map edge weight to", common_edge_weight[1], ". Available mapping options are: exec_cyles."))
            quit("no", 1)
        } else {
            if (cl_args$verbose) {
                my_print("[subprocess] Setting edge weight attribute ...")
                num_vertices <- length(V(grain_graph))
                progress_bar <- txtProgressBar(min = 0, max = num_vertices, style = 3)
                ctr <- 0
            }
            top_sort_graph <- topological.sort(grain_graph)
            for(node in top_sort_graph[-1])
            {
                incident_edges <- incident(grain_graph, node, mode="in")
                incident_edge_vertices <- get.edges(grain_graph, incident_edges)
                incident_edge_weights <- V(grain_graph)[incident_edge_vertices[,1]]$exec_cycles
                grain_graph <- set.edge.attribute(grain_graph, name="weight", index=incident_edges, value=incident_edge_weights)
                if (cl_args$verbose) {
                    ctr <- ctr + 1
                    setTxtProgressBar(progress_bar, ctr)
                }
            }
            if (cl_args$verbose) {
                close(progress_bar)
            }
        }
    } else {
        top_sort_graph <- topological.sort(grain_graph)
        for(node in top_sort_graph[-1])
        {
            incident_edges <- incident(grain_graph, node, mode="in")
            grain_graph <- set.edge.attribute(grain_graph, name="weight", index=incident_edges, value=common_edge_weight[1])
        }
    }

    # Set edge type and color
    grain_graph <- set.edge.attribute(grain_graph, name="type", index=E(grain_graph), value="scope")
    grain_graph <- set.edge.attribute(grain_graph, name="color", index=E(grain_graph), value=scope_edge_color)
} else {
    # Set task grain attributes
    task_index <- match(as.character(prof_data$task), V(grain_graph)$name)
    grain_graph <- set.vertex.attribute(grain_graph, name="shape", index=task_index, value=task_shape)
    grain_graph <- set.vertex.attribute(grain_graph, name="type", index=task_index, value="task")

    # Set annotations
    for (annot in colnames(prof_data)) {
        temp <- prof_data[,annot]
        if (is.numeric(temp)) {
            values <- as.numeric(temp)
        } else if (is.logical(temp)) {
            values <- as.logical(temp)
        } else {
            values <- as.character(temp)
        }
        grain_graph <- set.vertex.attribute(grain_graph, name=annot, index=task_index, value=values)
    }
}

# Set attributes of start grain
start_index <- V(grain_graph)$name == '0'
grain_graph <- set.vertex.attribute(grain_graph, name="color", index=start_index, value=start_color)
grain_graph <- set.vertex.attribute(grain_graph, name="label", index=start_index, value='S')
grain_graph <- set.vertex.attribute(grain_graph, name="shape", index=start_index, value=start_shape)
grain_graph <- set.vertex.attribute(grain_graph, name="type", index=start_index, value="start")
grain_graph <- set.vertex.attribute(grain_graph, name="width", index=start_index, value=start_width)
grain_graph <- set.vertex.attribute(grain_graph, name="height", index=start_index, value=start_height)

# Set attributes of end grain
end_index <- V(grain_graph)$name == "E"
grain_graph <- set.vertex.attribute(grain_graph, name="color", index=end_index, value=end_color)
grain_graph <- set.vertex.attribute(grain_graph, name="label", index=end_index, value='E')
grain_graph <- set.vertex.attribute(grain_graph, name="shape", index=end_index, value=end_shape)
grain_graph <- set.vertex.attribute(grain_graph, name="type", index=end_index, value="end")
grain_graph <- set.vertex.attribute(grain_graph, name="width", index=end_index, value=end_width)
grain_graph <- set.vertex.attribute(grain_graph, name="height", index=end_index, value=end_height)

# Set fork grain attributes
fork_nodes_index <- startsWith(V(grain_graph)$name, 'f')
grain_graph <- set.vertex.attribute(grain_graph, name="color", index=fork_nodes_index, value=fork_color)
grain_graph <- set.vertex.attribute(grain_graph, name="label", index=fork_nodes_index, value='^')
grain_graph <- set.vertex.attribute(grain_graph, name="shape", index=fork_nodes_index, value=fork_shape)
grain_graph <- set.vertex.attribute(grain_graph, name="type", index=fork_nodes_index, value="fork")
grain_graph <- set.vertex.attribute(grain_graph, name="width", index=fork_nodes_index, value=fork_width)
grain_graph <- set.vertex.attribute(grain_graph, name="height", index=fork_nodes_index, value=fork_height)

# Set join grain attributes
join_nodes_index <- startsWith(V(grain_graph)$name, 'j')
grain_graph <- set.vertex.attribute(grain_graph, name="color", index=join_nodes_index, value=join_color)
grain_graph <- set.vertex.attribute(grain_graph, name="label", index=join_nodes_index, value='*')
grain_graph <- set.vertex.attribute(grain_graph, name="shape", index=join_nodes_index, value=join_shape)
grain_graph <- set.vertex.attribute(grain_graph, name="type", index=join_nodes_index, value="join")
grain_graph <- set.vertex.attribute(grain_graph, name="width", index=join_nodes_index, value=join_width)
grain_graph <- set.vertex.attribute(grain_graph, name="height", index=join_nodes_index, value=join_height)

# Set edge attributes
# Set weight to zero for edges not assigned weight before (essential for critical path calculation)
grain_graph <- set.edge.attribute(grain_graph, name="weight", index=which(is.na(E(grain_graph)$weight)), value=0)

if (cl_args$timing) toc("Setting attributes")

#
# Calculate critical path
#
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
        if (cl_args$unreduced) {
            incident_edge_weights <- V(grain_graph)[get.edges(grain_graph, incident_edges)[,1]]$exec_cycles
        } else {
            incident_edge_weights <- -E(grain_graph)[incident_edges]$weight
        }
        # Get distance from root to node's predecessors
        adjacent_nodes <- neighbors(grain_graph, node, mode="in")
        adjacent_nodes_root_dist <- graph_vertices$root_dist[adjacent_nodes]
        # Add distances (assuming one-one corr.)
        root_dist <- incident_edge_weights + adjacent_nodes_root_dist
        # Set node's distance from root to max of added distances
        max_root_dist <- max(root_dist)
        graph_vertices$root_dist[node] <- max_root_dist
        # Set node's path from root to path of max of added distances
        nodes_on_root_path <- as.vector(adjacent_nodes)[match(max_root_dist,root_dist)]
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
    grain_graph <- set.edge.attribute(grain_graph, name="color", index=critical_edges, value="#FF0000")
    grain_graph <- set.vertex.attribute(grain_graph, name="border-color", index=critical_nodes, value="#FF0000")
    # Cleanup
    grain_graph <- remove.vertex.attribute(grain_graph,"root_path")
    if (cl_args$verbose) {
        ctr <- ctr + 1
        setTxtProgressBar(progress_bar, ctr)
        close(progress_bar)
    }
}
#Rprof(NULL)

if (cl_args$timing) toc("Calculating critical path")

if (cl_args$unreduced) {
    #
    # Calculate instantaneous parellelism
    # Note: Instantaneous parallelism is found by counting overlapping ranges.
    # See question http://stackoverflow.com/questions/30978837/histogram-like-summary-for-interval-data
    #

    if (cl_args$verbose) my_print("Calculating instantaneous parellelism ...")
    if (cl_args$timing) tic(type="elapsed")

    # Select only fragments
    graph_fragments <- graph_vertices[with(graph_vertices, grepl("^[0-9]+.[0-9]+$", name)),]

    # Compute execution ranges
    # Each fragment has exection range [root_dist, root_dist + exec_cycles], based on the premise of earliest possible execution.
    graph_fragments$root_dist_exec_cycles <- graph_fragments$root_dist + graph_fragments$exec_cycles

    # Overlap boundaries
    overlap_interval_width <- median(graph_fragments$exec_cycles)
    stopifnot(overlap_interval_width > 0)
    overlap_bins_lower <- seq(0, max(graph_fragments$root_dist_exec_cycles) + 1 + overlap_interval_width, by=overlap_interval_width)
    stopifnot(length(overlap_bins_lower) > 1)
    overlap_bins_upper <- overlap_bins_lower + (overlap_bins_lower[2] - overlap_bins_lower[1] - 1)

    # Get overlapping ranges
    # Use IRanges package to countOverlaps using the fast NCList data structure.
    #suppressMessages(library(IRanges, quietly=TRUE, warn.conflicts=FALSE))
    #subject <- IRanges(graph_fragments$root_dist, graph_fragments$root_dist_exec_cycles)
    #query <- IRanges(overlap_bins_lower, overlap_bins_upper)
    #inst_par <- data.frame(low=overlap_bins_lower, count=countOverlaps(query, subject))

    # Get overlapping ranges
    # Use foverlaps from data.table.
    # See Arun's answer here:  http://stackoverflow.com/questions/30978837/histogram-like-summary-for-interval-data
    subject <- data.table(interval =  graph_fragments$name,
                          start = graph_fragments$root_dist,
                          end = graph_fragments$root_dist_exec_cycles)
    query <- data.table(start = overlap_bins_lower,
                        end = overlap_bins_upper)
    setkey(subject, start, end)
    overlaps <- foverlaps(query, subject, type=cl_args$overlap)
    overlaps <- overlaps[,
                         .(count = sum(!is.na(start)), fragment = paste(interval, collapse=";")),
                         by = .(i.start, i.end)]
    inst_par <- data.frame(low=overlap_bins_lower, count=overlaps$count)

    # Write instantaneous parallelism as bar graph of root distance
    temp_out_file <- paste(gsub(". $", "", cl_args$out), "-instantaneous-parallelism.pdf", sep="")
    pdf(temp_out_file)
    plot(inst_par, xlab="Distance from start node in execution cycles", ylab="Fragments", main="Instantaneous parallelism (black bars)\n Cilk average parallelism (red line)\n Number of cores (blue line)", col="black", type='h', yaxs='i')
    abline(h = length(unique(prof_data$cpu_id)), col = "blue", lty=2)
    if (!is.na(common_edge_weight[2])) {
        work <- sum(as.numeric(prof_data[,common_edge_weight[1]]))
    } else {
        work <- nrow(prof_data)*common_edge_weight[1]
    }
    abline(h = work/critical_path, col = "red", lty=1)
    res <- dev.off()
    my_print(paste("Wrote file:", temp_out_file))

    # Write instantaneous parallelism raw data to file
    temp_out_file <- paste(gsub(". $", "", cl_args$out), "-instantaneous-parallelism-raw.csv", sep="")
    res <- write.csv(overlaps, file=temp_out_file, row.names=F)
    my_print(paste("Wrote file:", temp_out_file))

    # Compute contribution of each task to instantaneous parallelism
    # Remove overlaps where count is NA
    overlaps[fragment == "NA"] <- NA
    overlaps <- overlaps[complete.cases(overlaps),]
    # Convert fragment list into rows
    overlaps <- overlaps[, list(count, fragment = as.numeric(unlist(strsplit(as.character(fragment), ';' )))), by = .(i.start, i.end)]
    # Get task of fragment
    overlaps$task <- as.integer(overlaps$fragment)
    overlaps_agg <- overlaps[, j=list(inst_par_median = as.numeric(median(count, na.rm = TRUE)), inst_par_min = as.numeric(min(count, na.rm = TRUE)), inst_par_max = as.numeric(max(count, na.rm = TRUE))), by=list(task)]
    overlaps_agg <- as.data.frame(overlaps_agg)

    # Write to graph
    for (annot in c("inst_par_median", "inst_par_min", "inst_par_max")) {
        values <- as.numeric(overlaps_agg[match(fragment_tasks, overlaps_agg$task),annot])
        grain_graph <- set.vertex.attribute(grain_graph, name=annot, index=fragment_index, value=values)
    }

    # Write per task contribution to instantaneous parallelism to file
    temp_out_file <- paste(gsub(". $", "", cl_args$out), "-instantaneous-parallelism-contribution.csv", sep="")
    res <- write.csv(overlaps_agg, file=temp_out_file, row.names=F)
    my_print(paste("Wrote file:", temp_out_file))

    if (cl_args$timing) toc("Calculating instantaneous parallelism")
}

#
# Set attributes (post graph computation)
#
if (cl_args$verbose) my_print("Setting attributes post graph computation...")
if (cl_args$timing) tic(type="elapsed")

if (!(cl_args$unreduced)) {
    # Set task grain attributes
    task_index <- match(as.character(prof_data$task), V(grain_graph)$name)

    # Set size
    if (!is.na(task_width[2])) {
        temp <- apply_task_size_mapping(as.numeric(prof_data[,task_width[1]]), task_width[2])
        grain_graph <- set.vertex.attribute(grain_graph, name=annot_name, index=task_index, value=temp)
    } else {
        grain_graph <- set.vertex.attribute(grain_graph, name="width", index=task_index, value=task_width[1])
    }
    if (!is.na(task_height[2])) {
        temp <- apply_task_size_mapping(as.numeric(prof_data[,task_height[1]]), task_height[2])
        grain_graph <- set.vertex.attribute(grain_graph, name="height", index=task_index, value=temp)
    } else {
        grain_graph <- set.vertex.attribute(grain_graph, name="height", index=task_index, value=task_height[1])
    }

    # Set color
    if (!is.na(task_color[2])) {
        temp <- apply_task_color_mapping(as.numeric(prof_data[,task_color[1]]), task_color[2], paste("task-", task_color[1], "-", task_color[2], ".colormap", sep=""))
        grain_graph <- set.vertex.attribute(grain_graph, name="color", index=task_index, value=temp)
    } else {
        grain_graph <- set.vertex.attribute(grain_graph, name="color", index=task_index, value=task_color[1])
    }

    # TODO: Map task size linearly based on "ins_count", "work_cycles", "overhead_cycles", "exec_cycles"
    # TODO: Map task color linearly based on "mem_fp", "-compute_int", "PAPI_RES_STL_sum", "-mem_hier_util", "work_deviation", "overhead_deviation", "-parallel_benefit", "-inst_par_median", "-inst_par_max","-inst_par_min", "sibling_work_balance"
    # "-" higher is better
    # TODO: Map task color linearly based on "sibling_scatter" for task-based profiling data
    # TODO: Map task color linearly based on "chunk_work_balance", "chunk_work_cpu_balance" for for-loop based profiling data
    # TODO: Map task color using linear-step mapping for "cpu_id", "outl_func", "tag", "outline_function"
}

if (cl_args$timing) toc("Setting attributes post graph computation")

#
# Compute basic information about grain graph
#
sink(grain_graph_info_out_file, append=T)
my_print("# Structure")
my_print(paste("Number of nodes =", length(V(grain_graph))))
my_print(paste("Number of edges =", length(E(grain_graph))))
my_print(paste("Number of tasks =", length(prof_data$task)))
if (cl_args$unreduced) {
    my_print(paste("Number of fragments =", length(which(graph_vertices$type == "fragment"))))
}
my_print(paste("Number of forks =", length(unique(as.character(get.vertex.attribute(grain_graph, name="name", index=fork_nodes_index))))))
my_print("# Out-degree distribution of forks")
degree.distribution(grain_graph, v=fork_nodes_index, mode="out")
if (!is.na(common_edge_weight[2])) {
    my_print(paste("# Cilk theory parallelism (metric = ", common_edge_weight[1], ")", sep=""))
    my_print(paste("Span (critical path) =", critical_path))
    work <- sum(as.numeric(prof_data[,common_edge_weight[1]]))
    my_print(paste("Work =", work))
    my_print(paste("Parallelism (Work/Span) =", work/critical_path))
} else {
    my_print(paste("# Cilk theory parallelism (metric = constant", common_edge_weight[1], ")", sep=""))
    my_print(paste("Span (critical path) =", critical_path))
    work <- nrow(prof_data)*common_edge_weight[1]
    my_print(paste("Work =", work))
    my_print(paste("Parallelism (Work/Span) =", work/critical_path))
}
if (cl_args$enumcriticalpath) {
    if (cl_args$unreduced) {
        my_print(paste("Number of critical fragments =", length(graph_vertices$name[graph_vertices$on_crit_path == 1 & graph_vertices$type == "fragment"])))
    } else {
        my_print(paste("Number of critical tasks =", length(graph_vertices$name[graph_vertices$on_crit_path == 1 & graph_vertices$type == "task"])))
    }
}
sink()

# Write information to file
if (cl_args$verbose) my_print("Writing grain graph information ...")
my_print(paste("Wrote file:", grain_graph_info_out_file))

#
# Write grain graph to file
#
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

#
# Cleanup
#

# Warn
wa <- warnings()
if (class(wa) != "NULL")
    print(wa)
