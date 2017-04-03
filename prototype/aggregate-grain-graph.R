#
# Setup
#

# Clean slate
rm(list=ls())

# Include support functions
mir_root <- Sys.getenv("GRAIN_GRAPHS_ROOT")
source(paste(mir_root,"/prototype/common.R",sep=""))

# Import
library(XML)

# Group size
group_width <- 50
group_height <- 50

# Read arguments
Rstudio_mode <- F
if (Rstudio_mode) {
  cl_args <- list(graph="grain-graph.graphml",
                 verbose=T,
                 timing=F,
                 outdata="aggregated-grain-graph-data.csv",
                 outgraph="aggregated-grain-graph.graphml")
} else {
  option_list <- list(
    make_option(c("--graph"), help = "Grain graph (GRAPHML).", metavar="FILE"), # Cannot add -g as option since it is parsed by Rscript as "gui" option
    make_option(c("--verbose"), action="store_true", default=TRUE, help="Print output [default]."),
    make_option(c("--timing"), action="store_true", default=FALSE, help="Print timing information."),
    make_option(c("--outdata"), default="aggregated-grain-graph-data.csv", help = "Aggregated graph data output file [default \"%default\"]", metavar="STRING"),
    make_option(c("--outgraph"), default="aggregated-grain-graph.graphml", help = "Aggregated graph output file [default \"%default\"]", metavar="STRING"),
    make_option(c("--quiet"), action="store_false", dest="verbose", help="Print little output."))

  cl_args <- parse_args(OptionParser(option_list = option_list), args = commandArgs(TRUE))

  if (!exists("graph", where=cl_args)) {
    my_print("Error: Graph argument missing. Check help (-h).")
    quit("no", 1)
  }
}

if (cl_args$verbose) my_print(paste("Reading", cl_args$graph, "..."))

# Read data from graph
g_graphml <- read.graph(cl_args$graph, format="graphml")
g_data <- get.data.frame(g_graphml, what="vertices")
if (!("task" %in% colnames(g_data))) {
    my_print("Error: Graph does not have task annotation. Aborting!")
    quit("no", 1)
}
g_data[g_data == "NA"] <- NA
is.na(g_data) <- is.na(g_data)
g_data <- subset(g_data, !is.na(task))

# Read graph as XML
g <-xmlParse(cl_args$graph)

# Graph variables
# Namespace
ns <- c(graphml = "http://graphml.graphdrawing.org/xmlns")
# Key IDs for vertices
key_ids <- sapply(getNodeSet(g,"//graphml:key", namespaces = ns), xmlGetAttr, "id")
key_ids <- key_ids[grep("^v_", key_ids)]
stopifnot(length(key_ids) > 0)
# # Initializer for key values
# #... DO NOT USE.
# #... yEd complains about invalid format
# key_init <- function(attr_val) { n_n <- xmlNode("data", attrs = c("key" = attr_val), "NA"); n_n}
# # Index for keys of interest
# v_label_ind <- which(key_ids == "v_label")
# v_name_ind <- which(key_ids == "v_name")
# v_task_ind <- which(key_ids == "v_task")
# v_parent_ind <- which(key_ids == "v_parent")
# v_joins_at_ind <- which(key_ids == "v_joins_at")
# Edges
edges_xml <- getNodeSet(g, "//graphml:edge", namespaces = ns)
# Nodes
nodes_xml <- getNodeSet(g, "//graphml:node", namespaces = ns)
stopifnot(xmlSize(nodes_xml) > 0)
node_ids <- sapply(nodes_xml, xmlGetAttr, "id")
node_names <- sapply(getNodeSet(g, "//graphml:node/graphml:data[@key='v_name']", namespaces = ns),
                     xmlValue)
# Node table
node_ids_names <- data.frame(id = node_ids, name = node_names)
# Groups
# Root node to hold groups
groups_xml <-xmlNode("root")
# Group table
group_ids_names <- data.frame(id = character(0), name = character(0))

# Counter to create unique group IDs
# Starts at max_task_id + 1
# Note: <<- is a scoping assignment allowing the function to maintain state
max_task_id <- max(g_data$task)
get_group_id = (function(){gid = max_task_id; function() gid <<- gid + 1 })()

# Set group variables
# Unique ID of the group
# Each task belongs to its own group inititally
# During aggregation, tasks are nested under explicitly created groups
# A group is also added as a task
g_data$group_id <- g_data$task
# Indicates if the task is a member of (i.e., nested immediately under) an explicit group
g_data$grouped <- F
# The number of task nested under the group including sub-groups
g_data$num_tasks <- 1
# The number of tasks nested *immediately* under the group
g_data$num_members <- 1
# Group types: task, sibling, family
# Task means an implicit group
g_data$group_type <- "task"

# Rough sketch of iterative algorithm for grouping
# 1. Mark all leaf siblings
# ... Leaf siblings are tasks without children and sharing the same fork/join points
# ... Single leaf siblings are possible too
# ... Since groups are added as tasks, groups that share the same join point are leaf siblings too
# 2. Group leaf siblings, join point, fork point together (members) into a leaf sibling group
# 3. Assign unique id, accumulate attributes and performance of group members sensibly, and add as leaf task
# 4. Mark all families
# ... Families are headed by parents whose children completely belong to leaf sibling groups
# 5. Group families i.e., parent and leaf sibling groups of children together (members) into a family group
# 6. Assign unique ID, accumulate attributes and performance of group members sensibly, and add as leaf task
# 7. Repeat 1-6 until all tasks except the first implicit task are grouped

if (cl_args$verbose) my_print("Aggregating ...")

itr_count <- 0
while(any(!g_data$grouped))
{
  itr_count  <- itr_count + 1
  if (cl_args$verbose) my_print(paste("In grouping iteration", itr_count))

  # Break out if only first implicit task is ungrouped
  ungrouped <- which(!g_data$grouped)
  if(length(ungrouped) == 1)
  {
    if(g_data[ungrouped,]$parent == 0)
    {
      if (cl_args$verbose) my_print("Only implicit task is ungrouped. Stopping.")
      break
    }
  }

  # Group leaf siblings
  #if (cl_args$verbose) my_print("Grouping leaf siblings ...")
  if (cl_args$timing) tic(type="elapsed")

  # Mark leaf siblings
  e0 <- g_data %>% group_by(parent, joins_at) %>% filter(leaf == T & grouped == F)
  # Group and compute aggregated attributes
  #... Group has same parent and join parent as members
  #... Assiging unique ID row-wise is essential else we obtain inconsistent results
  #... TODO: Understand why row-wise is required to ensure consistency
  e <- e0 %>% summarize(num_tasks = sum(num_tasks), num_members = n(), work_cycles = sum(as.numeric(work_cycles)), child_number = NA, leaf = T, group_id = NA, grouped = F, group_type = "sibling") %>% rowwise() %>% mutate(task = get_group_id())
  # Ensure groups exist
  stopifnot(nrow(e) > 0)

  # For each group,
  #... compute aggregated attributes
  #... mark members as grouped and update group_id
  #... FIXME: Loop is slow. Make fast!
  for(i in seq_len(nrow(e)))
  {
    # For quick retrieval
    ei_task <-  e[i,]$task
    ei_parent <- e[i,]$parent
    ei_joins_at <- e[i,]$joins_at

    # Update members attributes
    #e1 <- e0 %>% filter(parent == e[i,]$'parent' & joins_at == e[i,]$'joins_at') %>% select(task, parent, joins_at)
    e1 <- subset(e0, parent == ei_parent & joins_at == ei_joins_at)
    matches <- which(g_data$task %in% e1$task)
    g_data[matches, ]$group_id <- ei_task
    g_data[matches, ]$grouped <- T

    # Update group attributes
    e[i,]$group_id <- ei_task
    # Save parent's child count for use during family grouping
    match <- which(g_data$task == ei_parent)
    stopifnot(length(match) == 1)
    # Using child_number as a placeholder for number of children of parent
    e[i,]$child_number <- g_data[match,]$num_children

    # Add group node to graph
    # Create new group node
    new_group <- xmlNode("graph", attrs = c("id" = paste("g", as.character(ei_task), sep = ""),
                                            "edgedefault" = "directed"))
    # Add nodes of members as children
    parent_join <- paste(as.character(ei_parent), ".", as.character(ei_joins_at), sep = "")
    fork_join_nodes <- c(paste("f.", parent_join, sep = ""), paste("j.", parent_join, sep = ""))
    member_nodes <- which(node_ids_names$name %in% c(as.character(e1$task), fork_join_nodes))
    new_group <- append.xmlNode(new_group, lapply(member_nodes, function(member_node) nodes_xml[[member_node]]))
    member_groups <- which(group_ids_names$name %in% as.character(e1$task))
    new_group <- append.xmlNode(new_group, lapply(member_groups, function(member_group) groups_xml[[member_group]]))
    # Node wrapper
    node_wrapper <- xmlNode("node", attrs = c("id" = paste("g", as.character(ei_task), "n", sep = ""), "yfiles.foldertype" = "folder"))
    # Add specific key nodes as children
    node_wrapper <- append.xmlNode(node_wrapper,
                                  xmlNode("data", attrs = c("key" = "v_label"), as.character(ei_task)),
                                  xmlNode("data", attrs = c("key" = "v_name"), as.character(ei_task)),
                                  xmlNode("data", attrs = c("key" = "v_task"), as.character(ei_task)),
                                  xmlNode("data", attrs = c("key" = "v_parent"), as.character(ei_parent)),
                                  xmlNode("data", attrs = c("key" = "v_joins_at"), as.character(ei_joins_at)),
                                  xmlNode("data", attrs = c("key" = "v_work_cycles"), as.character(e[i,]$work_cycles)),
                                  xmlNode("data", attrs = c("key" = "v_width"), as.character(group_width)),
                                  xmlNode("data", attrs = c("key" = "v_height"), as.character(group_height)),
                                  xmlNode("data", attrs = c("key" = "v_shape"), "round rectangle"),
                                  xmlNode("data", attrs = c("key" = "v_num_tasks"), as.character(e[i,]$num_tasks)),
                                  xmlNode("data", attrs = c("key" = "v_num_members"), as.character(e[i,]$num_members)),
                                  xmlNode("data", attrs = c("key" = "v_group_type"), as.character(e[i,]$group_type))
                                  )
    new_group <- append.xmlNode(node_wrapper, new_group)
    # Keep only yet to be grouped
    if(length(member_groups) > 0)
    {
        group_ids_names <- group_ids_names[-c(member_groups),]
        groups_xml <- removeChildren(groups_xml, kids = as.numeric(member_groups), free = T)
    }
    # Add group as child
    groups_xml <- append.xmlNode(groups_xml, new_group)
    # Update table
    group_ids_names <- bind_rows(group_ids_names, data.frame(id = paste("g", as.character(ei_task), "n", sep = ""), name = as.character(ei_task)))
    # Make sure group table mirrors order in XML root
    # TODO: Make fast
    stopifnot(all(group_ids_names$id == as.character(lapply(xmlChildren(groups_xml), xmlGetAttr, "id"))))

    # Release unused nodes
    rm(new_group)
    rm(node_wrapper)
  }

  # Add leaf sibling groups as tasks
  g_data <- bind_rows(g_data, e)

  # Garbage collect
  invisible(gc(reset=T))

  if (cl_args$timing) toc("Group leaf siblings")

  # Group families
  #if (cl_args$verbose) my_print("Grouping families ...")
  if (cl_args$timing) tic(type="elapsed")

  # Mark groups that completely contain children
  #... TODO: Pick the first instead of maximum child_number
  #... ... All child_numbers are the same since they are proxies for child count of parent
  f <- g_data %>% group_by(parent) %>% filter(task > max_task_id & grouped == F) %>% summarize(num_tasks = sum(num_tasks), num_members = sum(num_members), num_sibling_groups = n(), work_cycles = sum(as.numeric(work_cycles)), child_number = max(child_number), leaf = T, group_id = NA, grouped = F) %>% filter(num_members == child_number)
  stopifnot(nrow(f) > 0)
  f$task <- NA
  f$joins_at <- NA
  f$group_type <- "family"

  # For each group,
  #... compute aggregated attributes
  #... mark members as grouped and update group_id
  #... FIXME: Loop is slow. Make fast!
  for(i in seq_len(nrow(f)))
  {
    # Update group attributes
    f[i,]$task <- get_group_id()
    f[i,]$num_members <- f[i,]$num_sibling_groups + 1
    f[i,]$group_id <- f[i,]$task

    # For quick retrieval
    fi_task <-  f[i,]$task

    # Update member attributes
    # Update parent member attributes
    match <- which(g_data$task == f[i,]$parent)
    stopifnot(length(match) == 1)

    g_data[match, ]$group_id <- fi_task
    g_data[match, ]$grouped <- T

    # Update group attributes
    f[i,]$joins_at <- g_data[match,]$joins_at
    temp <- f[i,]$parent
    f[i,]$parent <- g_data[match,]$parent
    f[i,]$work_cycles <- as.numeric(f[i,]$work_cycles) + as.numeric(g_data[match,]$work_cycles)
    f[i,]$child_number <- g_data[match,]$num_children
    f[i,]$num_tasks <-  f[i,]$num_tasks + g_data[match,]$num_tasks

    # Update child member attributes
    matches <- which(g_data$parent == temp & g_data$task > max_task_id & g_data$grouped == F)
    stopifnot(length(matches) > 0)

    g_data[matches, ]$group_id <- fi_task
    g_data[matches, ]$grouped <- T

    # Add group node to graph
    # Create new group node
    new_group <- xmlNode("graph", attrs = c("id" = paste("g", as.character(fi_task), sep = ""),
                                            "edgedefault" = "directed"))
    # Add nodes of members as children
    children <- as.character(c(g_data[matches,]$task, g_data[match,]$task))
    member_nodes <- which(node_ids_names$name %in% children)
    new_group <- append.xmlNode(new_group, lapply(member_nodes, function(member_node) nodes_xml[[member_node]]))
    member_groups <- which(group_ids_names$name %in% children)
    new_group <- append.xmlNode(new_group, lapply(member_groups, function(member_group) groups_xml[[member_group]]))
    # Node wrapper
    node_wrapper <- xmlNode("node", attrs = c("id" = paste("g", as.character(fi_task), "n", sep = ""), "yfiles.foldertype" = "folder"))
    # Add specific key nodes as children
    node_wrapper <- append.xmlNode(node_wrapper,
                                   xmlNode("data", attrs = c("key" = "v_label"), as.character(fi_task)),
                                   xmlNode("data", attrs = c("key" = "v_name"), as.character(fi_task)),
                                   xmlNode("data", attrs = c("key" = "v_task"), as.character(fi_task)),
                                   xmlNode("data", attrs = c("key" = "v_parent"), as.character(f[i,]$parent)),
                                   xmlNode("data", attrs = c("key" = "v_joins_at"), as.character(f[i,]$joins_at)),
                                   xmlNode("data", attrs = c("key" = "v_work_cycles"), as.character(f[i,]$work_cycles)),
                                   xmlNode("data", attrs = c("key" = "v_size"), as.character(group_size)),
                                   xmlNode("data", attrs = c("key" = "v_width"), as.character(group_size)),
                                   xmlNode("data", attrs = c("key" = "v_height"), as.character(group_size)),
                                   xmlNode("data", attrs = c("key" = "v_shape"), "round rectangle"),
                                   xmlNode("data", attrs = c("key" = "v_num_tasks"), as.character(f[i,]$num_tasks)),
                                   xmlNode("data", attrs = c("key" = "v_num_members"), as.character(f[i,]$num_members)),
                                   xmlNode("data", attrs = c("key" = "v_group_size"), as.character(group_size)),
                                   xmlNode("data", attrs = c("key" = "v_group_type"), as.character(f[i,]$group_type))
                                   )
    new_group <- append.xmlNode(node_wrapper, new_group)
    # Keep only yet to be grouped
    if(length(member_groups) > 0)
    {
        group_ids_names <- group_ids_names[-c(member_groups),]
        groups_xml <- removeChildren(groups_xml, kids = as.numeric(member_groups), free = T)
    }
    # Add group as child
    groups_xml <- append.xmlNode(groups_xml, new_group)
    # Update table
    group_ids_names <- bind_rows(group_ids_names, data.frame(id = paste("g", as.character(fi_task), "n", sep = ""), name = as.character(fi_task)))
    # Make sure group table mirrors order in XML root
    # TODO: Make fast
    stopifnot(all(group_ids_names$id == as.character(lapply(xmlChildren(groups_xml), xmlGetAttr, "id"))))

    # Release unused nodes
    rm(new_group)
    rm(node_wrapper)
  }

  # Add families as tasks
  g_data <- bind_rows(g_data, f %>% select(-num_sibling_groups))

  # Garbage collect
  invisible(gc(reset=T))

  if (cl_args$timing) toc("Group families")
}

#if (cl_args$verbose) my_print("Post-grouping processes ...")
if (cl_args$timing) tic(type="elapsed")

# Explicitly mark groups for covnenience
g_data$group <- ifelse(g_data$task > max_task_id, T, F)

# Add graphml wrappers to graph
graph_wrapper <- xmlNode("graphml", attrs = c("xmlns" = "http://graphml.graphdrawing.org/xmlns",
                                              "xmlns:xsi" = "http://www.w3.org/2001/XMLSchema-instance",
                                              "xsi:schemaLocation" = "http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd"
                                              ))
graph_wrapper <- append.xmlNode(graph_wrapper, xmlCommentNode("Created by MIR profiling infrastructure"))
# Key id information
graph_wrapper <- append.xmlNode(graph_wrapper, getNodeSet(g, "//graphml:key", namespaces = ns))
# Group specific keys
kn_num_tasks <- xmlNode("key", attrs = c("id" = "v_num_tasks",
                         "for" = "node",
                         "attr.name" = "num_tasks",
                         "attr.type" = "double"))
graph_wrapper <- append.xmlNode(graph_wrapper, kn_num_tasks)
kn_num_members <- xmlNode("key", attrs = c("id" = "v_num_members",
                         "for" = "node",
                         "attr.name" = "num_members",
                         "attr.type" = "double"))
graph_wrapper <- append.xmlNode(graph_wrapper, kn_num_members)
kn_group_size <- xmlNode("key", attrs = c("id" = "v_group_size",
                         "for" = "node",
                         "attr.name" = "group_size",
                         "attr.type" = "double"))
graph_wrapper <- append.xmlNode(graph_wrapper, kn_group_size)
kn_group_type <- xmlNode("key", attrs = c("id" = "v_group_type",
                         "for" = "node",
                         "attr.name" = "group_type",
                         "attr.type" = "string"))
graph_wrapper <- append.xmlNode(graph_wrapper, kn_group_type)
# Set top graph
top_graph <- xmlNode("graph", attrs = c("id" = "G", "edgedefault" = "directed"))
# Special nodes S, E, f.0.0, and j.0.0.0
# S has name = 0
special_nodes <- which(node_ids_names$name %in% c("0", "E", "f.0.0", "j.0.0"))
top_graph <- append.xmlNode(top_graph, lapply(special_nodes, function(special_node) nodes_xml[[special_node]]))
# The last added group node contains the full graph
top_graph <- append.xmlNode(top_graph, groups_xml[[nrow(group_ids_names)]])
# Edges
top_graph <- append.xmlNode(top_graph, edges_xml)
# Add top graph
graph_wrapper <- append.xmlNode(graph_wrapper, top_graph)
# Release unused nodes and garbage collect
rm(top_graph)
invisible(gc(reset=T))

if (cl_args$timing) toc("Post-grouping processes")

# Write out aggregated data
if (cl_args$verbose) my_print("Writing aggregated data ...")
if (cl_args$timing) tic(type="elapsed")

out_file <- cl_args$outdata
sink(out_file)
write.csv(g_data, out_file, row.names=F)
sink()
my_print(paste("Wrote file:", out_file))

out_file <- cl_args$outgraph
junk <- saveXML(graph_wrapper, out_file)
my_print(paste("Wrote file:", out_file))

if (cl_args$timing) toc("Write aggregated data")

# Warn
w <- warnings()
if (class(w) != "NULL")
  print(w)
