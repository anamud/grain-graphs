#
# Setup
#

# Clean slate
rm(list=ls())

# Include support functions
grain_graphs_root_dir <- Sys.getenv("GRAIN_GRAPHS_ROOT")
source(paste(grain_graphs_root_dir,"/prototype/common.R", sep=""))

# Read arguments
Rstudio_mode <- F
if (Rstudio_mode) {
  cl_args <- list(graph="grain-graph.graphml",
                 data="aggregated-grain-graph-data.csv",
                 groupname="1",
                 verbose=T,
                 timing=F,
                 out="sub-graph-<group_name>.graphml")
} else {
  option_list <- list(
    make_option(c("--graph"), default="grain-graph.graphml", help="Non-aggregated grain graph [default \"%default\"].", metavar="FILE"), # Cannot add -g as option since it is parsed by Rscript as "gui" option
    make_option(c("--data"), default="aggregated-grain-graph-data.csv", help="Aggregated grain graph data [default \"%default\"].", metavar="FILE"),
    make_option(c("--groupname"), default="1", help="Group at root of subgraph [default \"%default\"].", metavar="INT"),
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
g_data <- read.csv(cl_args$data, header=TRUE, comment.char='#', na.strings="NA")
if (!("task" %in% colnames(g_data))) {
    my_print("Error: Graph data does not have task annotation. Aborting!")
    quit("no", 1)
}
g_data[g_data == "NA"] <- NA
is.na(g_data) <- is.na(g_data)

if (cl_args$verbose) my_print("Initializing ...")
if (cl_args$timing) tic(type="elapsed")

# Function to return members of group
# Return value depends on group type
# ... "task" => Return task
# ... "family" => Return parent and all members of sibling groups
# ... "sibling" => Return siblings with fork and join nodes
# ... "non-problematic-sibling" => Return non-problematic siblings with fork and join nodes
# ... abort on unknown group type
# Function recurses into members that are groups
get_members <- function(group_name, recursive=F)
{
    #my_print(paste("Processing group", group_name, "..."))
    ioe <- g_data %>% filter(group_id == group_name) %>% select(task, group_type)

    if (nrow(ioe) == 0)
        return (group_name)

    if (recursive) {
        temp <- ioe %>% filter(group_type == "family" | group_type == "sibling" | group_type == "non-problematic-sibling")
        if (nrow(temp) > 0) {
            sub_groups <- unname(unlist(temp %>% select(task)))
            members <- unname(unlist(ioe %>% filter(task != sub_groups) %>% select(task)))
            for (sub_group in sub_groups)
            {
                members <- append(members, get_members(sub_group, recursive))
            }
        } else {
            members <- unname(unlist(ioe %>% select(task)))
        }
    } else {
        members <- unname(unlist(ioe %>% select(task)))
    }

    print(members)

    ioe <- g_data %>% filter(task == group_name) %>% select(group_type, group_leader)
    if (ioe$group_type %in% c("sibling", "non-problematic-sibling")) {
        leader <- g_data %>% filter(task == ioe$group_leader) %>% select(parent, joins_at)
        fork_node <- paste("f.", leader$parent, ".", leader$joins_at, sep="")
        join_node <- paste("j.", leader$parent, ".", leader$joins_at, sep="")
        members <- append(members, c(fork_node, join_node))
    }

    members <- unique(members)

    return (members)
}

if (cl_args$verbose) my_print(paste("Getting members of group", cl_args$groupname, "..."))
if (cl_args$timing) tic(type="elapsed")

if (cl_args$timing) toc("Getting members of group")

# Get members of group provided
# ... group exists in data
# ... else abort
if (cl_args$groupname %in% g_data$task) {
    members <- get_members(cl_args$groupname, T)
} else {
    my_print(paste("Error: Group ", group_name, " not found in graph data. Aborting!", sep=""))
    quit("no", 1)
}

if (cl_args$verbose) my_print("Plotting subgraph ...")
if (cl_args$timing) tic(type="elapsed")

# TODO: Check if mebers exist in graph

# Get subgraph
members_sg <- induced_subgraph(g_graphml, members, impl="auto")

# Plot subgraph
temp_out_file <- paste(gsub(". $", "", cl_args$out), paste("-",cl_args$groupname, ".graphml", sep=""), sep="")
res <- write.graph(members_sg, file=temp_out_file, format="graphml")
my_print(paste("Wrote file:", temp_out_file))

if (cl_args$timing) toc("Plotting subgraph")

# Warn
w <- warnings()
if (class(w) != "NULL")
  print(w)
