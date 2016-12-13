# Clean slate
rm(list=ls())

# Treat warnings as errors
#options(warn=2)

# Strings as factors
options(stringsAsFactors = F)

# Include
mir_root <- Sys.getenv("MIR_ROOT")
source(paste(mir_root,"/scripts/profiling/task/common.R",sep=""))

# Library
suppressMessages(library(data.table, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(library(dplyr))
suppressMessages(library(igraph, quietly=TRUE))
#library(bit64)

# Parse arguments
# TODO: Understand how to capture if not running inside RStudio.
running_outside_rstudio <- 1
if (running_outside_rstudio) {
    library(optparse, quietly=TRUE)
    option_list <- list(
                        make_option(c("-d","--data"), help = "Task performance data file.", metavar="FILE"),
                        make_option(c("-p","--palette"), default="color", help = "Color palette for graph elements [default \"%default\"]."),
                        make_option(c("-o","--out"), default="full-task-graph", help = "Output file prefix [default \"%default\"].", metavar="STRING"),
                        make_option(c("--verbose"), action="store_true", default=TRUE, help="Print output [default]."),
                        make_option(c("--quiet"), action="store_false", dest="verbose", help="Print little output."),
                        make_option(c("--timing"), action="store_true", default=FALSE, help="Print processing time."),
                        make_option(c("--overlap"), default="any", help = "Overlap type for instantaneous parallelism calculation. Choose one among: any, full. [default \"%default\"].", metavar="STRING"))
    parsed <- parse_args(OptionParser(option_list = option_list), args = commandArgs(TRUE))
    if (!exists("data", where=parsed))
    {
        my_print("Error: Data argument missing. Check help (-h)")
        quit("no", 1)
    }
    if (!(parsed$overlap == "any" | parsed$overlap == "full"))
    {
        my_print("Error: Invalid overlap argument. Check help (-h)")
        quit("no", 1)
    }
    if (parsed$overlap == "full")
    {
        parsed$overlap <- "within"
    }

    # Set argments into placeholders
    # Placeholders help while testing in RStudio where command line arguments are difficult to pass.

    arg_data <- parsed$data
    arg_palette <- parsed$pal
    arg_outfileprefix <- parsed$out
    arg_verbose <- parsed$verbose
    arg_timing <- parsed$timing
} else {
    arg_data <- "mir-task-stats"
    arg_palette <- "color"
    arg_overlap <- "any"
    arg_outfileprefix <- "full-task-graph"
    arg_verbose <- 1
    arg_timing <- 1
}

# Read data
if (arg_timing) tic(type="elapsed")
tg_file_in <- arg_data
if (arg_verbose) my_print(paste("Reading file:", tg_file_in, sep=" "))
tg_data <- fread(tg_file_in, header=TRUE)
if ("integer64" %in% sapply(tg_data, class))
{
    my_print("Error: Data contains unsupported integer64 class.")
    quit("no", 1)
}
if (arg_timing) toc("Read data")

# Join frequeny
if (arg_timing) tic(type="elapsed")
join_freq <- tg_data %>% arrange(parent, joins_at) %>% group_by(parent, joins_at) %>% summarise(count = n())
if (arg_timing) toc("Compute join frequency")

# Funciton returns edges of fragment chain for input task
fragmentize <- function (task, num_children, parent, child_number, joins_at)
{
    num_fragments <- 1
    # Connect to parent fork
    new_fragment <- paste(task, num_fragments, sep='.')
    parent_fork <- paste('f', parent, child_number, sep='.')
    edg <- c(parent_fork, new_fragment)
    last_fragment <- new_fragment
    # Connect fragments
    if (num_children > 0)
    {
        joins <- join_freq[parent == task]$count
        joins_ind <- 1
        num_forks <- 1
        num_joins <- 0
        while(1)
        {
            stopifnot(joins[joins_ind] > 0)
            for(i in 1:joins[joins_ind])
            {
                num_fragments <- num_fragments + 1
                new_fragment <- paste(task, num_fragments, sep='.')
                fork <- paste('f', task, num_forks, sep=".")
                edg <- c(edg, last_fragment, fork, fork, new_fragment)
                num_forks <- num_forks + 1
                last_fragment <- new_fragment
                if (i == joins[joins_ind])
                {
                    num_fragments <- num_fragments + 1
                    new_fragment <- paste(task, num_fragments, sep='.')
                    join <- paste('j', task, num_joins, sep=".")
                    edg <- c(edg, last_fragment, join, join, new_fragment)
                    num_joins <- num_joins + 1
                    last_fragment <- new_fragment
                }
            }
            joins_ind <- joins_ind + 1
            if (joins_ind > length(joins))
                break
        }
    }
    # Connect to parent join
    parent_join <- paste('j', parent, joins_at, sep='.')
    edg <- c(edg, last_fragment, parent_join)
    return(edg)
}

# Apply fragmentize on all tasks
if (arg_timing) tic(type="elapsed")
tg_edges <- unlist(mapply(fragmentize,
                          task=tg_data$task,
                          num_children=tg_data$num_children,
                          parent=tg_data$parent,
                          child_number=tg_data$child_number,
                          joins_at=tg_data$joins_at))
tg_edges <- matrix(tg_edges, nc=2, byrow=TRUE)
if (arg_timing) toc("Fragmentize")

# Construct graph
if (arg_timing) tic(type="elapsed")
tg <- graph.edgelist(tg_edges, directed = TRUE)
if (arg_timing) toc("Construct graph")

# Visual attributes
node_min_size <- 10
node_max_size <- 50

# Assign weights to nodes
# Read graph as table
if (arg_timing) tic(type="elapsed")
tg_vertices <- data.table(node=get.data.frame(tg, what="vertices")$name, exec_cycles=0, kind="fragment")
pesky_factors <- sapply(tg_vertices, is.factor)
tg_vertices[pesky_factors] <- lapply(tg_vertices[pesky_factors], as.character)
# Get helpful indexes
fork_ind <- with(tg_vertices, grepl("^f.[0-9]+.[0-9]+$", node))
join_ind <- with(tg_vertices, grepl("^j.[0-9]+.[0-9]+$", node))
fragment_ind <- with(tg_vertices, grepl("^[0-9]+.[0-9]+$", node))
start_ind <- with(tg_vertices, grepl("^f.0.1$", node))
end_ind <- with(tg_vertices, grepl("^j.0.0$", node))
# Assign kind
tg_vertices[fork_ind]$kind <- "fork"
tg_vertices[join_ind]$kind <- "join"
tg <- set.vertex.attribute(tg, name="kind", index=V(tg), value=tg_vertices$kind)
if (arg_timing) toc("Assign other weights")
# Assign exec cycles
compute_fragment_duration <- function(task, wait, exec_cycles, choice)
{
    # Each task has breaks at these instants: (execution start, child creation, child wait, execution end)
    # A fragments executes upto the next break.
    wait_instants <- as.numeric(unlist(strsplit(substring(wait, 2, nchar(wait)-1), ";", fixed = TRUE)))
    create_instants <- as.numeric(tg_data$create_instant[tg_data$parent == task])
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
if (arg_timing) tic(type="elapsed")
fd <- data.table(fragment=unlist(mapply(compute_fragment_duration,
                                        task=tg_data$task,
                                        wait=tg_data$wait_instants,
                                        exec_cycles=tg_data$exec_cycles,
                                        choice=1)),
                 duration=unlist(mapply(compute_fragment_duration,
                                        task=tg_data$task,
                                        wait=tg_data$wait_instants,
                                        exec_cycles=tg_data$exec_cycles,
                                        choice=2)))
if (arg_timing) toc("Assign execution cycles [step 1]")
if (arg_timing) tic(type="elapsed")
tg_vertices[match(fd$fragment, tg_vertices$node)]$exec_cycles <- fd$duration
if (arg_timing) toc("Assign execution cycles [step 2]")
if (arg_timing) tic(type="elapsed")
tg <- set.vertex.attribute(tg, name="exec_cycles", index=V(tg), value=tg_vertices$exec_cycles)
scale_fn <- function(x) { node_min_size + (x * 100/max(x, na.rm = TRUE)) }
tg <- set.vertex.attribute(tg, name="scaled_exec_cycles", index=V(tg), value=scale_fn(tg_vertices$exec_cycles))
if (arg_timing) toc("Assign execution cycles [step 3]")

# Calculate critical path
if (arg_verbose) my_print("Calculating critical path ...")
if (arg_timing) tic(type="elapsed")
#Rprof("profile-critpathcalc.out")
# Progress bar
lntg <- length(V(tg))
if (arg_verbose) {
    pb <- txtProgressBar(min = 0, max = lntg, style = 3)
    ctr <- 0
}
# Topological sort
tsg <- topological.sort(tg)
# Set root path attributes
V(tg)[tsg[1]]$rdist <- 0
V(tg)[tsg[1]]$depth <- 0
V(tg)[tsg[1]]$rpath <- tsg[1]
# Get data frame of task graph.
# It is much faster to work on data frame than the task graph.
# TODO: Convert to data.table for speed.
# TODO: Write loop in C for speed.
tg_vertices_df <- get.data.frame(tg, what="vertices")
# Get longest paths from root
for(node in tsg[-1])
{
    # Get distance from node's predecessors
    ni <- incident(tg, node, mode="in")
    w <- V(tg)[get.edges(tg, ni)[,1]]$exec_cycles
    # Get distance from root to node's predecessors
    nn <- neighbors(tg, node, mode="in")
    d <- tg_vertices_df$rdist[nn]
    # Add distances (assuming one-one corr.)
    wd <- w+d
    # Set node's distance from root to max of added distances
    mwd <- max(wd)
    tg_vertices_df$rdist[node] <- mwd
    # Set node's path from root to path of max of added distances
    mwdn <- as.vector(nn)[match(mwd,wd)]
    nrp <- list(c(unlist(tg_vertices_df$rpath[mwdn]), node))
    tg_vertices_df$rpath[node] <- nrp
    # Set node's depth as one greater than the largest depth its predecessors
    tg_vertices_df$depth[node] <- max(tg_vertices_df$depth[nn]) + 1
    # Progress report
    if (arg_verbose) {
        ctr <- ctr + 1;
        setTxtProgressBar(pb, ctr);
    }
}
# Longest path is the largest root distance
lpl <- max(tg_vertices_df$rdist)
# Enumerate longest path
lpm <- unlist(tg_vertices_df$rpath[match(lpl,tg_vertices_df$rdist)])
tg_vertices_df$on_crit_path <- 0
tg_vertices_df$on_crit_path[lpm] <- 1
# Assign to task graph
tg <- set.vertex.attribute(tg, name="on_crit_path", index=V(tg), value=tg_vertices_df$on_crit_path)
tg <- set.vertex.attribute(tg, name="rdist", index=V(tg), value=tg_vertices_df$rdist)
tg <- set.vertex.attribute(tg, name="depth", index=V(tg), value=tg_vertices_df$depth)
critical_edges <- E(tg)[V(tg)[on_crit_path==1] %--% V(tg)[on_crit_path==1]]
tg <- set.edge.attribute(tg, name="on_crit_path", index=critical_edges, value=1)
# Progress report
if (arg_verbose) {
    ctr <- ctr + 1;
    setTxtProgressBar(pb, ctr);
    close(pb)
}
#Rprof(NULL)
# Print critical path info
tg_file_out <- paste(gsub(". $", "", arg_outfileprefix), ".info", sep="")
sink(tg_file_out)
my_print("Cilk Theory Parallelism (Unit = Cycles):")
my_print("Span (critical path):")
my_print(lpl)
my_print("Work:")
total_work <- sum(as.numeric(tg_data$exec_cycles))
my_print(total_work)
my_print("Parallelism (Work/Span):")
parallelism <- total_work/lpl
my_print(parallelism)
sink()
my_print(paste("Wrote file:", tg_file_out))
# Clear rpath since dot/table writing complains
tg <- remove.vertex.attribute(tg,"rpath")
if (arg_timing) toc("Critical path calculation")

# Calc shape
if (arg_timing) tic(type="elapsed")
tg_vertices_df <- get.data.frame(tg, what="vertices")
# Select only fragments
tg_vertices_df <- tg_vertices_df[with(tg_vertices_df, grepl("^[0-9]+.[0-9]+$", name)),]
# Shape is found by counting overlapping ranges.
# See question http://stackoverflow.com/questions/30978837/histogram-like-summary-for-interval-data
# Each fragment has exection range [rdist, rdist + exec_cycles], based on the premise of earliest possible execution.
tg_vertices_df$rdist_exec_cycles <- tg_vertices_df$rdist + tg_vertices_df$exec_cycles
if (arg_timing) toc("Shape calculation [Step 1]")
if (arg_timing) tic(type="elapsed")
# Calculate breaks based on summary lengths of fragment
#shape_breaks <- total_work/(length(unique(tg_data$cpu_id))*median(tg_vertices_df$exec_cycles))
shape_interval_width <- median(tg_vertices_df$exec_cycles)
stopifnot(shape_interval_width > 0)
shape_bins_lower <- seq(0, max(tg_vertices_df$rdist_exec_cycles) + 1 + shape_interval_width, by=shape_interval_width)
stopifnot(length(shape_bins_lower) > 1)
shape_bins_upper <- shape_bins_lower + (shape_bins_lower[2] - shape_bins_lower[1] - 1)
if (arg_timing) toc("Shape calculation [Step 2]")
if (arg_timing) tic(type="elapsed")
# Use IRanges package to countOverlaps using the fast NCList data structure.
#suppressMessages(library(IRanges, quietly=TRUE, warn.conflicts=FALSE))
#subject <- IRanges(tg_vertices_df$rdist, tg_vertices_df$rdist_exec_cycles)
#query <- IRanges(shape_bins_lower, shape_bins_upper)
# Use foverlaps from data.table.
# See Arun's answer here:  http://stackoverflow.com/questions/30978837/histogram-like-summary-for-interval-data
subject <- data.table(interval =  tg_vertices_df$name,
                      start = tg_vertices_df$rdist,
                      end = tg_vertices_df$rdist_exec_cycles)
query <- data.table(start = shape_bins_lower,
                    end = shape_bins_upper)
setkey(subject, start, end)
overlaps <- foverlaps(query, subject, type=parsed$overlap)
overlaps <- overlaps[,
                     .(count = sum(!is.na(start)), fragment = paste(interval, collapse=";")),
                     by = .(i.start, i.end)]
if (arg_timing) toc("Shape calculation [Step 3]")
if (arg_timing) tic(type="elapsed")
#tg_shape <- data.frame(low=shape_bins_lower, count=countOverlaps(query, subject))
tg_shape <- data.frame(low=shape_bins_lower, count=overlaps$count)
if (arg_timing) toc("Shape calculation [Step 4]")
# Write shape to file
tg_file_out <- paste(gsub(". $", "", arg_outfileprefix), "-shape.pdf", sep="")
pdf(tg_file_out)
plot(tg_shape, xlab="Distance from START in execution cycles", ylab="Fragments", main="Instantaneous task parallelism", col="black", type='h', yaxs='i')
abline(h = length(unique(tg_data$cpu_id)), col = "blue", lty=2)
abline(h = parallelism , col = "red", lty=1)
write_res <- dev.off()
my_print(paste("Wrote file:", tg_file_out))
# Write shape raw data to file
tg_file_out <- paste(gsub(". $", "", arg_outfileprefix), "-shape.csv", sep="")
write_res <- write.csv(overlaps, file=tg_file_out, row.names=F)
my_print(paste("Wrote file:", tg_file_out))

# Task shape contribution
if (arg_timing) tic(type="elapsed")
# Remove intervals where count is NA
overlaps[fragment == "NA"] <- NA
overlaps <- overlaps[complete.cases(overlaps),]
# Melt csv fragment list into rows
overlaps <- overlaps[, list(count, fragment = as.numeric(unlist(strsplit(as.character(fragment), ';' )))), by = .(i.start, i.end)]
# Get task of fragment
overlaps$task <- as.integer(overlaps$fragment)
overlaps_agg <- overlaps[, j=list(median_shape_contrib = as.numeric(median(count, na.rm = TRUE)), min_shape_contrib = as.numeric(min(count, na.rm = TRUE)), max_shape_contrib = as.numeric(max(count, na.rm = TRUE))), by=list(task)]
if (arg_timing) toc("Shape calculation [Step 5]")
# Write shape contribution to file
tg_file_out <- paste(gsub(". $", "", arg_outfileprefix), "-shape-contrib.csv", sep="")
write_res <- write.csv(overlaps_agg, file=tg_file_out, row.names=F)
my_print(paste("Wrote file:", tg_file_out))

# Write task graph to file
tg_file_out <- paste(gsub(". $", "", arg_outfileprefix), ".graphml", sep="")
write_res <- write.graph(tg, file=tg_file_out, format="graphml")
my_print(paste("Wrote file:", tg_file_out))

# Warn
wa <- warnings()
if (class(wa) != "NULL")
    print(wa)
