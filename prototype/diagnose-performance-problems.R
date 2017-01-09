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
    cl_args <- list(graph="grain-graph.graphml",
                   problemthresholdconfig="problem-thresholds.cfg",
                   out="grain-graph-view",
                   layout=F,
                   verbose=T,
                   timing=F)
} else {
    option_list <- list(make_option(c("--graph"), help = "Grain graph in GRAPHML format.", metavar="FILE"),
                        make_option(c("--grainproblemconfig"), default="grain-problems.cfg", help = "Grain problem configuration file [default \"%default\"].", metavar="FILE"),
                        make_option(c("-o","--out"), default="grain-graph-diagnosis", help = "Output file suffix [default \"%default\"].", metavar="STRING"),
                        make_option(c("--layout"), action="store_true", default=FALSE, help="Layout using Sugiyama style and plot to PDF."),
                        make_option(c("--verbose"), action="store_true", default=TRUE, help="Print output [default]."),
                        make_option(c("--quiet"), action="store_false", dest="verbose", help="Print little output."),
                        make_option(c("--timing"), action="store_true", default=FALSE, help="Print processing time."))

    cl_args <- parse_args(OptionParser(option_list = option_list), args = commandArgs(TRUE))

    if (!exists("graph", where=cl_args)) {
        my_print("Error: Graph argument missing. Check help (-h)")
        quit("no", 1)
    }
}

#
# Initialize variables
#
if (cl_args$verbose) my_print("Initializing ...")
if (cl_args$timing) tic(type="elapsed")

# Prepare information output file
graph_diagnosis_info_out_file <- paste(gsub(". $", "", cl_args$out), ".info", sep="")
sink(graph_diagnosis_info_out_file)
my_print("# Problems")
sink()

# Read graph
grain_graph <- read.graph(cl_args$graph, format="graphml")
task_data <- get.data.frame(grain_graph, what="vertices")
if (!("task" %in% colnames(task_data))) {
    my_print("Error: Graph does not have task annotation. Aborting!")
    quit("no", 1)
}
task_data[task_data == "NA"] <- NA
is.na(task_data) <- is.na(task_data)
task_data <- subset(task_data, !is.na(task))

# Dim elements
dimming_alpha <- 0.2
vertex_colors_dimmed <- add.alpha(get.vertex.attribute(grain_graph, name='color'), alpha=dimming_alpha)
grain_graph <- set.vertex.attribute(grain_graph, name='color', value=vertex_colors_dimmed)
edge_colors_dimmed <- add.alpha(get.edge.attribute(grain_graph, name='color'), alpha=dimming_alpha)
grain_graph <- set.edge.attribute(grain_graph, name='color', value=edge_colors_dimmed)
#grain_graph <- set.edge.attribute(grain_graph, name='color', value="#c0c0c0")

# Mark vertices as non-problematic by default
grain_graph <- set.vertex.attribute(grain_graph, name='problematic', value=0)

# Read problem threshold config
grain_problems <- read.csv(cl_args$grainproblemconfig, header=TRUE, comment.char='#', na.strings="NA")

if (cl_args$timing) toc("Initializing")

#
# Diagnose graph for problems
#
if (cl_args$verbose) my_print("Diagnosing problems ...")
if (cl_args$timing) tic(type="elapsed")

diagnose_problem <- function(problem_type, problem_variable, problem_condition)
{
    if (problem_variable %in% colnames(task_data)) {
        problem_graph <- grain_graph
        problem_condition_full <- paste("!is.na(", problem_variable, ")", sep="")
        problem_condition_full <- paste(problem_condition_full, "&",problem_variable, problem_condition)
        problem_tasks <- subset(task_data, eval(parse(text=problem_condition_full)), select=task)

        problem_info_text <- paste(length(problem_tasks$task), "grains have", problem_condition_full)

        if ("on_crit_path" %in% colnames(task_data)) {
            problem_tasks_critical <- subset(task_data, eval(parse(text=paste(problem_condition_full, "& on_crit_path == 1"))), select=task)
            problem_info_text <- paste(problem_info_text, ", ", sep="")
            problem_info_text <- paste(problem_info_text, length(problem_tasks_critical$task), "of which are on the critical path.")
        }

        sink(graph_diagnosis_info_out_file, append=T)
        my_print(problem_info_text)
        sink()

        problem_tasks_index <- match(as.character(problem_tasks$task), V(problem_graph)$task)
        problem_graph <- set.vertex.attribute(problem_graph, name='color', index=problem_tasks_index, value="#FF0000")
        problem_graph <- set.vertex.attribute(problem_graph, name='problematic', index=problem_tasks_index, value=1)

        temp_out_file <- paste(gsub(". $", "", cl_args$out), paste("-",problem_type, "-problem.graphml", sep=""), sep="")
        res <- write.graph(problem_graph, file=temp_out_file, format="graphml")
        my_print(paste("Wrote file:", temp_out_file))

        # Layout in Sugiyama style and write to PDF
        if (cl_args$layout) {
            temp_out_file <- paste(gsub(". $", "", cl_args$out), paste("-",problem_type, "-problem.pdf", sep=""), sep="")
            sugiyama_layout <- layout_with_sugiyama(problem_graph, attributes="all")
            pdf(temp_out_file)
            res <- plot(problem_graph, layout=sugiyama_layout$layout)
            res <- dev.off()
            my_print(paste("Wrote file:", temp_out_file))
        }
    } else {
        my_print(paste("Warning: Problem", problem_type, "is not diagnosed since", problem_variable, "attribute is missing in graph!"))
    }

    return(NA)
}

ret_val <- mapply(diagnose_problem, grain_problems[,"type"], grain_problems[,"variable"], grain_problems[,"condition"])

if (cl_args$timing) toc("Diagnosing for problems")

#
# Write information about grain graph diagnosis
#
if (cl_args$verbose) my_print("Writing grain graph diagnosis information ...")
my_print(paste("Wrote file:", graph_diagnosis_info_out_file))

#
# Cleanup
#

# Warn
wa <- warnings()
if (class(wa) != "NULL")
    print(wa)
