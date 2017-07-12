#!/bin/bash

# Outputs
working_dir=new
out_dir_worker_1=worker-1
out_dir_worker_2=worker-2
out_dir_graph=graph

mkdir -p $working_dir

cd $working_dir

# Program parameters
PROG=$MIR_ROOT/examples/OMP/for-loop/for-loop-opt.out
INPUT="8 42 4"

# MIR parameters
STACK=128
SCHED=ws
MEMPOL=system

echo "Profiling on 1 worker ..."
if [ -d $out_dir_worker_1 ]; then
    echo "Directory $out_dir_worker_1 exists. Making a backup copy ..."
    if [ -d $out_dir_worker_1.backup ]; then
        echo "Directory $out_dir_worker_1.backup exists. Aborting!"
        exit
    fi
    mv $out_dir_worker_1 $out_dir_worker_1.backup
    echo "Moved $out_dir_worker_1 to $out_dir_worker_1.backup"
fi

mkdir -p $out_dir_worker_1

touch timestamp

MIR_CONF="--workers=1 --schedule=$SCHED --stack-size=$STACK --memory-policy=$MEMPOL --chunks-are-tasks --idle-task --worker-stats --task-stats -r" ${PROG} ${INPUT}

$MIR_ROOT/scripts/profiling/thread/get-events-per-task.sh mir-recorder-trace-*.rec

Rscript $MIR_ROOT/scripts/profiling/task/merge-task-stats.R -l mir-task-stats -r events-per-task-summary.csv -k task -c left

Rscript $MIR_ROOT/scripts/profiling/task/process-task-stats.R -d task-stats.merged --lineage --forloop

find . -maxdepth 1 -type f -newer timestamp -exec mv {} $out_dir_worker_1 \;

echo "Profiling on 2 workers ..."
if [ -d $out_dir_worker_2 ]; then
    echo "Directory $out_dir_worker_2 exists. Making a backup copy ..."
    if [ -d $out_dir_worker_2.backup ]; then
        echo "Directory $out_dir_worker_2.backup exists. Aborting!"
        exit
    fi
    mv $out_dir_worker_2 $out_dir_worker_2.backup
    echo "Moved $out_dir_worker_2 to $out_dir_worker_2.backup"
fi

mkdir -p $out_dir_worker_2

touch timestamp

MIR_CONF="--workers=2 --schedule=$SCHED --stack-size=$STACK --memory-policy=$MEMPOL --chunks-are-tasks --idle-task --worker-stats --task-stats -r" ${PROG} ${INPUT}

$MIR_ROOT/scripts/profiling/thread/get-events-per-task.sh mir-recorder-trace-*.rec

Rscript $MIR_ROOT/scripts/profiling/task/merge-task-stats.R -l mir-task-stats -r events-per-task-summary.csv -k task -c left

Rscript $MIR_ROOT/scripts/profiling/task/process-task-stats.R -d task-stats.merged --lineage --forloop --linenumbers -e ${PROG}

Rscript $MIR_ROOT/scripts/profiling/task/compare-task-stats.R -l task-stats.processed -r ${out_dir_worker_1}/task-stats.processed -k chunk_lineage

Rscript $MIR_ROOT/scripts/profiling/task/merge-task-stats.R -l task-stats.processed -r task-stats.compared -k chunk_lineage

find . -maxdepth 1 -type f -newer timestamp -exec mv {} $out_dir_worker_2 \;

echo "Plotting grain graph ..."
if [ -d $out_dir_graph ]; then
    echo "Directory $out_dir_graph exists. Making a backup copy ..."
    if [ -d $out_dir_graph.backup ]; then
        echo "Directory $out_dir_graph.backup exists. Aborting!"
        exit
    fi
    mv $out_dir_graph $out_dir_graph.backup
    echo "Moved $out_dir_graph to $out_dir_graph.backup"
fi

mkdir -p $out_dir_graph

touch timestamp

echo "Plotting graph ..."
Rscript $GRAIN_GRAPHS_ROOT/prototype/plot-grain-graph.R -d $out_dir_worker_2/task-stats.merged --grainpropertyconfig=$GRAIN_GRAPHS_ROOT/prototype/configs/grain-properties-for-loop.cfg --edgepropertyconfig=$GRAIN_GRAPHS_ROOT/prototype/configs/edge-properties.cfg --enumcriticalpath

echo "Diagnosing problems ..."
Rscript $GRAIN_GRAPHS_ROOT/prototype/diagnose-performance-problems.R --graph=grain-graph.graphml --grainproblemconfig=$GRAIN_GRAPHS_ROOT/prototype/configs/grain-problems.cfg

echo "Converting colormaps to PDF ..."
find . -maxdepth 1 -name "*.colormap" | xargs -n1 basename | xargs -n1 Rscript $GRAIN_GRAPHS_ROOT/prototype/plot-colormap.R -d

find . -maxdepth 1 -type f -newer timestamp -exec mv {} $out_dir_graph \;

rm timestamp
