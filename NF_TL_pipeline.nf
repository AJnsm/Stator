#!/usr/bin/env nextflow

cellTypes_ch = Channel.from( 7, 8, 10, 13, 18 ) 

process makeData {
    executor 'sge'
    cores 1
    memory '64 GB'
    time '10h'

    conda '/Users/s1855283/anaconda3/envs/NF_TL_env'
    publishDir "${projectDir}/plots", mode: 'copy'
        
    input:
    path dataScript from "${projectDir}/pipelineScripts/makeTrainingData.py"
    path rawData from "${projectDir}/cellData/1M_neurons_20k.h5"
    path clusters from "${projectDir}/cellData/clusters_20k.csv"
    path bcDoublets from "${projectDir}/cellData/bcDoublets_20k.csv"
    val cellType from cellTypes_ch
    
    output:
    path "trainingData_CL*.csv" into dataSets mode flatten
    path "*.png" into plots
    
    """
    python ${dataScript} ${rawData} ${clusters} ${params.nGenes} ${params.nCells} ${cellType} ${bcDoublets}
    """
}

process estimatePCgraph {
    executor 'sge'
    cores 16
    memory '2 GB'
    time '140h'

    conda '/Users/s1855283/anaconda3/envs/rEnv'
    publishDir "${projectDir}/output", mode: 'copy'

    input:
    path PCgraphEstScript from "${projectDir}/pipelineScripts/parallelPCscript.R"
    path dataSet from dataSets

    output:
    tuple path(dataSet), path('PCgraph*.csv') into PCgraphs_forMCMC_ch mode flatten
    tuple path(dataSet), path('*graph*.csv') into PC_and_ctrl_graphs_ch mode flatten

    """
    Rscript ${PCgraphEstScript} ${dataSet}
    """
 }

process iterMCMCscheme {
    executor 'sge'
    cores 1
    memory '8 GB'
    time '8h'

    conda '/Users/s1855283/anaconda3/envs/rEnv'
    publishDir "${projectDir}/output", mode: 'copy'
    
    input:
    path MCMCscript from "${projectDir}/pipelineScripts/iterMCMCscript.R"
    tuple path(dataSet), path(PCgraph) from PCgraphs_forMCMC_ch

    output:
    tuple path(dataSet), path('MCMCgraph*.csv') into MCMCgraphs_ch mode flatten

    """
    Rscript ${MCMCscript} ${PCgraph} ${dataSet} ${params.nGenes} 
    """
}

data_and_graphs_ch = PC_and_ctrl_graphs_ch.mix(MCMCgraphs_ch)


process estimateCoups {
    executor 'sge'
    cores 12
    memory '4 GB'
    time '10h'

    conda '/Users/s1855283/anaconda3/envs/NF_TL_env'
    publishDir "${projectDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/estimateTLcoups.py"
    tuple path(dataSet), path(graph) from data_and_graphs_ch
    
    output:
    path 'interactions*.npy' into interaction_ch
    

    """
    python ${estimationScript} ${dataSet} ${graph} 2 10 2
    """

}





