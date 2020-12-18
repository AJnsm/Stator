#!/usr/bin/env nextflow

cellTypes_ch = Channel.from( 7, 8, 10, 13, 18 ) 

process makeData {

    
    publishDir "${launchDir}/plots", mode: 'copy', pattern: '*.png'
        
    input:
    path dataScript from "${projectDir}/pipelineScripts/makeTrainingData.py"
    path rawData from params.h5File
    path clusters from params.clusterFile
    path bcDoublets from params.doubletFile
    val cellType from cellTypes_ch
    
    output:
    path "trainingData_CL*.csv" into dataSets mode flatten
    path "*.png" into plots
    
    """
    python ${dataScript} ${rawData} ${clusters} ${params.nGenes} ${params.nCells} ${cellType} ${bcDoublets}
    """
}

process estimatePCgraph {

    
    publishDir "${launchDir}/output", mode: 'copy'

    input:
    path PCgraphEstScript from "${projectDir}/pipelineScripts/parallelPCscript.R"
    path dataSet from dataSets

    output:
    tuple path(dataSet), path('PCgraph*.csv') into PCgraphs_forMCMC_ch mode flatten
    tuple path(dataSet), path('*graph*.csv') into PC_and_ctrl_graphs_ch mode flatten

    """
    Rscript ${PCgraphEstScript} ${dataSet} ${params.PCcores}
    """
 }

process iterMCMCscheme {

    
    publishDir "${launchDir}/output", mode: 'copy'
    
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

    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/estimateTLcoups.py"
    tuple path(dataSet), path(graph) from data_and_graphs_ch
    
    output:
    path 'interactions*.npy' into interaction_ch
    

    """
    python ${estimationScript} ${dataSet} ${graph} 2 10 2
    """

}





