#!/usr/bin/env nextflow

cellTypes_ch = Channel.from(7) 

process makeData {

    
    publishDir "${launchDir}/plots", mode: 'copy', pattern: '*.png'
        
    input:
    path dataScript from "${projectDir}/pipelineScripts/makeTrainingData.py"
    path rawData from params.h5File
    path clusters from params.clusterFile
    path bcDoublets from params.doubletFile
    val cellType from cellTypes_ch
    
    output:
    path "trainingData_CL*Genes.csv" into dataSets mode flatten
    path "*.png" into plots
    path "*coords.csv" into embeddings
    
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
    Rscript ${PCgraphEstScript} ${dataSet} ${params.PCcores} ${params.PCalpha}
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
data_and_graphs_ch.into {data_and_graphs_1pts; data_and_graphs_2pts; data_and_graphs_3pts}


process estimateCoups_1pts {
    label 'interactionEstimation'    
    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/estimateTLcoups.py"
    path dipTest_pVals from "${projectDir}/dipPvals.csv"
    tuple path(dataSet), path(graph) from data_and_graphs_1pts
    
    output:
    path 'interactions*.npy' into interaction_1pts_ch
    

    """
    python ${estimationScript} ${dataSet} ${graph} 1 ${params.bsResamps} ${params.coupCores} ${dipTest_pVals} ${params.estimationMethod}
    """

}


process estimateCoups_2pts {
    label 'interactionEstimation'
    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/estimateTLcoups.py"
    path dipTest_pVals from "${projectDir}/dipPvals.csv"
    tuple path(dataSet), path(graph) from data_and_graphs_2pts
    
    output:
    path 'interactions*.npy' into interaction_2pts_ch
    

    """
    python ${estimationScript} ${dataSet} ${graph} 2 ${params.bsResamps} ${params.coupCores} ${dipTest_pVals} ${params.estimationMethod}
    """

}


process estimateCoups_3pts {
    label 'interactionEstimation'
    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/estimateTLcoups.py"
    path dipTest_pVals from "${projectDir}/dipPvals.csv"
    tuple path(dataSet), path(graph) from data_and_graphs_3pts
    
    output:
    path 'interactions*.npy' into interaction_3pts_ch
    

    """
    python ${estimationScript} ${dataSet} ${graph} 3 ${params.bsResamps} ${params.coupCores} ${dipTest_pVals} ${params.estimationMethod}
    """

}



