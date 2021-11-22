#!/usr/bin/env nextflow

cellTypes_ch = Channel.from(params.clusterArray) 

process makeData {

    publishDir "${launchDir}/plots", mode: 'copy', pattern: '*.png'
    publishDir "${launchDir}/embeddings", mode: 'copy', pattern: '*coords.csv'

    input:
    path dataScript from "${projectDir}/pipelineScripts/makeTrainingData.py"
    val cellType from cellTypes_ch
    path rawData from params.rawDataPath
    path clusters from params.clusterFile
    path userGenes from params.userGenes
    
    output:
    path "trainingData_*Genes.csv" into dataSets mode flatten
    path "*.png" optional true into plots
    path "*coords.csv" optional true into embeddings
    

    script:
    if( params.dataType == 'agnostic' )
        """
        python ${dataScript} --dataType ${params.dataType} --rawData ${rawData} --clusters ${clusters} --nGenes ${params.nGenes} --nCells ${params.nCells} --cluster ${cellType} --bcDoublets ${params.doubletFile} --userGenes ${userGenes} --twoReplicates ${params.twoReplicates}
        """

    else if( params.dataType == 'expression' )
        """
        python ${dataScript} --dataType ${params.dataType} --rawData ${rawData} --clusters ${clusters} --nGenes ${params.nGenes} --nCells ${params.nCells} --cluster ${cellType} --bcDoublets ${params.doubletFile} --userGenes ${userGenes} --twoReplicates ${params.twoReplicates}
        """

    else if( params.dataType == '10X' )
        """
        python ${dataScript} --dataType ${params.dataType} --rawData ${rawData} --clusters ${clusters} --nGenes ${params.nGenes} --nCells ${params.nCells} --cluster ${cellType} --bcDoublets ${params.doubletFile}
        """
    else if( params.dataType == 'Zeisel' )
        """
        python ${dataScript} --dataType ${params.dataType} --rawData ${rawData} --nCells ${params.nCells} --nGenes ${params.nGenes}
        """
    else
        error "Invalid data type"
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
    Rscript ${PCgraphEstScript} ${dataSet} ${params.cores_PC} ${params.PCalpha}
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
    python ${estimationScript} ${dataSet} ${graph} 1 ${params.bsResamps} ${params.cores_1pt} ${dipTest_pVals} ${params.estimationMethod} ${params.edgeListAlpha}
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
    path 'edgeList*.csv' into interaction_2pts_ch_edgeList

    """
    python ${estimationScript} ${dataSet} ${graph} 2 ${params.bsResamps} ${params.cores_2pt} ${dipTest_pVals} ${params.estimationMethod} ${params.edgeListAlpha}
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
    path 'edgeList*.csv' into interaction_3pts_ch_edgeList

    """
    python ${estimationScript} ${dataSet} ${graph} 3 ${params.bsResamps} ${params.cores_3pt} ${dipTest_pVals} ${params.estimationMethod} ${params.edgeListAlpha}
    """

}



