#!/usr/bin/env nextflow

cellTypes_ch = Channel.from(params.clusterArray) 

process makeData {

    publishDir "${launchDir}/plots", mode: 'copy', pattern: '*.png'
    publishDir "${launchDir}/embeddings", mode: 'copy', pattern: '*coords.csv'

    input:
    path dataScript from "${projectDir}/pipelineScripts/makeTrainingData.py" 
    val cellType from cellTypes_ch
    path rawData from params.rawDataPath
    path clusterFile from params.clusterFile
    path doubletFile from params.doubletFile
    path userGenes from params.userGenes
    
    output:
    path "trainingData_*Genes.csv" into dataSets mode flatten
    path "*.png" optional true into plots
    path "*PCAcoords.csv" optional true into PCAembeddings
    path "*UMAPcoords.csv" optional true into UMAPembeddings
    

    script:
    if( params.dataType == 'agnostic' )
        """
        python ${dataScript} --dataType ${params.dataType} --rawData ${rawData} --clusters ${clusterFile} --nGenes ${params.nGenes} --nCells ${params.nCells} --cluster ${cellType} --bcDoublets ${doubletFile} --userGenes ${userGenes} --twoReplicates ${params.twoReplicates}
        """

    else if( params.dataType == 'expression' )
        """
        python ${dataScript} --dataType ${params.dataType} --rawData ${rawData} --clusters ${clusterFile} --nGenes ${params.nGenes} --nCells ${params.nCells} --cluster ${cellType} --bcDoublets ${doubletFile} --userGenes ${userGenes} --twoReplicates ${params.twoReplicates}
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
    tuple path(dataSet), path('CTRLgraph*.csv') into CTRLgraphs_ch mode flatten

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
    path 'CPDAGgraph*.csv' into CPDAGgraphs_ch
    tuple path(dataSet), path('MCMCgraph*.csv') into MCMCgraphs_ch mode flatten

    """
    Rscript ${MCMCscript} ${PCgraph} ${dataSet} ${params.nGenes} 
    """
}

MCMCgraphs_ch.into {MCMCgraphs_ch1; MCMCgraphs_ch2}
data_and_graphs_ch = CTRLgraphs_ch.mix(MCMCgraphs_ch1)
data_and_graphs_ch.into {data_and_graphs_1pts; data_and_graphs_2pts; data_and_graphs_3pts; data_and_graphs_HOIs_MB; data_and_graphs_HOIs_6n7}


process estimateCoups_1pts {
    label 'interactionEstimation'    
    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/estimateTLcoups.py" 
    path utilities from "${projectDir}/pipelineScripts/utilities.py" 
    path genesToOne from params.genesToOne
    tuple path(dataSet), path(graph) from data_and_graphs_1pts
    
    output:
    path 'interactions*.npy' into interaction_1pts_ch
    

    """
    python ${estimationScript} --dataPath ${dataSet} --graphPath ${graph} --intOrder 1 --nResamps ${params.bsResamps} --nCores ${params.cores_1pt} --estimationMethod ${params.estimationMethod} --edgeListAlpha ${params.edgeListAlpha} --genesToOne ${genesToOne} --dataDups ${params.dataDups} --boundBool ${params.boundBool}
    """

}


process estimateCoups_2pts {
    label 'interactionEstimation'
    
    publishDir "${launchDir}/coupling_output", mode: 'copy', pattern: '*.npy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/estimateTLcoups.py" 
    path utilities from "${projectDir}/pipelineScripts/utilities.py" 
    path genesToOne from params.genesToOne
    tuple path(dataSet), path(graph) from data_and_graphs_2pts
    
    output:
    path 'interactions_order2_MCMCgraph*CI_F.npy' into interaction_2pts_CI_F_ch
    path 'interactions_order2_MCMCgraph*_undef.npy' into interaction_2pts_undef_ch
    path 'interactions_order2_MCMCgraph*_inf.npy' into interaction_2pts_inf_ch
    path 'interactions_order2_MCMCgraph*.npy' into interaction_2pts_ch
    path 'edgeList*.csv' into interaction_2pts_ch_edgeList

    """
    python ${estimationScript} --dataPath ${dataSet} --graphPath ${graph} --intOrder 2 --nResamps ${params.bsResamps} --nCores ${params.cores_2pt} --estimationMethod ${params.estimationMethod} --edgeListAlpha ${params.edgeListAlpha} --genesToOne ${genesToOne} --dataDups ${params.dataDups} --boundBool ${params.boundBool}
    """

}


process estimateCoups_3pts {
    label 'interactionEstimation'
    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/estimateTLcoups.py" 
    path utilities from "${projectDir}/pipelineScripts/utilities.py" 
    path genesToOne from params.genesToOne
    tuple path(dataSet), path(graph) from data_and_graphs_3pts
    
    output:
    path 'interactions*.npy' into interaction_3pts_ch
    path 'edgeList*.csv' into interaction_3pts_ch_edgeList

    """
    python ${estimationScript} --dataPath ${dataSet} --graphPath ${graph} --intOrder 3 --nResamps ${params.bsResamps} --nCores ${params.cores_3pt} --estimationMethod ${params.estimationMethod} --edgeListAlpha ${params.edgeListAlpha} --genesToOne ${genesToOne} --dataDups ${params.dataDups} --boundBool ${params.boundBool}
    """

}


process estimateCoups_345pts_WithinMB {
    label 'interactionEstimation'
    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/calcHOIsWithinMB.py" 
    path genesToOne from params.genesToOne
    path utilities from "${projectDir}/pipelineScripts/utilities.py" 
    tuple path(dataSet), path(graph) from data_and_graphs_HOIs_MB
    
    output:
    path 'interactions_withinMB_3pts*.npy' into interaction_withinMB_3pts
    path 'interactions_withinMB_4pts*.npy' into interaction_withinMB_4pts
    path 'interactions_withinMB_5pts*.npy' into interaction_withinMB_5pts

    """
    python ${estimationScript} --dataPath ${dataSet} --graphPath ${graph} --nResamps ${params.bsResamps} --nCores ${params.cores_HOIs_MB} --nRandoms ${params.nRandomHOIs} --genesToOne ${genesToOne} --dataDups ${params.dataDups} --boundBool ${params.boundBool}
    """

}

interaction_withinMB_5pts.into {interaction_withinMB_5pts_ch1; interaction_withinMB_5pts_ch2}

process estimateCoups_6n7pts {
    label 'interactionEstimation'
    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/calcHOIs_6n7pts.py" 
    path genesToOne from params.genesToOne
    path withinMB_5pts from interaction_withinMB_5pts_ch1
    path utilities from "${projectDir}/pipelineScripts/utilities.py" 
    tuple path(dataSet), path(graph) from data_and_graphs_HOIs_6n7
        
    output:
    path 'interactions*.npy' optional true into interaction_6n7pts

    """
    python ${estimationScript} --dataPath ${dataSet} --graphPath ${graph} --pathTo5pts ${withinMB_5pts} --nResamps ${params.bsResamps} --nCores ${params.cores_HOIs_6n7} --nRandoms ${params.nRandomHOIs} --genesToOne ${genesToOne} --dataDups ${params.dataDups} --boundBool ${params.boundBool}
    """

}


process createHOIsummaries {
    
    publishDir "${launchDir}/HOIplots", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/createHOIsummaries.py" 
    path utilities from "${projectDir}/pipelineScripts/utilities.py" 
    tuple path(dataSet), path(MCMCgraph) from MCMCgraphs_ch2
    path CPDAGgraph from CPDAGgraphs_ch
    path path2pts from interaction_2pts_ch
    path path2pts_CI_F from interaction_2pts_CI_F_ch
    path path2pts_undef from interaction_2pts_undef_ch
    path path2pts_inf from interaction_2pts_inf_ch
    path path3pts from interaction_withinMB_3pts
    path path4pts from interaction_withinMB_4pts
    path path5pts from interaction_withinMB_5pts_ch2
    path pcaCoords from PCAembeddings

    output:
    path '*.png' optional true into HOIsummaries

    """
    python ${estimationScript} --dataPath ${dataSet} --PCApath ${pcaCoords} --CPDAGgraphPath ${CPDAGgraph} --MCMCgraphPath ${MCMCgraph} --pathTo2pts ${path2pts} --pathTo2pts_CI_F ${path2pts_CI_F} --pathTo2pts_undef ${path2pts_undef} --pathTo2pts_inf ${path2pts_inf} --pathTo3pts ${path3pts} --pathTo4pts ${path4pts} --pathTo5pts ${path5pts}
    """

}





















