#!/usr/bin/env nextflow
nextflow.enable.dsl=1

process makeData {

    publishDir "${launchDir}/output", mode: 'copy'

    input:
    path dataScript from "${projectDir}/pipelineScripts/makeTrainingData.py" 
    path rawData from params.rawDataPath
    
    output:
    path "*.h5ad" into unbinarizedData
    path "trainingData_*Genes.csv" into dataSets mode flatten
    path "*.png" optional true into plots
    path "*PCAcoords.csv" optional true into PCAembeddings
    path "*UMAPcoords.csv" optional true into UMAPembeddings
    

    script:
    """
    python ${dataScript} \
    --dataType ${params.dataType} \
    --rawData ${rawData} \
    --nGenes ${params.nGenes} \
    --nCells ${params.nCells} \
    --bcDoublets ${params.doubletFile} \
    --userGenes ${params.userGenes} \
    --fracMito ${params.fracMito} \
    --minGenes ${params.minGenes} \
    --minCells ${params.minCells}
    """
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

MCMCgraphs_ch.into {MCMCgraphs_ch1; MCMCgraphs_ch2; MCMCgraphs_ch3; MCMCgraphs_ch4}
data_and_graphs_ch = CTRLgraphs_ch.mix(MCMCgraphs_ch1)
data_and_graphs_ch.into {data_and_graphs_1pts; data_and_graphs_2pts}


process estimateCoups_1pts {
    label 'interactionEstimation'    
    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/estimateTLcoups.py" 
    path utilities from "${projectDir}/pipelineScripts/utilities.py" 
    path genesToOne from params.genesToOne
    tuple path(dataSet), path(graph) from data_and_graphs_1pts
    
    output:
    path 'interactions*.npy'
    
    """
    python ${estimationScript} \
    --dataPath ${dataSet} \
    --graphPath ${graph} \
    --intOrder 1 \
    --nResamps ${params.bsResamps} \
    --nCores ${params.cores_1pt} \
    --estimationMethod ${params.estimationMethod} \
    --genesToOne ${genesToOne} \
    --dataDups ${params.dataDups} \
    --boundBool ${params.boundBool} \
    --asympBool ${params.asympBool}
    """

}


process estimateCoups_2pts {
    label 'interactionEstimation'
    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/estimateTLcoups.py" 
    path utilities from "${projectDir}/pipelineScripts/utilities.py" 
    path genesToOne from params.genesToOne
    tuple path(dataSet), path(graph) from data_and_graphs_2pts
    
    output:
    path 'interactions*.npy'

    script:
    """
    python ${estimationScript} \
    --dataPath ${dataSet} \
    --graphPath ${graph} \
    --intOrder 2 \
    --nResamps ${params.bsResamps} \
    --nCores ${params.cores_2pt} \
    --estimationMethod ${params.estimationMethod} \
    --genesToOne ${genesToOne} \
    --dataDups ${params.dataDups} \
    --boundBool ${params.boundBool} \
    --asympBool ${params.asympBool}
    """
}



process estimateCoups_2345pts_WithinMB {
    label 'interactionEstimation'
    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/calcHOIsWithinMB.py" 
    path genesToOne from params.genesToOne
    path utilities from "${projectDir}/pipelineScripts/utilities.py" 
    tuple path(dataSet), path(graph) from MCMCgraphs_ch2
    
    output:
    path 'interactions_withinMB_2pts*.npy' into interaction_withinMB_2pts
    path 'interactions_withinMB_3pts*.npy' into interaction_withinMB_3pts
    path 'interactions_withinMB_4pts*.npy' into interaction_withinMB_4pts
    path 'interactions_withinMB_5pts*.npy' into interaction_withinMB_5pts

    """
    python ${estimationScript} \
    --dataPath ${dataSet} \
    --graphPath ${graph} \
    --nResamps ${params.bsResamps} \
    --nCores ${params.cores_HOIs_MB} \
    --nRandoms ${params.nRandomHOIs} \
    --genesToOne ${genesToOne} \
    --dataDups ${params.dataDups} \
    --boundBool ${params.boundBool} \
    --asympBool ${params.asympBool}
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
    tuple path(dataSet), path(graph) from MCMCgraphs_ch3
        
    output:
    path 'interactions*.npy' optional true into interaction_6n7pts

    """
    python ${estimationScript} \
    --dataPath ${dataSet} \
    --graphPath ${graph} \
    --pathTo5pts ${withinMB_5pts} \
    --nResamps ${params.bsResamps} \
    --nCores ${params.cores_HOIs_6n7} \
    --nRandoms ${params.nRandomHOIs} \
    --genesToOne ${genesToOne} \
    --dataDups ${params.dataDups} \
    --boundBool ${params.boundBool} \
    --asympBool ${params.asympBool}
    """

}


process createHOIsummaries {
    
    publishDir "${launchDir}/HOIsummaries", mode: 'copy', pattern: '*.csv'
    publishDir "${launchDir}/HOIsummaries", mode: 'copy', pattern: '*_summary.png'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/createHOIsummaries.py" 
    path utilities from "${projectDir}/pipelineScripts/utilities.py" 
    tuple path(dataSet), path(MCMCgraph) from MCMCgraphs_ch4
    path CPDAGgraph from CPDAGgraphs_ch
    // path path2pts from interaction_2pts_ch
    // path path2pts_CI_F from interaction_2pts_CI_F_ch
    // path path2pts_undef from interaction_2pts_undef_ch
    // path path2pts_inf from interaction_2pts_inf_ch
    path path2pts from interaction_withinMB_2pts
    path path3pts from interaction_withinMB_3pts
    path path4pts from interaction_withinMB_4pts
    path path5pts from interaction_withinMB_5pts_ch2
    path pcaCoords from PCAembeddings

    output:
    path '*.png' optional true into HOIsummaries
    path 'top_DTuples.csv' optional true into topDeviators
    path 'all_DTuples.csv' optional true into allDeviators_csv
    path 'DTuples_binaryReps.csv' optional true into binaryReps_csv
    path dataSet into dataSet_forPlots
    path pcaCoords into PCAembeddings_forPlots

    """
    python ${estimationScript} \
    --dataPath ${dataSet} \
    --PCApath ${pcaCoords} \
    --CPDAGgraphPath ${CPDAGgraph} \
    --MCMCgraphPath ${MCMCgraph} \
    --pathTo2pts ${path2pts} \
    --pathTo3pts ${path3pts} \
    --pathTo4pts ${path4pts} \
    --pathTo5pts ${path5pts} \
    --sigHOIthreshold ${params.sigHOIthreshold} \
    --minStateDeviation ${params.minStateDeviation} \
    --stateDevAlpha ${params.stateDevAlpha}  \
    --plotPairwiseUpsets ${params.plotPairwiseUpsets}
    """

}


process identifyStates {
    
    publishDir "${launchDir}/states_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/pipelineScripts/identifyStates.py" 
    path devStates from topDeviators
    path dataSet from dataSet_forPlots
    path pcaCoords from PCAembeddings_forPlots

    output:
    path '*.png' optional true into stateID_dendograms
    path '*.csv' optional true into stateID_outputs
    """
    python ${estimationScript} \
    --dataPath ${dataSet} \
    --PCApath ${pcaCoords} \
    --devStates ${devStates} \
    --diffCutoff ${params.dendCutoff} \
    --bsResamps ${params.bsResamps_HC} \
    --auThreshold ${params.auThreshold}
    """

} 

















