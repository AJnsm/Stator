#!/usr/bin/env nextflow
nextflow.enable.dsl=1

formattedNCells = String.format( "%05d", params.nCells )
formattedNGenes = String.format( "%04d", params.nGenes )
dataSetID = "${formattedNCells}Cells_${formattedNGenes}Genes"


process makeData {

    publishDir "${launchDir}/output", mode: 'copy'

    input:
    path dataScript from "${projectDir}/scripts/makeTrainingData.py" 
    path rawData from params.rawDataPath
    
    output:
    path "unbinarised_cell_data.h5ad"
    path "trainingData_${dataSetID}.csv" into dataSets mode flatten
    path "*.png" optional true
    path "*PCAcoords.csv" into PCAembeddings
    path "*UMAPcoords.csv" into UMAPembeddings
    

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
    path PCgraphEstScript from "${projectDir}/scripts/parallelPCscript.R" 
    path dataSet from dataSets

    output:
    tuple path(dataSet), path("PCgraph_${dataSetID}.csv") into PCgraphs_forMCMC_ch mode flatten
    tuple path(dataSet), path("CTRLgraph_${dataSetID}.csv") into CTRLgraphs_ch mode flatten

    """
    Rscript ${PCgraphEstScript} ${dataSet} ${params.cores_PC} ${params.PCalpha}
    """
 }

process iterMCMCscheme {

    
    publishDir "${launchDir}/output", mode: 'copy'
    
    input:
    path MCMCscript from "${projectDir}/scripts/iterMCMCscript.R" 
    tuple path(dataSet), path(PCgraph) from PCgraphs_forMCMC_ch

    output:
    path "CPDAGgraph_${dataSetID}.csv" into CPDAGgraphs_ch
    tuple path(dataSet), path("MCMCgraph_${dataSetID}.csv") into MCMCgraphs_ch mode flatten

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
    path estimationScript from "${projectDir}/scripts/estimateTLcoups.py" 
    path utilities from "${projectDir}/scripts/utilities.py" 
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
    path estimationScript from "${projectDir}/scripts/estimateTLcoups.py" 
    path utilities from "${projectDir}/scripts/utilities.py" 
    path genesToOne from params.genesToOne
    tuple path(dataSet), path(graph) from data_and_graphs_2pts
    
    output:
    path 'interactions*.npy' optional true

    script:
    if( params.calcAll2pts == 1 )
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
    else
        """
        echo skipping calculation of all 2-point interactions
        """
}



process estimateCoups_2345pts_WithinMB {
    label 'interactionEstimation'
    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/scripts/calcHOIsWithinMB.py" 
    path genesToOne from params.genesToOne
    path utilities from "${projectDir}/scripts/utilities.py" 
    tuple path(dataSet), path(graph) from MCMCgraphs_ch2
    
    output:
    path "interactions_withinMB_2pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy" into interaction_withinMB_2pts
    path "interactions_withinMB_3pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy" into interaction_withinMB_3pts
    path "interactions_withinMB_4pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy" into interaction_withinMB_4pts
    path "interactions_withinMB_5pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy" into interaction_withinMB_5pts

    path "interactions_random_2pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy" into interaction_random_2pts
    path "interactions_random_3pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy" into interaction_random_3pts
    path "interactions_random_4pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy" into interaction_random_4pts
    path "interactions_random_5pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy" into interaction_random_5pts

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
    --asympBool ${params.asympBool} \
    --estimationMode ${params.estimationMode}
    """

}

interaction_withinMB_5pts.into {interaction_withinMB_5pts_ch1; interaction_withinMB_5pts_ch2}

process estimateCoups_6n7pts {
    label 'interactionEstimation'
    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript from "${projectDir}/scripts/calcHOIs_6n7pts.py" 
    path genesToOne from params.genesToOne
    path withinMB_5pts from interaction_withinMB_5pts_ch1
    path utilities from "${projectDir}/scripts/utilities.py" 
    tuple path(dataSet), path(graph) from MCMCgraphs_ch3
        
    output:
    path 'interactions*.npy' optional true

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
    path estimationScript from "${projectDir}/scripts/createHOIsummaries.py" 
    path utilities from "${projectDir}/scripts/utilities.py" 
    tuple path(dataSet), path(MCMCgraph) from MCMCgraphs_ch4
    path CPDAGgraph from CPDAGgraphs_ch

    path path2pts from interaction_withinMB_2pts
    path path3pts from interaction_withinMB_3pts
    path path4pts from interaction_withinMB_4pts
    path path5pts from interaction_withinMB_5pts_ch2
    path pcaCoords from PCAembeddings

    output:
    path '*.png' optional true
    path 'top_DTuples.csv' into topDeviators
    path 'all_DTuples.csv' into allDeviators_csv
    path 'DTuples_binaryReps.csv' optional true
    path dataSet into dataSet_forPlots
    path pcaCoords into PCAembeddings_forPlots

    """
    python ${estimationScript} \
    --dataPath ${dataSet} \
    --PCApath ${pcaCoords} \
    --CPDAGgraphPath ${CPDAGgraph} \
    --MCMCgraphPath ${MCMCgraph} \
    --estimationMode ${params.estimationMode} \
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
    path estimationScript from "${projectDir}/scripts/identifyStates.py" 
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

















