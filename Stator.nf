#!/usr/bin/env nextflow
nextflow.enable.dsl=2

formattedNCells = String.format( "%05d", params.nCells )
formattedNGenes = String.format( "%04d", params.nGenes )
dataSetID = "${formattedNCells}Cells_*Genes"


process makeData {
    label 'python'
    publishDir "${launchDir}/output", mode: 'copy'

    input:
    path dataScript 
    path rawData
    path userGenes
    path doubletFile
    
    output:
    path "unbinarised_cell_data.h5ad"
    path "trainingData_${dataSetID}.csv", emit: trainingData
    path "*.png" optional true
    path "*PCAcoords.csv", emit: PCAembeddings
    path "*UMAPcoords.csv", emit: UMAPembeddings
    

    script:
    """
    python ${dataScript} \
    --dataType ${params.dataType} \
    --rawData ${rawData} \
    --nGenes ${params.nGenes} \
    --nCells ${params.nCells} \
    --bcDoublets ${doubletFile} \
    --userGenes ${userGenes} \
    --fracMito ${params.fracMito} \
    --minGenes ${params.minGenes} \
    --minCells ${params.minCells}
    """
}

process estimatePCgraph {
    label 'R'
    publishDir "${launchDir}/output", mode: 'copy'

    input:
    path PCgraphEstScript
    path dataSet

    output:
    path("PCgraph_${dataSetID}.csv"), emit: PCgraph
    
    """
    Rscript ${PCgraphEstScript} ${dataSet} ${params.cores_PC ?: 1} ${params.PCalpha}
    """
 }


process iterMCMCscheme {
    label 'R'
    
    publishDir "${launchDir}/output", mode: 'copy'
    
    input:
    path MCMCscript
    path PCgraph
    path dataSet

    output:
    path "MCMCgraph_${dataSetID}.csv", emit: MCMCgraph
    path "CPDAGgraph_${dataSetID}.csv", emit: CPDAGgraph

    """
    Rscript ${MCMCscript} ${PCgraph} ${dataSet} ${params.nGenes} 
    """
}

process estimateCoups_2345pts_WithinMB {
    label 'python'
    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript
    path genesToOne
    path utilities
    path MCMCgraph
    path dataSet
    
    output:
    path "interactions_withinMB_2pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy", emit: interactions_withinMB_2pts
    path "interactions_withinMB_3pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy", emit: interactions_withinMB_3pts
    path "interactions_withinMB_4pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy", emit: interactions_withinMB_4pts
    path "interactions_withinMB_5pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy", emit: interactions_withinMB_5pts

    path "interactions_random_2pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy", emit: interactions_random_2pts
    path "interactions_random_3pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy", emit: interactions_random_3pts
    path "interactions_random_4pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy", emit: interactions_random_4pts
    path "interactions_random_5pts_${params.estimationMode}_MCMCgraph_${dataSetID}.npy", emit: interactions_random_5pts

    """
    python ${estimationScript} \
    --dataPath ${dataSet} \
    --graphPath ${MCMCgraph} \
    --nResamps ${params.bsResamps} \
    --nCores ${params.cores_HOIs_MB ?: 1} \
    --nRandoms ${params.nRandomHOIs} \
    --genesToOne ${genesToOne} \
    --dataDups ${params.dataDups} \
    --boundBool ${params.boundBool} \
    --asympBool ${params.asympBool} \
    --estimationMode ${params.estimationMode}
    """

}

process estimateCoups_6n7pts {
    label 'python'
    
    publishDir "${launchDir}/coupling_output", mode: 'copy'

    input:
    path estimationScript
    path genesToOne
    path withinMB_5pts
    path utilities
    path dataSet
    path MCMCgraph
        
    output:
    path 'interactions*.npy' optional true

    """
    python ${estimationScript} \
    --dataPath ${dataSet} \
    --graphPath ${MCMCgraph} \
    --pathTo5pts ${withinMB_5pts} \
    --nResamps ${params.bsResamps} \
    --nCores ${params.cores_HOIs_6n7 ?: 1} \
    --nRandoms ${params.nRandomHOIs} \
    --genesToOne ${genesToOne} \
    --dataDups ${params.dataDups} \
    --boundBool ${params.boundBool} \
    --asympBool ${params.asympBool}
    """

}

process identifyDTuples {
    label 'python'
    publishDir "${launchDir}/dtuples_output", mode: 'copy', pattern: '*.csv'
    publishDir "${launchDir}/dtuples_output", mode: 'copy', pattern: '*_summary.png'

    input:
    path estimationScript
    path utilities
    path dataSet
    path MCMCgraph
    path CPDAGgraph

    path path2pts 
    path path3pts 
    path path4pts 
    path path5pts 
    path pcaCoords

    output:
    path '*.png' optional true
    // path 'top_DTuples.csv', emit: topDeviators
    path 'all_DTuples.csv', emit: allDeviators
    path 'DTuples_binaryReps.csv' optional true

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

// This process is deprecated and replaced by the Shiny R app
// process identifyStates {
//     label 'python'
//     publishDir "${launchDir}/states_output", mode: 'copy'

//     input:
//     path estimationScript
//     path devStates
//     path pcaCoords
//     path dataSet

//     output:
//     path '*.png' optional true
//     path '*.csv' optional true
//     """
//     python ${estimationScript} \
//     --dataPath ${dataSet} \
//     --PCApath ${pcaCoords} \
//     --devStates ${devStates} \
//     --diffCutoff ${params.dendCutoff} \
//     --bsResamps ${params.bsResamps_HC} \
//     --auThreshold ${params.auThreshold}
//     """

// } 


// NOTE: Nextflow does not yet allow optional inputs, so we create some temporary empty files to pass to the process when there are no userGenes, doubletFile, or genesToOne files specified. 
workflow {
    script_makeTrainingData = "${projectDir}/scripts/makeTrainingData.py" 
    script_calcHOIsWithinMB = "${projectDir}/scripts/calcHOIsWithinMB.py"
    script_identifyDTuples = "${projectDir}/scripts/identifyDTuples.py" 
    script_identifyStates = "${projectDir}/scripts/identifyStates.py"
    script_iterMCMC = "${projectDir}/scripts/iterMCMC.R"
    script_parallelPC = "${projectDir}/scripts/parallelPC.R"
    script_identifyStates = "${projectDir}/scripts/identifyStates.py"
    script_calcHOIs_6n7pts = "${projectDir}/scripts/calcHOIs_6n7pts.py"
    utils = "${projectDir}/scripts/utilities.py"

    ch_empty_userGenes = file("${projectDir}/EMPTY_userGenes")
    ch_empty_doubletFile = file("${projectDir}/EMPTY_doubletFile")
    ch_empty_genesToOne = file("${projectDir}/EMPTY_genesToOne")

    makeData(script_makeTrainingData,
            params.rawDataPath,
            params.userGenes ? file(params.userGenes) : ch_empty_userGenes,
            params.doubletFile ? file(params.doubletFile) : ch_empty_doubletFile)

    estimatePCgraph(script_parallelPC, makeData.out.trainingData)

    iterMCMCscheme(script_iterMCMC, estimatePCgraph.out.PCgraph, 
                        makeData.out.trainingData)
    estimateCoups_2345pts_WithinMB(script_calcHOIsWithinMB, 
                        params.genesToOne ? file(params.genesToOne) : ch_empty_genesToOne,
                        utils,
                        iterMCMCscheme.out.MCMCgraph,
                        makeData.out.trainingData)

    estimateCoups_6n7pts(script_calcHOIs_6n7pts,
                        params.genesToOne ? file(params.genesToOne) : ch_empty_genesToOne,
                        estimateCoups_2345pts_WithinMB.out.interactions_withinMB_5pts,
                        utils,
                        makeData.out.trainingData,
                        iterMCMCscheme.out.MCMCgraph)

    identifyDTuples(script_identifyDTuples,
                        utils,
                        makeData.out.trainingData,
                        iterMCMCscheme.out.MCMCgraph,
                        iterMCMCscheme.out.CPDAGgraph,
                        estimateCoups_2345pts_WithinMB.out.interactions_withinMB_2pts,
                        estimateCoups_2345pts_WithinMB.out.interactions_withinMB_3pts,
                        estimateCoups_2345pts_WithinMB.out.interactions_withinMB_4pts,
                        estimateCoups_2345pts_WithinMB.out.interactions_withinMB_5pts,
                        makeData.out.PCAembeddings)

    // identifyStates(script_identifyStates,  
    //                     identifyDTuples.out.topDeviators,
    //                     makeData.out.PCAembeddings,
    //                     makeData.out.trainingData)
}

