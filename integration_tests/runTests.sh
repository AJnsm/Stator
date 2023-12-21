#!/bin/bash

/Users/s1855283/Desktop/PhD/RBMs/MCMC/magneto-master/magneto.exe -L=3 -J=1 -TMin=2.26 -TSteps=1 -N1=10000 -N2=100000 -N3=100 -states=isingTestData -record=main


source /Users/s1855283/anaconda3/bin/activate NF_TL_env
python reshapeMagneto.py

python ../scripts/makeTrainingData.py \
    --dataType agnostic \
    --rawData isingTestData.csv \
    --nGenes 9 \
    --nCells 100000 \

source /Users/s1855283/anaconda3/bin/activate rEnv
Rscript ../scripts/parallelPCscript.R trainingData_CL01_100000Cells_0009Genes.csv 4 0.1


source /Users/s1855283/anaconda3/bin/deactivate
Rscript ../scripts/iterMCMCscript.R PCgraph_CL01_100000Cells_0009Genes.csv trainingData_CL01_100000Cells_0009Genes.csv 9


source /Users/s1855283/anaconda3/bin/activate NF_TL_env
python ../scripts/estimateTLcoups.py \
    --dataPath trainingData_CL01_100000Cells_0009Genes.csv \
    --graphPath MCMCgraph_CL01_100000Cells_0009Genes.csv \
    --intOrder 1 \
    --nResamps 1000 \
    --nCores 4 \
    --estimationMethod expectations \
    --dataDups 0 --boundBool 0

python ../scripts/estimateTLcoups.py \
    --dataPath trainingData_CL01_100000Cells_0009Genes.csv \
    --graphPath MCMCgraph_CL01_100000Cells_0009Genes.csv \
    --intOrder 2 \
    --nResamps 1000 \
    --nCores 4 \
    --estimationMethod expectations \
    --dataDups 0 --boundBool 0

python ../scripts/estimateTLcoups.py \
    --dataPath trainingData_CL01_100000Cells_0009Genes.csv \
    --graphPath isingTestData_trueAdjMat.csv \
    --intOrder 1 \
    --nResamps 1000 \
    --nCores 4 \
    --estimationMethod expectations \
    --dataDups 0 --boundBool 0

python ../scripts/estimateTLcoups.py \
    --dataPath trainingData_CL01_100000Cells_0009Genes.csv \
    --graphPath isingTestData_trueAdjMat.csv \
    --intOrder 2 \
    --nResamps 1000 \
    --nCores 4 \
    --estimationMethod expectations \
    --dataDups 0 --boundBool 0
