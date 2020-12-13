#first run:
#conda create -n NF_TL_env 
#conda activate NF_TL_env
# See also /Users/s1855283/anaconda3/envs/NF_TL_env/conda-meta/history|grep 'cmd'

conda install python=3.6
conda install numpy=1.19
conda install -c anaconda scipy
conda install -c conda-forge matplotlib
conda install -c anaconda pandas
conda install -c bioconda scanpy
conda install -c bioconda scrublet
conda install jupyter
conda install -c anaconda ipykernel
python -m ipykernel install --user --name=NF_TL_env
conda install -c anaconda mkl
conda install -c conda-forge python-igraph
