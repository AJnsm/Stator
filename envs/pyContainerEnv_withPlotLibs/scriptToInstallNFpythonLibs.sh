# could first create and activate env, but does not seem necessary, can also just install in base env.
# See also cat /Users/s1855283/anaconda3/envs/NF_TL_env/conda-meta/history|grep 'cmd'

conda config --add channels anaconda
conda config --add channels bioconda
conda config --add channels conda-forge

conda install python=3.6
conda install numpy=1.19

conda install -c anaconda scipy
conda install -c conda-forge matplotlib
conda install -c anaconda pandas
conda install -c bioconda scanpy
conda install -c conda-forge seaborn
conda install -c conda-forge upsetplot

pip install python-igraph
pip install hypernetx
pip install Pillow

