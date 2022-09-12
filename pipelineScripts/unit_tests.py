'''
A collection of unit tests for the interaction estimation utility functions. 
Running `pytest unit_tests.py` from the CL generates a report. 
Compiling the numba functions can take a few minutes. 
'''

from utilities import *
import igraph as ig
import numpy as np
import pandas as pd
import pytest

# ****** test graph with chain, fork, and collider ******
# Graph structure: 0 -> 1 <- 2 -> 3 -> 4
@pytest.fixture
def graphData():
    g = ig.Graph(directed=True)
    g.add_vertices(5)
    for i in range(len(g.vs())):
        g.vs()[i]['label']=i
    g.add_edges([(0, 1), (2, 1), (2, 3), (3, 4)])
    pytest.graph = g

# ****** test Markov blanket identification ******
def test_markovBlanket_collider(graphData):
    assert findMarkovBlanket(0, pytest.graph)==[1, 2]

def test_markovBlanket_collider_chain(graphData):
    assert findMarkovBlanket(2, pytest.graph)==[0, 1, 3]

def test_markovBlanket_chain(graphData):
    assert findMarkovBlanket(3, pytest.graph)==[2, 4]

# ****** test Markov blanket conditioning ******
def test_conditionOnMB(graphData):
    vec = np.array([0, 0, 0, 1, 1, 1])
    dataSet = pd.DataFrame(np.array([vec, vec, vec[::-1], vec, vec]).T)
    assert all(conditionOnMB([2], pytest.graph, dataSet)[2].values==np.array([1, 1, 1]))

def test_conditionOnMB_conditionOnOne(graphData):
    vec = np.array([0, 0, 0, 1, 1, 1])
    dataSet = pd.DataFrame(np.array([vec, vec, vec[::-1], vec, vec]).T)
    assert all(conditionOnMB([2], pytest.graph, dataSet, genesToOne=[0, 1, 3])[2].values==np.array([0, 0, 0]))



# Testing the different estimation functions. 
# It is asserted that all of these have correct and identical behaviour

# ****** test baseline estimation ******
def test_calcInteraction_expectations_pInf():
    vals = pd.DataFrame([1, 1, 1])
    assert calcInteraction_expectations(vals)==np.inf

def test_calcInteraction_expectations_mInf():
    vals = pd.DataFrame([0, 0, 0])
    assert calcInteraction_expectations(vals)==-np.inf

def test_calcInteraction_expectations_zero():
    vals = pd.DataFrame([1, 1, 0, 0])
    assert calcInteraction_expectations(vals)==0.

def test_calcInteraction_expectations_nan():
    vals = pd.DataFrame([])
    assert np.isnan(calcInteraction_expectations(vals))

def test_calcInteraction_expectations_I_1():
    vals = pd.DataFrame(np.hstack([np.ones(10), np.zeros(10)]))
    assert (1 - np.exp(calcInteraction_expectations(vals)))<1e-5

def test_calcInteraction_expectations_I_2():
    vals = pd.DataFrame(np.hstack([np.ones(10), np.zeros(5)]))
    assert (2 - np.exp(calcInteraction_expectations(vals)))<1e-5



# ****** test numpy estimation ******
def test_calcInteraction_expectations_np_pInf():
    vals = pd.DataFrame([1, 1, 1])
    assert calcInteraction_expectations_np(vals)==np.inf

def test_calcInteraction_expectations_np_mInf():
    vals = pd.DataFrame([0, 0, 0])
    assert calcInteraction_expectations_np(vals)==-np.inf

def test_calcInteraction_expectations_np_zero():
    vals = pd.DataFrame([1, 1, 0, 0])
    assert calcInteraction_expectations_np(vals)==0.

def test_calcInteraction_expectations_np_nan():
    vals = pd.DataFrame([])
    assert np.isnan(calcInteraction_expectations_np(vals))

def test_calcInteraction_expectations_np_I_1():
    vals = pd.DataFrame(np.hstack([np.ones(10), np.zeros(10)]))
    assert (1 - np.exp(calcInteraction_expectations_np(vals)))<1e-5

def test_calcInteraction_expectations_np_I_2():
    vals = pd.DataFrame(np.hstack([np.ones(10), np.zeros(5)]))
    assert (2 - np.exp(calcInteraction_expectations_np(vals)))<1e-5



# ****** test numba estimation ******
# (compiling the numba function can take a few minutes)
def test_calcInteraction_expectations_numba_pInf():
    vals = pd.DataFrame(np.ones((4))).values
    assert calcInteraction_expectations_numba(vals)==np.inf

def test_calcInteraction_expectations_numba_mInf():
    vals = pd.DataFrame(np.zeros((4))).values
    assert calcInteraction_expectations_numba(vals)==-np.inf

def test_calcInteraction_expectations_numba_zero():
    vals = pd.DataFrame(np.array([1, 1, 0, 0])).values
    assert calcInteraction_expectations_numba(vals)==0.

def test_calcInteraction_expectations_numba_nan():
    vals = pd.DataFrame(np.array([])).values
    assert np.isnan(calcInteraction_expectations_numba(vals))

def test_calcInteraction_expectations_numba_I_1():
    vals = pd.DataFrame(np.hstack([np.ones(10), np.zeros(10)])).values
    assert (1 - np.exp(calcInteraction_expectations_numba(vals)))<1e-5

def test_calcInteraction_expectations_numba_I_2():
    vals = pd.DataFrame(np.hstack([np.ones(10), np.zeros(5)])).values
    assert (2 - np.exp(calcInteraction_expectations_numba(vals)))<1e-5