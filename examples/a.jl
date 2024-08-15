#using DataDeps
#using CodecZlib, Tar
#import PublicationDatasets as DS

using LinearAlgebra
using Serialization
#using Test

#import JSON3
#import BSON
import Random

import PythonPlot as PLT # do Pkg.add("PythonPlot") if this is missing.

using Revise

import NMRHamiltonian as HAM
#const HAM = NMRHamiltonian
#const SL = HAM.SL
#const Graphs = HAM.Graphs