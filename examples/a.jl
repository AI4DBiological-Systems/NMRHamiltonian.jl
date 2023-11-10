using DataDeps
using CodecZlib, Tar

using LinearAlgebra

using Test

import JSON3
import BSON
import Random


using Revise

import NMRHamiltonian
const HAM = NMRHamiltonian
const SL = HAM.SL