abstract AbstractSDPModel <: MathProgBase.AbstractMathProgModel
export AbstractSDPModel

abstract SDPSolver <: MathProgBase.AbstractMathProgSolver

verbose(solver::SDPSolver) = solver.verbose


abstract SDPAGEN <: SDPSolver

executable(solver::SDPAGEN) = solver.executable


immutable SDPA <: SDPAGEN
    speed::Int
    executable::String
    verbose::Bool
end

SDPA(; speed::Integer=0, 
       executable::String="sdpa", 
       verbose::Bool=false) = SDPA(speed, executable, verbose)


immutable SDPAQD <: SDPAGEN
    speed::Int
    executable::String
    verbose::Bool
end

SDPAQD(; speed::Integer=0, 
         executable::String="sdpa_qd", 
         verbose::Bool=false) = SDPAQD(speed, executable, verbose)


immutable SDPAGMP <: SDPAGEN
    speed::Int
    executable::String
    verbose::Bool
end

SDPAGMP(; speed::Integer=0, 
          executable::String="sdpa_gmp", 
          verbose::Bool=false) = SDPAGMP(speed, executable, verbose)


immutable CSDPSolver <: SDPSolver
    executable::String
    verbose::Bool
end

CSDPSolver(; executable="csdp", verbose=false) = CSDPSolver(executable, verbose)

executable(solver::CSDPSolver) = solver.executable

type CSDPSDPModel <: AbstractSDPModel
    solver::CSDPSolver
    sdp
    sol
    function CSDPSDPModel(solver::CSDPSolver)
        new(solver, nothing, nothing)
    end
end

function loadproblem!(m::CSDPSDPModel, sdp::SparseSDP)
    m.sdp = sdp
end

SDPModel(solver::CSDPSolver) = CSDPSDPModel(solver)

MathProgBase.ConicModel(s::CSDPSolver) = SDPtoConicBridge(SDPModel(s))

function MathProgBase.optimize!(m::CSDPSDPModel)
    m.sol = solve(m.sdp, m.solver)
end
MathProgBase.status(m::CSDPSDPModel) = :Optimal # FIXME
MathProgBase.getsolution(m::CSDPSDPModel) = primalmatrix(m.sol)
MathProgBase.getobjval(m::CSDPSDPModel) = obj(m.sol)
MathProgBase.getvartype(m::CSDPSDPModel) = error("Not implemented yet")
MathProgBase.getobjbound(m::CSDPSDPModel) = error("Not implemented yet")
MathProgBase.getobjgap(m::CSDPSDPModel) = error("Not implemented yet")
MathProgBase.getrawsolver(m::CSDPSDPModel) = error("Not implemented yet")
MathProgBase.getsolvetime(m::CSDPSDPModel) = error("Not implemented yet")
MathProgBase.getsimplexiter(m::CSDPSDPModel) = error("Not implemented yet")
MathProgBase.getbarrieriter(m::CSDPSDPModel) = error("Not implemented yet")
MathProgBase.getnodecount(m::CSDPSDPModel) = error("Not implemented yet")


function solve{T<:Real}(sdp::SparseSDP{T}, solver::SDPAGEN; io::IO=STDOUT, outputfname="output.temp")
    datafname, dataio = mktemp()
    if verbose(solver)
        println(io, "Data filename: $datafname")
    end
    writesdpasparse(dataio, sdp)
    flush(dataio)
    primalobjective = NaN
    dualobjective = NaN
   
    for l in eachline(`$(executable(solver)) -ds $datafname -o $outputfname`)
        if verbose(solver)
            print(io, l)
        end
        if startswith(l, "objValPrimal = ")
            f = float(strip(split(l, " = ")[2]))
            primalobjective = ismaximizationproblem(sdp) ? f : -f
        end
        if startswith(l, "objValDual ")
            f = float(strip(split(l, " = ")[2]))
            dualobjective = ismaximizationproblem(sdp) ? f : -f
         end
    end
    close(dataio)
    SparseSDPSolution(primalobjective, dualobjective) 
end

function solve(sdp::SparseSDP, solver::CSDPSolver; io::IO=STDOUT, outputfname=nothing)
    sdp, cm, bm, ems = normalize(sdp)

    datafname, dataio = mktemp()
    removeoutputfname = false
    if outputfname == nothing
        outputfname, outputio = mktemp()
        removeoutputfname = true
    else
        outputio = open(outputfname, "w")
    end
    if verbose(solver)
        println(io, "Data filename: $datafname")
        println(io, "Output filename: $outputfname")
    end
    writesdpasparse(dataio, sdp)
    close(dataio)
    
    primalobjective = NaN
    dualobjective = NaN
    for l in eachline(`$(executable(solver)) $datafname $outputfname`)
        if verbose(solver)
            print(io, l)
        end
        if startswith(l, "Primal objective value: ")
            f = float(strip(split(l, ": ")[2]))
            primalobjective = ismaximizationproblem(sdp) ? f : -f
        end
        if startswith(l, "Dual objective value: ")
            f = float(strip(split(l, ": ")[2]))
            dualobjective = ismaximizationproblem(sdp) ? f : -f
        end
    end
    sol = SparseSDPSolution(primalobjective, dualobjective)
    readcsdpoutput!(outputio, sol, cm, bm, ems)
    close(outputio)
    if removeoutputfname 
        #rm(outputfname)
    end
    #rm(datafname)
    sol   
end
