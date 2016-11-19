# wrapper to convert SDP solver into Conic solver

# To enable Conic support from an SDP solver, define, e.g.,
# ConicModel(s::CSDPSolver) = SDPtoConicBridge(SDPModel(s))

# This file is adapted from lpqp_to_conic.jl file from MathProgBase

type SDPtoConicBridge <: MathProgBase.AbstractConicModel
    sdpmodel::AbstractSDPModel
    varmap
    c
    A
    b
    constr_cones
    var_cones
end

# FIXME implements supportedcones

SDPtoConicBridge(m::AbstractSDPModel) = SDPtoConicBridge(m, nothing, nothing, nothing, nothing, nothing, nothing)

export SDPtoConicBridge

MathProgBase.numvar(m::SDPtoConicBridge) = size(m.A,2)
MathProgBase.numconstr(m::SDPtoConicBridge) = size(m.A,1)

function getmatdim(k)
    # n*(n+1)/2 = k
    # n^2+n-2k = 0
    # (-1 + sqrt(1 + 8k))/2
    n = (-1 + sqrt(1 + 8k)) / 2
    if n * (n+1) != 2*k
        error("sdp dim not consistent")
    end
    convert(Int, n)
end

# To transform Conic problems into SDP problems
function MathProgBase.loadproblem!(model::SDPtoConicBridge, c, A, b, constr_cones, var_cones)
    m, n = size(A)
    model.c = c
    model.A = A
    model.b = b
    model.constr_cones = constr_cones
    model.var_cones = var_cones

    # Conic form        LP form
    # min  c'x          min      c'x
    #  st b-Ax ∈ K_1     st lb <= Ax <= b
    #        x ∈ K_2         l <=  x <= u

    # If a cone is anything other than [:Free,:Zero,:NonNeg,:NonPos,:SOC,:SOCRotated,:SDP], give up.
    bad_cones = [:ExpPrimal, :ExpDual]
    for (cone,idxs) in var_cones
        cone in bad_cones && error("Cone type $(cone) not supported")
    end

    blk = 0
    # For a variable at column index `col' in the conic model,
    # varmap[col] gives an array such that each coefficient A[.,col] should be replaced by the sum,
    # over each element (blk, i, j, coef) of the array of
    # X[blk][i,j] * (A[.,col] * coef)
    # Where X[blk] is the blk'th block of X
    model.varmap = varmap = Vector{Vector{Tuple{Int,Int,Int,Float64}}}(n)
    for (cone,idxs) in var_cones
        # If a cone is anything other than [:Free,:Zero,:NonNeg,:NonPos,:SOC,:SOCRotated,:SDP], give up.
        if cone == :Free
            for i in idxs
                blk += 2
                # x free transformed into x = y - z with y, z >= 0
                varmap[i] = [(blk-1,1,1,1.), (blk,1,1,-1.)]
            end
        elseif cone == :Zero
            for i in idxs
                varmap[i] = []
            end
        elseif cone == :NonNeg
            for i in idxs
                blk += 1
                varmap[i] = [(blk,1,1,1.)]
            end
        elseif cone == :NonPos
            for i in idxs
                blk += 1
                varmap[i] = [(blk,1,1,-1.)]
            end
        elseif cone == :SOC
            error("not supported yet")
        elseif cone == :SOCRotated
            error("not supported yet")
        elseif cone == :SDP
            d = getmatdim(length(idxs))
            k = 0
            blk += 1
            for i in 1:d
                for j in i:d
                    k += 1
                    # In the MPB conic model, those are scaled by sqrt(2)
                    coef = i == j ? 1. : 1/sqrt(2)
                    varmap[idxs[k]] = [(blk,i,j,coef)]
                end
            end
        end
    end
    constr = 0
    # For the constraint at row index `row' in the conic model,
    # constrmap[row] gives the index of the constraint in the SDP model,
    # a value of 0 meaning that it does not correspond to any constraint
    constrmap = Vector{Int}(m)
    # slackmap[row] gives (blk,i,j,coef) indicating that a slack variable has been created at X[blk][i,j] with coefficient coef
    # blk=0 corresponds to no slack
    slackmap = Vector{Tuple{Int,Int,Int,Float64}}(m)
    for (cone,idxs) in constr_cones
        if cone == :Free
            constrmap[idxs] = 0
            slackmap[idxs] = 0
        elseif cone == :SOC
            error("not supported yet")
        elseif cone == :SOCRotated
            error("not supported yet")
        else
            for idx in idxs
                constr += 1
                constrmap[idx] = constr
            end
            if cone == :Zero
                slackmap[idxs] = (0,0,0,0.)
            elseif cone == :NonNeg
                for idx in idxs
                    blk += 1
                    slackmap[idx] = (blk,1,1,-1.)
                end
            elseif cone == :NonPos
                for idx in idxs
                    blk += 1
                    slackmap[idx] = (blk,1,1,1.)
                end
            elseif cone == :SDP
                d = getmatdim(length(idxs))
                k = 0
                blk += 1
                for i in 1:d
                    for j in i:d
                        k += 1
                        slackmap[idxs[k]] = (blk,i,j,-1.)
                    end
                end
            end
        end
    end

    # Writing the sparse block diagonal matrices
    sdp = SparseSDP(maximize=false)
    for row in 1:m
        if constrmap[row] != 0
            setrhs!(sdp, constrmap[row], b[row])
        end
        blk, i, j, coef = slackmap[row]
        if blk != 0
            setcon!(sdp, constrmap[row], blk, i, j, coef)
        end
    end
    rows = rowvals(A)
    vals = nonzeros(A)
    for col in 1:n
        for k in nzrange(A, col)
            row = rows[k]
            if constrmap[row] != 0 # Free constraint
                val = vals[k]
                for (blk, i, j, coef) in varmap[col]
                    setcon!(sdp, constrmap[row], blk, i, j, val*coef)
                end
            end
        end
    end
    for col in 1:n
        if c[col] != 0
            for (blk, i, j, coef) in varmap[col]
                setobj!(sdp, blk, i, j, coef*c[col])
            end
        end
    end

    loadproblem!(model.sdpmodel, sdp)
end

for f in [:optimize!, :status, :getobjval, :getvartype]
    @eval MathProgBase.$f(model::SDPtoConicBridge) = MathProgBase.$f(model.sdpmodel)
end

function MathProgBase.getdual(model::SDPtoConicBridge)
    error("Not implemented yet")
end

function MathProgBase.getsolution(model::SDPtoConicBridge)
    X = MathProgBase.getsolution(model.sdpmodel)
    n = size(model.A, 2)
    x = zeros(Float64, n)
    for col in 1:n
        for (blk, i, j, coef) in model.varmap[col]
            x[col] += X[blk][i,j] / coef
        end
    end
    x
end

function setvartype!(model::SDPtoConicBridge, vtype)
    if any(t->t != :Cont, vtype)
        error("Invalid variable type")
    end
end

for f in MathProgBase.SolverInterface.methods_by_tag[:rewrap]
    @eval MathProgBase.$f(model::SDPtoConicBridge) = MathProgBase.$f(model.sdpmodel)
end
