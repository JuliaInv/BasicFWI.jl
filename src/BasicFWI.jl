module BasicFWI
using jInv.Mesh
using jInv.Utils
using Distributed
using SparseArrays
using LinearAlgebra

export FWIparam
import jInv.ForwardShare.ForwardProbType
"""
	FWIparam <: ForwardProbType

	type defining Full Waveform Inversion (FWI) problem.

	Constructor: getFWIparam(omega,gamma,Q,P,Mesh)

	Fields:
		omega     - frequencies
		gamma     - attenuation
		Sources   - matrix of sources
		Receivers - receiver matrix
		Mesh      - forward mesh
		Fields    -
		Ainv      - for linear solver, LU factorizations
"""
mutable struct FWIparam <: ForwardProbType
    omega     # frequencies
    gamma     # attenuation
    Sources   # Sources
    Receivers # Receivers
    Mesh      # Mesh
	Fields    # Fields
	Ainv      # LU factorization
end

export getFWIparam
"""
	pFor = getFWIparam(omega,gamma,Q,P,Mesh,doDistribute=false)

	constructs FWI param

"""
function getFWIparam(omega,gamma,Q,P,Mesh,doDistribute=false)
    nfreq = length(omega)

    if doDistribute && (nfreq > 1)
        pFor     = Array{RemoteChannel}(undef, nfreq)
        probsMax = ceil(Integer,nfreq/nworkers())
        nprobs   = zeros(maximum(workers()))

        i=1; nextidx() = (idx = i; i+=1; idx)

        @sync begin
            for p=workers()
                @async begin
                    while true
                        if nprobs[p]>=probsMax
                            break
                        end
                        idx = nextidx()
                        if (idx>nfreq)
                            break
                        end
                        pFor[idx] = initRemoteChannel(getFWIparam,p,omega[idx],gamma,Q,P,Mesh)
                        nprobs[p] +=1
                        wait(pFor[idx])
                    end
                end
            end
        end
    else
        pFor = FWIparam(omega,gamma,Q,P,Mesh,[],[])
    end
    return pFor
end

include("getData.jl")
include("getSensMatVec.jl")
include("getSensTMatVec.jl")
include("getHelmholtzOperator.jl")
include("getMassMatrix.jl")

end
