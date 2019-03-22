import jInv.ForwardShare.getSensMatVec

function getSensMatVec(v::Vector,m::Vector,pFor::FWIparam)

    # extract pointers
    Mesh  = pFor.Mesh
    omega = pFor.omega
    gamma = pFor.gamma
    Q     = pFor.Sources
    P     = pFor.Receivers
    U     = pFor.Fields
    LU    = pFor.Ainv

    nsrc = size(Q,2)
    nfreq = length(omega)

    # allocate space for matvec product
    Jv   = zeros(size(P,2),nsrc,nfreq)

    # derivative of mass matrix
    An2cc = getNodalAverageMatrix(Mesh)
    dM   = repeat(An2cc'*((1 .- 1im*vec(gamma)).*v),1,nsrc)
    for i=1:nfreq
        R   = U[:,:,i].*dM
        Lam = LU[i]\R
        Jv[:,:,i] = real(omega[i]^2*P'*Lam)
    end

    return vec(Jv)
end
