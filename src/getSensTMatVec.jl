import jInv.ForwardShare.getSensTMatVec


function getSensTMatVec(v::Vector,m::Vector,pFor::FWIparam)

    # extract pointers
    Mesh  = pFor.Mesh
    omega = pFor.omega
    gamma = pFor.gamma
    Q     = pFor.Sources
    P     = pFor.Receivers
    LU    = pFor.Ainv
    U     = pFor.Fields

    nsrc = size(Q,2); nfreq = length(omega); nrec = size(P,2)

    # reshape v by receivers x sources x frequencies
    v    = reshape(v,nrec,nsrc,nfreq)

    # allocate space for result
    JTv  = zeros(length(m))

    # derivative of mass matrix
    An2cc = getNodalAverageMatrix(Mesh)
    dM    = repeat(1 .- 1im*vec(gamma),1,nsrc)
    for i=1:nfreq
        Lam  = LU[i]\(P*v[:,:,i])
        JTvi =  omega[i]^2*dM.*(An2cc*(U[:,:,i].*Lam))
        JTv +=  sum(real(JTvi),dims=2)
    end
    return vec(JTv)
end
