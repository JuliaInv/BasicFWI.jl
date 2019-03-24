using BasicFWI
using jInv.Mesh
using jInv.Utils
using jInv.ForwardShare

if nworkers()==1
	addprocs(2)
end


@everywhere begin
	using BasicFWI
	using jInv.Mesh
	using jInv.Utils
	using jInv.ForwardShare
end

# setup model and attenuation function
nx = 96; nz = 48
domain = [-1.,1,0,1]
Mr = getRegularMesh(domain,[nx,nz])

# the slowness model
m = ones(nx,nz);
m[20:70,10:32] .= 3;


# attenuation for BC
padx = 6; padz = 6
a    = 20;
xc = getCellCenteredGrid(Mr)
gammaxL = a*(xc[:,1] .- xc[padx,1]).^2;
gammaxL[padx+1:end,:] .= 0
gammaxR = a*(xc[:,1] .- xc[end-padx+1,1]).^2
gammaxR[1:end-padx,:] .= 0

gammax = gammaxL + gammaxR

gammaz = a*(xc[:,2] .- xc[end-padz+1,2]).^2
gammaz[:,1:end-padz] .= 0

gamma = gammax + gammaz

# parameters for the Helmholtz (units in km)
h = [50.0; 50.0;]/1000
n = [nx; nz;]
omega = 2*pi*[0.2;0.5;1;1.2;1.4;1.5;2.0;3.0]
nfreq = length(omega)

# generate sources
q = zeros(tuple(n.+1...)); q[padx+1:2:end-padx-1,1] .= 1e4
Q = sdiag(vec(q))
Q = Q[:,(LinearIndices(sum(Q,dims=2) .!= 0))[findall(sum(Q,dims=2) .!= 0)]]
nsrc = size(Q,2)

# receivers
p = zeros(tuple(n.+1...)); p[padx+1:end-padx-1,1] .= 1
P = sdiag(vec(p))
P = P[:,(LinearIndices(sum(P,dims=2) .!= 0))[findall(sum(P,dims=2) .!= 0)]]
nrec = size(P,2)

pFor  = getFWIparam(omega,gamma,Q,P,Mr)
pForp = getFWIparam(omega,gamma,Q,P,Mr,true)

# inversion mesh and forward mesh are the same here
M2Mp = ones(length(pForp))
