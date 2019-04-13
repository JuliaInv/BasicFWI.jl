include("setupFWItests.jl");

using jInvVis
using PyPlot
figure(20)
imshow(m');

figure(21);
imshow(reshape(gamma,(nx,nz))')

Dp,pFor = getData(vec(m),pFor,false)
nfreq = length(omega);

Dobs = fetch(Dp);
Wd = 1.0./(Dobs .+ 1e-2*maximum(abs(Dobs[:])));


# Wd   = Array{Array{Float64,2}}(undef,length(omega))
# Dobs = Array{Array{Float64,2}}(undef,length(omega))

# for k=1:nfreq
	# Dobs[k] = fetch(Dp[k]);
	# Wd[k] = 1.0./(Dobs[k] .+ 1e-2*maximum(abs(Dobs[k][:])));
	# figure(k);
	# imshow(Dobs[k])
# end
#||W_d*(D(m) - Dobs)||


N = prod(Mr.n);
Iact = sparse(I,N,N);
mback   = zeros(Float64,N);

########################################################################################################
##### Set up remote workers ############################################################################
########################################################################################################

Ainv = getJuliaSolver(); 
## Choose the workers for FWI (here, its a single worker)
misfun = SSDFun;
pMis = getMisfitParam(pFor, Wd, Dobs, misfun, Iact,mback);

