include("setupFWItests.jl")

using Test


println("\tTest getData serial and parallel")
Ds,pFor  = getData(vec(m),pFor)
Dp,pForp = getData(vec(m),pForp,M2Mp,true)
for k=1:length(Dp)
	@test norm(vec(Ds[:,:,k]-fetch(Dp[k])))/norm(vec(Ds[:,:,k])) < 1e-10
end

# check sensitivities
function dFun(m,pFor,v=[])
	D,pFor = getData(m,pFor)
	if isempty(v)
		return vec(D)
	else
		dD = getSensMatVec(v,m,pFor)
		return vec(D),vec(dD)
	end
end

println("\tderivative tests")
m0 = vec(m) + .001*randn(prod(n))
f(x,v=[]) = dFun(x,pFor,v)
chkDer1, = checkDerivative(f,m0,out=false)
@test chkDer1 == true

println("\tadjoint tests")
v1 = randn(length(Ds))
v2 = randn(length(m))

t1 = dot(v1,getSensMatVec(v2,vec(m),pFor))
t2 = dot(v2,getSensTMatVec(v1,vec(m),pFor))
@test norm(t1-t2)/norm(t1)<1e-8
