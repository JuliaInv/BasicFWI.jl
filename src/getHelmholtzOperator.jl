export getHelmholtzOperator, getABL

"""
	H = getHelmholtzOperator(m,gamma,omega,Mesh)

	build Helmholtz operator
		H = Laplace - omega^2 * MassMatrix(m,gamma)

	Input:
		m                  - model
		gamma              - attenuation
		omega::Number      - frequency
		Mesh::RegularMesh  - mesh

	Output
		H::SparseMatrixCSC
"""
function getHelmholtzOperator(m,gamma,omega::Number,Mesh::RegularMesh)

Lap   = getNodalLaplacianMatrix(Mesh)
M     = getMassMatrix(vec(m),vec(gamma),Mesh)
# Get the Helmholtz operator (note the sign)
H = Lap .- omega^2 * M
return H
end


function getABL(Msh::RegularMesh,NeumannAtFirstDim::Bool,ABLpad::Array{Int64},ABLamp::Float64)
  h = Msh.h;
  n = Msh.n;
  pad = ABLpad;
  ntup = tuple(n...);
  
  if Msh.dim==2
	# x1 = linspace(-1,1,n[1]);
	# x2 = linspace(0,1,n[2]);
	# X1,X2 = ndgrid(x1,x2);
	# padx1 = ABLpad[1];
	# padx2 = ABLpad[2];
	# gammaxL = (X1 - x1[padx1]).^2;
	# gammaxL[padx1+1:end,:] = 0
	# gammaxR = (X1 - x1[end-padx1+1]).^2
	# gammaxR[1:end-padx1,:] = 0

	# gammax = gammaxL + gammaxR
	# gammax = gammax/maximum(gammax);

	# gammaz = (X2 - x2[end-padx2+1]).^2
	# gammaz[:,1:end-padx2] = 0
	# gammaz = gammaz/maximum(gammaz);
	
	# gamma = gammax + gammaz
	# gamma *= ABLamp;
	# gamma[gamma.>=ABLamp] = ABLamp;

	gamma = zeros(ntup);
	b_bwd1 = ((pad[1]:-1:1).^2)./pad[1]^2;
	b_bwd2 = ((pad[2]:-1:1).^2)./pad[2]^2;
  
	b_fwd1 = ((1:pad[1]).^2)./pad[1]^2;
	b_fwd2 = ((1:pad[2]).^2)./pad[2]^2;
	I1 = (n[1] - pad[1] + 1):n[1];
	I2 = (n[2] - pad[2] + 1):n[2];
  
	if NeumannAtFirstDim==false
		gamma[:,1:pad[2]] += ones(n[1],1)*b_bwd2';
		gamma[1:pad[1],1:pad[2]] -= b_bwd1*b_bwd2';
		gamma[I1,1:pad[2]] -= b_fwd1*b_bwd2';
	end

	gamma[:,I2] +=  ones(n[1],1)*b_fwd2';
	gamma[1:pad[1],:] += b_bwd1*ones(1,n[2]);
	gamma[I1,:] += b_fwd1*ones(1,n[2]);
	gamma[1:pad[1],I2] -= b_bwd1*b_fwd2';
	gamma[I1,I2] -= b_fwd1*b_fwd2';
	gamma *= ABLamp;
	# figure()
	# imshow(gamma'); colorbar()
	
  else
	x1 = linspace(-1,1,n[1]);
	x2 = linspace(-1,1,n[2]);
	x3 = linspace( 0,1,n[3]);
	X1,X2,X3 = ndgrid(x1,x2,x3);
	padx1 = ABLpad[1];
	padx2 = ABLpad[2];
	padx3 = ABLpad[3];
	gammaL = (X1 - x1[padx1]).^2;
	gammaL[padx1+1:end,:,:] = 0.0
	gammaR = (X1 - x1[end-padx1+1]).^2
	gammaR[1:end-padx1,:,:] = 0.0
	
	gammat = gammaL + gammaR;
	gammat = gammat/maximum(gammat);
	gamma = copy(gammat);
	gammat[:] = 0.0;

	gammaL = (X2 - x2[padx2]).^2;
	gammaL[:,padx2+1:end,:] = 0.0
	gammaR = (X2 - x2[end-padx2+1]).^2
	gammaR[:,1:end-padx2,:] = 0.0
	
	gammat = gammaL + gammaR
	gammat = gammat/maximum(gammat);
	gamma += gammat;

	gammat = (X3 - x3[end-padx3+1]).^2
	gammat[:,:,1:end-padx3] = 0.0
	gammat = gammat/maximum(gammat);
	gamma += gammat;
	gamma *= ABLamp;
	gamma[gamma.>=ABLamp] = ABLamp;
  end
  return gamma;
end