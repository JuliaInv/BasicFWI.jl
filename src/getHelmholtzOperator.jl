export getHelmholtzOperator

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
