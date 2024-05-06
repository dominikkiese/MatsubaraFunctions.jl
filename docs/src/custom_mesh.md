# Custom mesh implementation

To implement a new mesh type `NewMesh <: AbstractMesh`, users need to provide a new value type 
`NewValue <: AbstractValue` as well as a `NewDomain <: AbstractDomain` which encodes implementation details.
As an example, consider the `IndexMesh` type, for which we have

| Abstract type            | Concrete type                            |
|--------------------------|------------------------------------------|
| AbstractMesh             | Mesh{MeshPoint{Index}, IndexDomain}      |
| AbstractValue            | Index                                    |
| AbstractDomain           | IndexDomain                              |

For `NewMesh` to be operational, the following functions need to be implemented/overloaded as well:
 * one or multiple outer constructors
 * `value`, for returning the bare value (e.g. the temperature dependent value of a Matsubara frequency) 
 * `is_inbounds`, for bounds checking
 * `save!`, `load_mesh`, for I/O to HDF5 files
 * `mesh_index`, for mapping `NewValue` instances to linear mesh indices
 * `mesh_index_bc`: for mapping `NewValue` instances to linear mesh indices (respecting boundary conditions)