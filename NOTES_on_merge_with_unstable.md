### Changes to unstable:
* We don't need MatsubaraIndex anymore => delete struct and related functions.
* The current structure of the meshes should allow operations (+,-) for both AbstractMeshPoints, AbstractValue and mixed.
  I also removed the return type of these operations since e.g. typeof(x + x) = MatsubaraFrequency{Boson} for typeof(x) = MatsubaraFrequency{Fermion}.

### Suggestions:
* change function names reciprocals(...) --> bz_indices(...)     euclideans(...) --> bz_euclideans(...)
* The code looks a bit strange when we iterate over mesh.points and need to obtain the underlying value by calling value(value(...)). It's just a cosmetic thing. I like the construction. I introduced the function    'plain_value(...)' which always returns the plain underlying value of a mesh point.