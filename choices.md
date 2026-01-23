# Design and implementation choices

## Overall goals

The primary goal of this project is educational: building a small but complete finite element pipeline in C++, with an emphasis on clarity, explicit data flow, and progressive design.

The code favors readability and separation of concerns over performance optimizations. All components are written without external dependencies to better understand the underlying numerical and software concepts.

## Mesh representation

### Node

A Node represents a geometric point in space and stores only:

- a unique identifier
- its spatial coordinates

Nodes do not store any connectivity information or physical values (like temperature).
This choice keeps the geometric description independent from the physics being solved.

### Element

An Element stores only the indices of its nodes, not references or copies of Node objects.

This indirection was chosen to:

- avoid data duplication
- keep elements lightweight
- allow the same mesh structure to be reused with different physical problems

Element orientation is assumed to be consistent and is handled explicitly when assembling local matrices.

### Mesh

The Mesh class owns the nodes and elements and acts as a central access point.

It is intentionally kept generic and physics-agnostic:

- no temperature field
- no material properties
- no solver logic

This allows the same mesh representation to be reused for other finite element problems (elasticity, diffusion).

## Separation of responsibilities

### Local physics vs global assembly

The computation of local element matrices (Ke) is separated from the global assembly process.

- Local matrices are computed in the physics layer, based only on geometry and material parameters.
- Assembly is handled in the solver layer, which maps local contributions into the global system.

This separation makes the physical formulation easier to test and modify independently of the global solver.

### Solver independence

The linear solver operates only on abstract Matrix and Vector types and does not depend on mesh or physics details.

This design allows:

- experimenting with different solvers
- changing the physical model without modifying the solver code

## Linear algebra layer

A minimal linear algebra layer (Matrix and Vector) was implemented to avoid external dependencies and to better understand the numerical algorithms involved.

The current implementation uses dense storage for simplicity.

## Boundary conditions handling

Dirichlet boundary conditions are imposed using a classical elimination method.

The contribution of constrained degrees of freedom is first transferred to the right-hand side before modifying the global matrix.
This ensures the resulting linear system remains consistent and solvable.
