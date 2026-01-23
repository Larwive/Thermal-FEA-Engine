# Thermal-FEA-Engine

A C++ finite element solver for steady-state heat conduction problems.

## Motivation

This project was initiated as a personal learning exercise to strengthen my C++ skills through a concrete numerical problem.

The goal is to build a small but complete finite element pipeline: mesh definition, matrix assembly, boundary conditions and numerical solving, with an emphasis on code clarity and progressive design, without relying on external libraries.

## What it does

The engine solves a steady-state heat equation on a 2D triangular mesh:

- Nodes are defined by their coordinates
- Elements connect nodes into triangular cells
- Boundary temperatures (Dirichlet conditions) can be imposed
- The global stiffness matrix is assembled
- The temperature field is solved using a conjugate gradient method

## Example

For a square domain discretized into two triangular elements, with a temperature of 100°C imposed on one corner and 15°C on the opposite corner, the solver computes a smooth temperature gradient across the domain.

## Project status

The project is currently functional and supports:

- 2D triangular meshes
- Heat conduction with constant material properties
- Dirichlet boundary conditions
- A conjugate gradient solver

Current work focuses on:

- Improving usability and API clarity
- Adding sparse matrix support to handle larger meshes
- Extending boundary condition support

## Design choices

Implementation and design decisions are documented in [choices.md](choices.md).
