# 2024 Mesh-Motion
This repository is a resource for the Mesh Motion test suite working group in the 2024 High-Fidelity CFD Verification Workshop.


# Cloning the repository

Note, the 2022 workshop repository has been added as a submodule. To clone the 2024 repository along with the 2022 repo submodule:
```
git clone --recurse-submodules https://github.com/HighFidelityCFDVerificationWorkshop/2024MeshMotion 
```

If you have already cloned the repository, you can retrieve the 2022 repository submodule using:
```
git submodule update --init --recursive
```

# TODO
- [x] Finalize test suite motions and test cases ([SciTech 2024 paper](https://highfidelitycfdverificationworkshop.github.io/papers/mesh_motion.pdf))
- [x] Provide mesh motion function implementations and their derivatives in C, Python, and Fortran.
- [ ] Generate family of cylinder meshes (Per)
- [ ] Set up input format + data-processing (Nathan)
- [ ] Document appropriate initial conditions for incompressible tests (James)
- [ ] Preliminary results (All)


# Working group notes

## 5 April, 2023

**Notes:**
- MIT group showed preliminary results of cylinder motion.
- James(MIT) discussed appropriate initial conditions for incompressible cylinder case.
  
**Action items:**
- Per generating + providing cylinder meshes in Gmsh format.
- Nathan working towards adding 2022 data and setting up standardized data-input/processing.
- James adding write-up substantiating incompressible initial conditions.

