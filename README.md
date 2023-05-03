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

# Data submission

## Data organization
Groups should submit time-histories for each case (Config,Motion,Order,Space-Res,Time-Res) in separate files following the file location and name convention:
```
<GroupName>/<Case>-M<MotionNumber>-h<GridIndex>-p<OrderIndex>-t<TimeIndex>.txt
```

An example data file location and name would be
```
AFRL/Airfoil-M1-h0-p0-t0.txt
```

## Data format
For each contributed Case/Motion/Resolution, we are requesting time-series data for a set of outputs. Time-integrated quantities will be computed in data-processing by the organizers. Time-histories should include the time-history bounds (i.e. data at initial time t=0 and also data at final time t=1,2, or 40 depending on the test case). Time-histories should include the time-value for each time-instance as well as the requested outputs at each time-instance (Outputs are described in Eqns. 14-17 in the test suite document). Each contributed data-file (representative of a particular Case/Motion/Resolution) should be submitted in comma-separated-value format that consists of a single-line header and time-series data on subsequent lines. Time-series data should be provided with at least 8-digits of precision. If a requested output is not able to be provided the entry should be filled with the junk-value 12345678.

The data-header should be the following:
```
Time, Y-Force, Work integrand, Mass, Mass error
```

An example data file contents for a submission that does not provide 'Mass error' would be:
```
Time, Y-Force, Work integrand, Mass, Mass error
0.0000000, 1.5438375, 3.4932846, 3.0829579, 12345678.
0.2000000, 1.5648394, 3.5349762, 3.0830752, 12345678.
0.4000000, 1.5740924, 3.8028847, 3.0840783, 12345678.
0.6000000, 1.5638740, 3.4397543, 3.0892051, 12345678.
0.8000000, 1.5503957, 3.4932846, 3.0913753, 12345678.
1.0000000, 1.5400933, 3.4932846, 3.0940148, 12345678.
```


# Working group notes

## 5 April, 2023

**Notes:**
- MIT group showed preliminary results of cylinder motion.
- James(MIT) discussed appropriate initial conditions for incompressible cylinder case.
  
**Action items:**
- Per generating + providing cylinder meshes in Gmsh format.
- Nathan working towards adding 2022 data and setting up standardized data-input/processing.
- James adding write-up substantiating incompressible initial conditions.

