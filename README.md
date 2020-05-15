# OpenFOAM solver: rhoCentralFoam for CO2
The rhoCentralFoam solver is coupled with the CO2 look-up table. The coupling is not based on the object-oriented way. It may be improved in the future.

---
## Files
- ***Main solver***: CO2FoamSchemeV1
- ***2 cases***: 1dshocktube, ejector
- ***Other related modules***: thermophysicalModelsCO2, TurbulenceModelCO2. This 2 modules are almost the same as the original ones.
- ***compute y+***: compressiblePlusPostRANS

---
## How to compile
- Compile CO2 table in CO2FoamSchemeV1/CO2 to create .a library
- Compile related libraries: thermophysicalModelCO2 and TurbulenceModelCO2
- Compile main solver with external library (need openblas library)
- When compiling main solver, you need to choose corresponding BC and initialization

## How to use
- Use gmsh to create a mesh with physical boundary, use gmshToFoam to convert mesh
- Change constant/polymesh/boundary
- Set system/ file
- Create 0 solution
- Initialise with setFields
- If interpolation use mapFields
- call solver name to run 
