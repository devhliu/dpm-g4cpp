# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Instructions

### Prerequisites
- Geant4 version >= 10.7.p01 must be installed on the system
- CMake 3.1+ required

### Build Commands
```bash
cd ../  # Navigate to the root dpm-g4cpp directory
mkdir build
cd build/
cmake ../ -DGeant4_DIR=<path_to_geant4_installation>
make
```

This creates two executables:
- `dpm_GenerateData`: Generates simulation data files (requires Geant4)
- `dpm_Simulate`: Runs the dose simulation using pre-generated data
- `test_brem`: Bremsstrahlung testing utility
- `modeltests`: Model validation tests using the validation library

### Running Simulations

1. **Generate Data** (one-time setup per configuration):
```bash
./dpm_GenerateData -c <config_index> -d <output_data_dir>
```

2. **Run Simulation**:
```bash
./dpm_Simulate -n <num_histories> -e <energy_MeV> -d <data_dir> -c <config_index>
```

3. **Run Tests**:
```bash
./modeltests
```

Common options:
- `-p/--primary-particle`: e- or gamma (default: e-)
- `-e/--primary-energy`: Energy in MeV (default: 0.1)
- `-n/--number-of-histories`: Number of primary events (default: 1.0E+5)
- `-d/--input-data-dir`: Data directory (default: ./data)
- `-b/--voxel-size`: Voxel size in mm (default: 1.0)
- `-c/--configuration-index`: Pre-defined geometry configuration (0-3)

## Architecture Overview

The project is divided into two main phases:

### 1. Data Generation Phase (`DataInit/`)
- Generates all physics data needed for simulation
- Relies on Geant4 for physics modeling
- Outputs data files to be used by the simulation phase

### 2. Simulation Phase (`Simulation/`)
- Pure C++ implementation without Geant4 dependency
- Optimized for deterministic execution (no rejection loops)
- Suitable for GPU/FPGA implementation

### Key Components

#### Core Simulation Classes
- **Geom**: Voxelized geometry representation with pre-defined configurations
  - Supports 4 pre-defined geometries (index 0-3)
  - Voxel-based navigation with fast boundary computation

- **Track**: Particle state container
  - Position, direction, energy, type (e-/e+/gamma)
  - Current voxel indices and material

- **TrackStack**: Global particle stack for managing primary and secondary particles

#### Physics Data Classes
- **SimMaterialData**: Material properties (density, cuts, names)
- **SimElectronData**: All electron/positron physics data
  - IMFP data for elastic, Moller, and bremsstrahlung interactions
  - Stopping power data
  - Sampling tables for energy and angular distributions
- **SimPhotonData**: Photon physics data
  - IMFP for photon interactions
  - Klein-Nishina sampling tables

#### Physics Interactions
- Elastic scattering (Goudsmit-Saunderson)
- Moller scattering (e-e interactions)
- Bremsstrahlung (with Seltzer-Berger DCS)
- Photon interactions (photoelectric, Compton, pair production)

### Geometry Configurations
The simulation supports pre-defined geometries:
- **0**: Homogeneous material (index 0)
- **1**: Two-layer slab (0-1cm: mat0, 1-3cm: mat1, >3cm: mat0)
- **2**: Two-layer slab (0-2cm: mat0, 2-4cm: mat1, >4cm: mat0)
- **3**: Three-layer slab (Water-Titanium-Bone-Water)

### Data Flow
1. Data generation creates physics tables for each material
2. Simulation loads these tables at runtime
3. Particles are transported voxel-by-voxel
4. Interactions are sampled using pre-computed tables (rejection-free)
5. Energy deposition is scored per voxel

## Important Implementation Details

- No dynamic memory allocation during simulation (all pre-allocated)
- Deterministic algorithms (no while loops with random conditions)
- All physics interactions use alias sampling for rejection-free sampling
- Units: Energy in MeV, distance in mm, density in g/cm³
- Voxel size determines geometry granularity and affects performance
- Track-stack based architecture for particle management
- Continuous tracking loop with distance-to-boundary computation

## Directory Structure

```
dpm-g4cpp/
├── CMakeLists.txt           # Main build configuration
├── dpm_GenerateData.cc      # Data generation main entry
├── dpm_Simulate.cc          # Simulation main entry
├── Configuration.hh         # Global configuration constants
├── DataInit/                # Data generation phase
│   ├── inc/                 # Headers for data generation
│   └── src/                 # Sources for data generation
├── Simulation/              # Simulation phase (this directory)
│   ├── inc/                 # Headers for simulation classes
│   └── src/                 # Sources for simulation classes
├── Utils/                   # Shared utilities
│   ├── inc/                 # Utility headers
│   └── src/                 # Utility sources
├── SBTables/                # Seltzer-Berger DCS tables
├── tests/                   # Test configurations
├── modelValidationTests/    # Validation framework
└── docs/                    # Architecture documentation
    └── dpm_architecture.md  # Detailed physics architecture
```

## Additional Resources

- `docs/dpm_architecture.md`: Comprehensive physics and transport architecture documentation
- `README.md`: Project overview and motivation
- `tests/WaterTiBoneWater_15MeV/`: Reference test data and results