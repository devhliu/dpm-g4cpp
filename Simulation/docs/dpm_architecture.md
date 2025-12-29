# DPM-g4cpp Physics Architecture Analysis

## Executive Summary

The DPM-g4cpp (Dose Planning Method in C++) is a highly optimized Monte Carlo simulation framework specifically designed for radiotherapy dose calculations. It implements deterministic, rejection-free particle transport algorithms that make it particularly suitable for GPU/FPGA implementations. The framework separates physics data generation from simulation execution, enabling fast dose computations with accurate modeling of photon and electron/positron interactions in voxelized geometries.

---

## 1. Particle Transport Architecture

### 1.1 Core Transport Framework

The simulation uses a **track-stack based architecture** where particles are transported one at a time in a **continuous tracking loop**:

1. **Track Management**:
   - `Track` structure stores complete particle state (position, direction, energy, type, voxel indices)
   - `TrackStack` (singleton) manages primary and secondary particles
   - particles are popped from stack, tracked to completion/death, then next particle processed

2. **Voxelized Geometry Navigation**:
   - Cubic voxels with configurable size (typically 1 mm)
   - Fast boundary computation using analytical geometry
   - Material indices assigned per voxel
   - Predefined geometries support homogeneous to multi-layer slab configurations

3. **Transport Algorithm**:
   - **Distance-to-boundary** computation determines maximum step
   - **Multiple interaction competition** between physics processes
   - **Stepping logic** with energy loss computation
   - **Boundary crossing** handling

#### 1.1.1 Track Stack Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    TrackStack (Singleton)                   │
├─────────────────────────────────────────────────────────────┤
│  Primary Track ──► Transport to completion                  │
│       │                │                                    │
│       ▼                ▼                                    │
│  [Pop from stack]  Track dies/leaves/energy < cut          │
│       │                │                                    │
│       └────────────────┼────────────────────────────────────┘
│                        ▼                                    │
│              Check secondaries generated?                   │
│                        │                                    │
│             Yes ┌──────┴──────┐                             │
│              │   Vacant?   ◄│───┐                          │
│              └──────┬──────┘   │                          │
│     [Push to stack]    Push     │                          │
│                     secondaries │                          │
└───────────────────────────────────┴──────────────────────────┘
```

#### 1.1.2 Voxel Navigation

```
           Z-axis (beam direction)
                ↑
                │
  ┌─────┬─────┬─────┬─────┬─────┐
  │  0  │  1  │  2  │  3  │  4  │ ← Voxels (iz index)
  ├─────┼─────┼─────┼─────┼─────┤
  │ H₂O │ Ti  │Bone │ H₂O │Vac  │ ← Material indices
  └─────┴─────┴─────┴─────┴─────┘
                │
     Particle path ◄───────────────── Y-axis
                │
     └─► Distance to boundary computation
         Fast: Δz = lbox - (z % lbox)
```

### 1.2 Deterministic Design Philosophy

Key innovation: **elimination of all rejection loops** through:
- Pre-computed sampling tables for all distributions
- Fixed number of random numbers per interaction
- Deterministic step size calculations
- No while-loops with random conditions

---

## 2. Photon Physics Implementation

### 2.1 Photon Transport Algorithm

**`KeepTrackingPhoton()`** implements photon transport with:

1. **Step Sampling**: Uses global maximum cross-section for delta-tracking
   ```
   step = -1/σ_max * ln(random)
   ```

2. **Interaction Selection**: Cumulative probability method:
   - Photoelectric absorption
   - Compton scattering
   - Pair production (E > 2mc²)

3. **Delta Interaction**: When no real interaction occurs

### 2.2 Photon Interactions

#### 2.2.1 Compton Scattering (Klein-Nishina)
- **Energy Sampling**: Uses pre-computed `SimKNTables`
- **3 random numbers**: For alias sampling without rejection
- **Kinematics**:
  - Energy fraction ε sampled from transformed distribution
  - Scattering angle: cosθ = 1 - (1-ε)/(εκ) where κ = Eγ/mc²
  - Azimuthal angle: φ = 2π × random

#### 2.2.2 Photoelectric Effect
- Complete photon energy absorption
- Energy deposited locally
- No secondary electron tracking (below cut)

#### 2.2.3 Pair Production
- **Energy sharing**: Uniform between e⁻ and e⁺
- **Total available energy**: Eγ - 2mc²
- **Secondary particles**: Both e⁻ and e⁺ created if above cut
- **Direction**: Same as incident photon (simplified)

---

## 3. Electron/Positron Physics Implementation

### 3.1 Electron Transport Algorithm

**`KeepTrackingElectron()`** implements complex electron transport with three competing interactions:

1. **Multiple Coulomb Scattering** (elastic)
2. **Moller/Bhabha Scattering** (inelastic)
3. **Bremsstrahlung** (radiative)

### 3.2 Multiple Coulomb Scattering (MCS)

#### 3.2.1 DPM Two-Step Method
Innovative approach avoiding continuous scattering:

1. **MSC Step Definition**:
   - Maximum step length: S_max(E) defined by parameters:
     - `s_low`: Minimum step (5 mm typical)
     - `s_high`: Maximum step (10 mm typical)
     - `e_cross`: Energy crossover (12 MeV typical)
   - Sigmoid-like function determines S_max vs energy

2. **Two-Step Transport**:
   - **Hinge point** at fraction of S_max
   - Angular deflection applied at hinge
   - Straight-line transport between deflections

3. **Goudsmit-Saunderson Sampling**:
   - Uses `SimGSTables` for rejection-free angular sampling
   - Based on numerical inversion of cumulative distribution
   - Reference material scaling for other materials

#### 3.2.2 First Transport Mean Free Path
- **tr1-mfp**: Specialized for DPM step control
- **Maximum scattering strength**: K₁(E) = S_max(E)/tr1-mfp(E)
- Units: dimensionless scaling parameter

#### 3.2.3 Two-Step MSC Visualization

```
DPM Multiple Scattering Implementation
────────────────────────────────────────────────────────────

Traditional Continuous MSC            DPM Two-Step MSC
───────────────────────              ─────────────────────
┌─────────────────────┐              ┌─────────────────────┐
│   • ──► • ──► •    │   ← Many    │ • ──► •            │ ← Step 1: to hinge
│   • ◄─► • ◄─► •    │   small    │   ↗                 │    (apply deflection)
│   • ──► • ──► •    │   steps    │  ○                  │
│   • ◄─► • ◄─► •    │            │ │                  │
│   • ──► • ──► •    │            │ ○                  │ ← Step 2: to end
└─────────────────────┘            │ ↘                 │    (no deflection)
                                   │   •                │
                                   └─────────────────────┘

Where:
• = particle position
→ = straight line transport
↗↘ = angular deflection at hinge
```

#### 3.2.4 S_max(E) Function Shape

```
S_max(E) vs Energy
      ↑
  S_h │───────────────┐
      │               │
      │               │
      │         . . .─┤───── Sigmoid transition
      │     . .
  S_l │───. .           │
      │
      └─────────────────┼─────→ Energy (log scale)
                E_cross
```

### 3.3 Moller Scattering (e⁻e⁻) / Bhabha (e⁺e⁻)

#### 3.3.1 Key Simplifications (DPM-style)
- **No energy dependence** of IMFP (approximate)
- **Material scaling** via Z/A ratio
- **Fast sampling** through `SimMollerTables`

#### 3.3.2 Energy Transfer Sampling
- **Alias method** with 3 random numbers
- **Kinematics**: Proper two-body final state
- **Secondary electron**: Tracked if above cut

### 3.4 Bremsstrahlung

#### 3.4.1 Implementation Details
- **Seltzer-Berger DCS** (differential cross section)
- **Material-dependent** IMFP and sampling
- **Energy sampling** via `SimSBTables`

#### 3.4.2 Angular Distribution
- **Simplified**: Forward peaked approximation
- **Photon direction**: Same as electron (performance optimization)

### 3.5 Energy Loss

#### 3.5.1 Continuous Slow Down Approximation
- **Restricted stopping power** from `SimStoppingPower`
- **Mid-point energy** evaluation for step
- **Linear approximation**: dE = s × dE/dx

#### 3.5.2 Threshold Tracking
- **Secondary production cut**: Typical 200 keV for e⁻
- **Below cut**: Energy deposited locally
- **Positron annihilation**: At end of track (2 × 511 keV photons)

---

## 4. Data Structures and Sampling Methods

### 4.1 Sampling Table Architecture

#### 4.1.1 Linear Alias Tables (`SimLinAliasData`)
- **Alias method** for O(1) sampling
- **Structure**:
  - Discrete values of the random variable
  - Alias probabilities
  - Alternative indices

#### 4.1.2 Energy Grid Structure
- **Logarithmic grid** for primary energies
- **Interpolation** between grid points
- **Range**: 0.1 keV to ~20 MeV

#### 4.1.3 Data Structure Hierarchies

```
Simulation Data Architecture
────────────────────────────────────────────────────────────

┌─────────────────────────────────────────────────────────────┐
│                Simulation Data Loaders                      │
├─────────────────────────────────────────────────────────────┤
│  SimElectronData          SimPhotonData                     │
│  ├─ fITr1MFPElastic        ├─ fIMFPTotal                   │
│  ├─ fMaxScatStrength       ├─ fIMFPCompton                 │
│  ├─ fIMFPMoller            ├─ fIMFPPairProd                │
│  ├─ fIMFPBrem              ├─ fIMFPPhoton                  │
│  ├─ fDEDX                  ├─ fTheKNTables                 │
│  └─ Sampling Tables:       └─ Sampling Tables:             │
│     ├─ fTheGSTables           └─ (none for photons)         │
│     ├─ fTheSBTables                                           │
│     └─ fTheMollerTables                                       │
└─────────────────────────────────────────────────────────────┘
                    │
                    ▼
┌─────────────────────────────────────────────────────────────┐
│                 Sampling Table Types                         │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  SimLinAliasData (General Alias Sampling)                   │
│  ┌─────────────────────────────────────────────────────┐    │
│  │ Value[0] │ Value[1] │ ... │ Value[N-1]               │    │
│  ├─────────────────────────────────────────────────────┤    │
│  │ Prob[0]  │ Prob[1]  │ ... │ Prob[N-1]                │    │
│  ├─────────────────────────────────────────────────────┤    │
│  │ Alias[0] │ Alias[1] │ ... │ Alias[N-1]               │    │
│  └─────────────────────────────────────────────────────┘    │
│                                                             │
│  SimGSTables (Goudsmit-Saunderson Angular)                  │
│  ┌─────────────────────────────────────────────────────┐    │
│  │ For each energy E:                                   │    │
│  │   TransformParam 'a'                                │    │
│  │   ┌─────────┬─────────┬────────┐                   │    │
│  │   │ VarU[i] │ ParmA[i]│ParmB[i]│  for i=0..N-1     │    │
│  │   └─────────┴─────────┴────────┘                   │    │
│  └─────────────────────────────────────────────────────┘    │
└─────────────────────────────────────────────────────────────┘
                    │
                    ▼
┌─────────────────────────────────────────────────────────────┐
│              Material Data Organization                     │
├─────────────────────────────────────────────────────────────┤
│  SimMaterialData                                            │
│  ├─ fNumMaterial = N                                        │
│  ├─ fMaterialName[0..N-1]  (e.g., "G4_WATER")             │
│  ├─ fMaterialDensity[0..N-1]  [g/cm³]                      │
│  ├─ fElectronCut, fGammaCut  [MeV]                         │
│  └─ fMollerIMFPScaling[0..N-1]  (Z/A ratio factors)       │
│                                                             │
│  Material Indexing in Geometry:                             │
│  ──► Fast lookup: material = geom.GetMaterialIndex(voxel)  │
│  ──► Density: density = geom.GetVoxelDensity(material)     │
└─────────────────────────────────────────────────────────────┘
```

### 4.2 Physics Data Files

| Process | Data File | Content |
|---------|-----------|---------|
| Compton | `compton_KNDtrData.dat` | Klein-Nishina sampling tables |
| Bremsstrahlung | `brem_SBDtrData.dat` | Seltzer-Berger DCS tables |
| Moller | `ioni_MollerDtrData.dat` | Moller energy transfer tables |
| MSC | `el_GSDtrData.dat` | Goudsmit-Saunderson angular tables |
| IMFPs | Various files | Inverse mean free paths |

### 4.3 Material Data Structure

- **Per-material properties**: Density, Z/A ratios, cuts
- **Scaling factors**: For cross-section interpolation
- **Index-based lookup**: Fast material identification

---

## 5. Transport Algorithm Detailed Flow

### 5.1 Initialization per Primary
1. Create primary track with initial conditions
2. Compute initial interaction distances:
   - numTr1MFP = S_max(E)/tr1-mfp(E) × random
   - numMollerMFP = -ln(random)
   - numBremMFP = -ln(random)

### 5.2 Step Algorithm (per iteration)
1. **Geometric limiting**: distance to voxel boundary
2. **Physics limiting**: distance to each interaction
3. **Select minimum**: actual step length
4. **Energy loss**: compute dE = s × dE/dx
5. **Update counters**: decrease remaining mfp fractions
6. **Apply interaction**: based on what happened

### 5.3 Interaction Handling
- **MSC hinge**: Sample and apply angular deflection
- **Discrete interaction**: Sample energy transfer, create secondaries
- **Boundary crossing**: Update voxel indices, restart competition

#### 5.3.1 Transport Algorithm Flow Diagram

```
           ┌─────────────────┐
           │   Start Track   │
           └─────────┬───────┘
                     │
           ┌─────────▼───────┐
           │ Compute Distance │
           │   to Boundary    │
           └─────────┬───────┘
                     │
        ┌────────────▼─────────────┐
        │ Compute Physics Distances│
        │ • Step_MSC = numTr1MFP/ │
        │   (1/tr1mfp)            │
        │ • Step_Moller = numMol/ │
        │   (1/IMFP)              │
        │ • Step_Brem = numBrem/  │
        │   (1/IMFP)              │
        └────────────┬────────────┘
                     │
        ┌────────────▼─────────────┐
        │ Find Minimum of:         │
        │ • Boundary distance      │
        │ • Step_MSC              │
        │ • Step_Moller           │
        │ • Step_Brem             │
        └────────────┬────────────┘
                     │
         ┌───────────▼───────────┐
         │ Take Step of Length L  │
         │ Update Position       │
         │ Compute Energy Loss    │
         └───────────┬───────────┘
                     │
         ┌───────────▼───────────┐
         │ What Happened?        │
         └─────┬─────────┬───────┘
               │         │
      ┌────────▼─┐ ┌─────▼──────┐
      │Boundary  │ │ Interaction│
      │Crossing  │ │ Type?      │
      └─────┬────┘ └─────┬──────┘
            │            │
   ┌────────▼────┐ ┌─────▼─────┐
   │Update Voxel │ │MSC/Moller/│
   │Indices     │ │Brem?      │
   └─────┬──────┘ └─────┬─────┘
         │              │
         └──────┬───────┘
                │
        ┌───────▼───────┐
        │Decrease Counters│
        │(numTr1, numMol,│
        │numBrem)         │
        └───────┬───────┘
                │
        ┌───────▼───────┐
        │Energy < Cut?  │◄─────────┐
        └───────┬───────┘          │
                │                  │ Yes
          ┌─────▼─────┐            │
          │   No      │            │
          └─────┬─────┘            │
                │                  │
        ┌───────▼───────┐          │
        │Continue Loop  │          │
        └───────┬───────┘          │
                │                  │
                └──────────────────┘
                         │
                         ▼
               ┌─────────────────┐
               │  Kill Track     │
               │ Deposition     │
               └─────────────────┘
```

#### 5.3.2 Photon Transport Flow

```
           ┌─────────────────┐
           │   Start Photon  │
           └─────────┬───────┘
                     │
           ┌─────────▼───────┐
           │ Sample Step:    │
           │ s = -ln(R)/σ_max│
           └─────────┬───────┘
                     │
           ┌─────────▼───────┐
           │ Move Position   │
           │ Update Track    │
           └─────────┬───────┘
                     │
           ┌─────────▼───────┐
           │ Inside Geometry?│◄─────┐
           └─────────┬───────┘      │
                     │              │ No
                ┌────▼────┐         │
                │   Yes   │         │
                └────┬────┘         │
                     │              │
           ┌─────────▼───────┐      │
           │ Real Interaction?│◄────┤
           │ R < 1-σ/σ_max?  │      │
           └─────────┬───────┘      │
                     │              │ Yes
                ┌────▼────┐         │
                │   No    │         │
                └────┬────┘         │
                     │              │
           ┌─────────▼───────┐      │
           │ Continue Loop   │      │
           └─────────────────┘      │
                     │              │
                     └──────────────┘
                         │ Yes
                         ▼
                ┌─────────────────┐
                │Which Interaction│
                │Compton/Photo/   │
                │Pair?           │
                └─────┬───────────┘
                      │
      ┌───────────────┼───────────────┐
      │               │               │
┌─────▼────┐   ┌─────▼────┐   ┌─────▼────────┐
│Compton   │   │Photo     │   │Pair        │
│Scatter   │   │Electric  │   │Production  │
│          │   │Effect    │   │            │
└─────┬────┘   └─────┬────┘   └─────┬───────┘
      │               │               │
Track γ     Deposit all      Create e⁻/e⁺
Create e⁻   energy          Secondaries
```

---

## 6. Performance Optimizations

### 6.1 Algorithmic Optimizations

1. **No conditional loops** (deterministic execution)
2. **Fixed random numbers** per interaction
3. **Fast geometry navigation** (voxel indexing)
4. **Prefetched data structures** (contiguous memory)

### 6.2 Physics Approximations

1. **Moller IMFP independence** from energy
2. **Forward bremsstrahlung** approximation
3. **Uniform pair production** energy sharing
4. **Simplified pair production** kinematics

### 6.3 Memory Access Patterns

- **Structure of Arrays** for better cache utilization
- **Sequential data access** for interpolation
- **Minimal dynamic allocation** during simulation

---

## 7. Medical Physics Considerations

### 7.1 Energy Range
- **Electrons**: 100 keV - 20 MeV
- **Photons**: 50 keV - 20 MeV
- **Cuts**: e⁻ cut = 200 keV, γ cut = 50 keV (typical)

### 7.2 Dosimetric Accuracy
- **Voxel scoring**: Energy deposition per depth bin
- **Track length**: Accurate path computation
- **Secondary particles**: All above threshold tracked

### 7.3 Clinical Relevance
- **Water equivalence** calibration
- **Heterogeneity** handling (bone, titanium)
- **Depth dose** curves validation against reference

---

## 7. Memory Management and Workflow

### 7.1 Memory Architecture

#### 7.1.1 Static vs Dynamic Memory Allocation

The DPM-g4cpp framework is designed with a **predominantly static memory model** to support deterministic execution:

```
Memory Allocation Strategy
────────────────────────────────────────────────────────────

┌─────────────────────────────────────────────────────────────┐
│                    STATIC MEMORY                            │
│  (Allocated during initialization)                          │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  Physics Tables (Loaded from files)                         │
│  ├─ Sampling tables (KN, SB, Moller, GS)                   │
│  ├─ IMFP data (per material)                                │
│  ├─ Stopping power data                                     │
│  └─ Material properties                                     │
│                                                             │
│  Geometry Data                                              │
│  ├─ Voxel arrays (material indices)                         │
│  ├─ Density arrays                                          │
│  └─ Energy deposition histograms                            │
│                                                             │
│  Track Stack (Fixed size pool)                              │
│  └─ Pre-allocated Track objects                            │
│                                                             │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│                 MINIMAL DYNAMIC MEMORY                      │
│  (Only during initialization)                               │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  File I/O buffers (temporary)                              │
│  └─ Data loading from files to static structures           │
│                                                             │
│  String storage (material names, etc.)                     │
│  └─std::string objects                                     │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

#### 7.1.2 Track Stack Implementation

The TrackStack uses a **fixed-size pool allocator** pattern:

```cpp
class TrackStack {
private:
    static const int kMaxTracks = 1000000;  // Configurable
    Track fTrackPool[kMaxTracks];           // Static array
    int fCurrentIndex;                       // Stack pointer
    int fActiveCount;                        // Number of active tracks
public:
    // Returns reference to next available Track
    Track& Insert() {
        // Assertion: fCurrentIndex < kMaxTracks
        return fTrackPool[fCurrentIndex++];
    }

    // Pops track from stack, returns true on success
    bool PopIntoThisTrack(Track& track) {
        if (fCurrentIndex > 0) {
            --fCurrentIndex;
            track = fTrackPool[fCurrentIndex];
            return true;
        }
        return false;
    }
};
```

#### 7.1.3 Data Locality Optimization

The framework prioritizes **cache-friendly data access patterns**:

1. **Structure of Arrays (SoA) Pattern**:
   - Physics data stored per process, not per particle
   - Energy grids contiguous in memory
   - Material properties aligned for SIMD access

2. **Sequential Table Access**:
   - Energy interpolation uses adjacent array entries
   - No pointer chasing or indirect access during simulation

3. **Voxel Data Layout**:
   - 3D geometry stored as 1D arrays with index calculation
   - Material indices: `idx = ix + nx*(iy + ny*iz)`

### 7.2 Execution Workflow

#### 7.2.1 Complete Simulation Timeline

```text
DPM Simulation Execution Flow
────────────────────────────────────────────────────────────

┌─────────────────────────────────────────────────────────────┐
│                    PHASE 0: Setup                          │
│                                                             │
│  1. Parse command line arguments                            │
│  2. Load geometry configuration                             │
│  3. Initialize material data                                │
│  4. Load physics tables from files                          │
│  5. Allocate static memory structures                       │
│  6. Initialize random number generator                     │
│                                                             │
└─────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────┐
│                PHASE 1: Primary Loop                       │
│                                                             │
│  For each primary particle (1..N_primaries):               │
│  ├─ 1. Initialize primary track                            │
│  │   - Position, direction, energy                         │
│  │   - Type (e-/γ), material index                         │
│  └─ 2. Push track to stack                                 │
│                                                             │
│  While stack not empty:                                     │
│  ├─ 1. Pop next track                                     │
│  ├─ 2. Transport until death/exit/cut                     │
│  ├─ 3. Score energy deposition                            │
│  └─ 4. Generate and push secondaries (if any)             │
│                                                             │
└─────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────┐
│                PHASE 2: Finalization                       │
│                                                             │
│  1. Normalize dose histograms                              │
│  2. Write results to file (hist.sim)                      │
│  3. Print statistics                                       │
│  4. Clean up (automatic - static memory)                  │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

#### 7.2.2 Per-Track Transport Workflow

```text
Individual Track Transport
───────────────────────────────────────────────

┌─────────────────────────────────────────────┐
│           Track Transport Loop               │
│                                              │
│  while track.Ekin > cut:                     │
│  │                                           │
│  │  // 1. Geometry Step                      │
│  │  step_boundary = DistToBoundary(pos, dir) │
│  │  material = GetMaterial(voxel)            │
│  │                                           │
│  │  // 2. Physics Competition                 │
│  │  step_physics = min(                       │
│  │    step_msc = numTr1MFP * tr1_mfp,         │
│  │    step_moller = numMolMFP * mol_mfp,     │
│  │    step_brem = numBremMFP * brem_mfp)     │
│  │                                           │
│  │  // 3. Select actual step                 │
│  │  step = min(step_boundary, step_physics)  │
│  │                                           │
│  │  // 4. Transport                         │
│  │  pos += dir * step                        │
│  │  Ekin_loss = step * dE_dx(Ekin, mat)      │
│  │  Ekin -= Ekin_loss                        │
│  │                                           │
│  │  // 5. Update counters                    │
│  │  numTr1MFP -= step / tr1_mfp              │
│  │  numMolMFP -= step / mol_mfp              │
│  │  numBremMFP -= step / brem_mfp            │
│  │                                           │
│  │  // 6. Interaction handling               │
│  │  if step == step_physics:                 │
│  │    - Sample and apply interaction         │
│  │    - Generate secondaries                 │
│  │                                           │
│  └───────────────────────────────────────────┘
```

#### 7.2.3 Secondary Particle Handling

```
Secondary Particle Management
────────────────────────────────────────────────────────────

┌─────────────────────────────────────────────────────────────┐
│                 Generation Flow                             │
│                                                             │
│  Primary Track ──► Interaction Event                        │
│                        │                                    │
│          ┌─────────────┼─────────────┐                      │
│          │             │             │                      │
│    Compton γ      Moller e⁻     Bremsstrahlung γ           │
│          │             │             │                      │
│   Create e⁻      Create e⁻      Create γ                    │
│          │             │             │                      │
│          └─────────────┼─────────────┘                      │
│                        │                                    │
│               Check Energy > Cut?                           │
│                        │                                    │
│                   ┌────┴────┐                               │
│                  Yes│        │No                           │
│                   ┌─▼────┐  ┌─▼─────────┐                    │
│                   │ Push │  │Deposit    │                    │
│                   │to    │  │Energy      │                    │
│                   │Stack │  │Locally     │                    │
│                   └──────┘  └────────────┘                    │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

### 7.3 Performance Characteristics

#### 7.3.1 Memory Footprint

Typical memory usage for a standard configuration:

| Component | Size | Description |
|-----------|------|-------------|
| Physics Tables | ~50-100 MB | All sampling tables for 3 materials |
| Geometry | ~1-10 MB | Voxel data for 1mm voxels, 20cm³ volume |
| Track Stack | ~64 MB | 1M tracks × 64 bytes per track |
| Total | ~115-174 MB | Fits easily in CPU cache |

#### 7.3.2 Cache Optimization Strategies

1. **Prefetching**: Physics tables accessed sequentially
2. **Data Alignment**: All structures aligned to cache line boundaries
3. **Temporal Locality**: Recent particles access same geometry regions
4. **Spatial Locality**: Related data stored contiguously

#### 7.3.3 Thread Safety Considerations

The current implementation is **single-threaded** but designed for easy parallelization:

- TrackStack could be made per-thread
- Physics tables are read-only (safe to share)
- Geometry is read-only during transport
- Only dose scoring needs atomic operations

---

## 8. Extension for Voxel-wise TIA and CT Data Support

### 8.1 Overview

To expand DPM-g4cpp for internal dosimetry applications, we propose adding support for:
- **Voxel-wise Time-Integrated Activity (TIA)** maps in Bq·s units
- **CT Hounsfield Unit (HU)** maps for tissue density and composition
- **Radionuclide configuration** for decay schemes and emission spectra
- **NIfTI format support** for medical image compatibility

### 8.2 Data Structures for Voxel-wise Support

#### 8.2.1 TIA Data Structure

```cpp
class SimTIAData {
private:
    // 3D TIA distribution (Bq·s)
    std::vector<double> fTIAMap;

    // Voxel dimensions
    int fNx, fNy, fNz;
    double fVoxelSize[3];  // mm

    // Radionuclide information
    struct Radionuclide {
        std::string name;           // e.g., "Lu-177", "I-131"
        double halfLife;           // seconds
        std::vector<Emission> emissions;
    } fRadionuclide;

    struct Emission {
        int type;                   // 0=beta-, 1=beta+, 2=gamma, 3=Auger
        double energy;              // MeV
        double intensity;           // decay fraction
        double endpointEnergy;      // for beta particles
    };

public:
    void LoadNIFTI(const std::string& filename);
    double GetTIA(int ix, int iy, int iz) const;
    const Radionuclide& GetRadionuclide() const { return fRadionuclide; }
    // Convert TIA to source term: decays per second
    double GetSourceStrength(int ix, int iy, int iz) const;
};
```

#### 8.2.2 CT Data Structure

```cpp
class SimCTData {
private:
    // CT Hounsfield units
    std::vector<short> fCTMap;

    // Derived material composition per HU range
    struct MaterialMapping {
        double HU_min, HU_max;
        int materialIndex;          // Index in material database
        double density;             // g/cm³
        double effectiveZ;          // Effective atomic number
    };
    std::vector<MaterialMapping> fHU_to_Material;

    // Voxel dimensions (must match TIA)
    int fNx, fNy, fNz;
    double fVoxelSize[3];

public:
    void LoadNIFTI(const std::string& filename);
    void SetupHUToMaterialMapping();
    int GetMaterialIndex(int ix, int iy, int iz) const;
    double GetDensity(int ix, int iy, int iz) const;
    void ConvertToMaterialData(SimMaterialData& matData) const;
};
```

#### 8.2.3 Integrated Geometry for Medical Data

```cpp
class MedicalGeometry : public Geom {
private:
    SimTIAData* fTIAData;
    SimCTData* fCTData;

    // Emission sampling cache
    mutable std::vector<double> fCumulativeEmission;

public:
    MedicalGeometry(SimCTData* ctData, SimTIAData* tiaData);

    // Source position sampling from TIA distribution
    SampleResult SampleSourcePosition() const;

    // Material index with CT-based density
    int GetMaterialIndex(int* iVoxel) const override;
    double GetVoxelMaterialDensity(int* iVoxel) const override;

    // Voxel source strength for decay sampling
    double GetVoxelSourceStrength(int* iVoxel) const;

    // Sample emission type and energy
    std::pair<int, double> SampleEmission() const;
};
```

### 8.3 NIfTI File Integration

#### 8.3.1 NIfTI Reader Implementation

```cpp
#include <nifti1_io.h>

class NIfTIReader {
public:
    struct Header {
        int dim[8];                // Image dimensions
        float pixdim[8];           // Voxel sizes (mm)
        float scl_slope;            // Data scaling
        float scl_inter;
        std::string description;   // Optional description
    };

    template<typename T>
    static std::vector<T> LoadNIFTI(const std::string& filename, Header& header) {
        nifti_image* nim = nifti_image_read(filename.c_str(), 1);
        if (!nim) throw std::runtime_error("Failed to read NIfTI");

        // Extract header info
        for (int i = 0; i < 8; i++) {
            header.dim[i] = nim->dim[i];
            header.pixdim[i] = nim->pixdim[i];
        }
        header.scl_slope = nim->scl_slope;
        header.scl_inter = nim->scl_inter;

        // Extract data
        std::vector<T> data(nim->nvox);
        T* ptr = static_cast<T*>(nim->data);
        std::copy(ptr, ptr + nim->nvox, data.begin());

        nifti_image_free(nim);
        return data;
    }
};
```

### 8.4 Radionuclide Configuration

#### 8.4.1 Configuration File Format (JSON)

```json
{
    "radionuclides": {
        "Lu-177": {
            "half_life": 6.647e6,
            "emissions": [
                {
                    "type": "beta-",
                    "endpoint_energy": 0.497,
                    "mean_energy": 0.133,
                    "intensity": 0.786
                },
                {
                    "type": "gamma",
                    "energy": 0.208,
                    "intensity": 0.11
                },
                {
                    "type": "gamma",
                    "energy": 0.113,
                    "intensity": 0.064
                }
            ]
        },
        "I-131": {
            "half_life": 6.93e6,
            "emissions": [
                {
                    "type": "beta-",
                    "endpoint_energy": 0.606,
                    "mean_energy": 0.182,
                    "intensity": 0.894
                },
                {
                    "type": "gamma",
                    "energy": 0.364,
                    "intensity": 0.812
                }
            ]
        }
    },
    "hu_to_material_mapping": [
        {
            "hu_min": -1000,
            "hu_max": -500,
            "material": "G4_AIR",
            "density": 0.001205
        },
        {
            "hu_min": -100,
            "hu_max": 50,
            "material": "G4_WATER",
            "density": 1.0
        },
        {
            "hu_min": 100,
            "hu_max": 300,
            "material": "G4_SOFT_TISSUE",
            "density": 1.05
        },
        {
            "hu_min": 700,
            "hu_max": 3000,
            "material": "G4_BONE_COMPACT_ICRU",
            "density": 1.85
        }
    ]
}
```

### 8.5 Internal Dosimetry Workflow

#### 8.5.1 Simulation Flow for Internal Sources

```text
Internal Dosimetry Simulation Flow
────────────────────────────────────────────────────────────

┌─────────────────────────────────────────────────────────────┐
│                    PHASE 0: Data Loading                    │
│                                                             │
│  1. Load TIA NIfTI file (Bq·s)                             │
│  2. Load CT NIfTI file (HU)                                │
│  3. Load radionuclide configuration                        │
│  4. Build tissue composition map from CT                   │
│  5. Generate physics tables for all tissues                │
│                                                             │
└─────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────┐
│                PHASE 1: Source Sampling                     │
│                                                             │
│  For each primary decay event:                             │
│  ├─ Sample voxel from TIA distribution                     │
│  │   - Normalize TIA to cumulative distribution            │
│  │   - Rejection sampling on 3D grid                       │
│  ├─ Sample emission type and energy                        │
│  │   - Based on radionuclide decay scheme                  │
│  ├─ Set initial particle state                             │
│  │   - Position: voxel center                              │
│  │   - Direction: isotropic (usually)                      │
│  │   - Energy: emission energy                             │
│  └─ Push particle to transport stack                       │
│                                                             │
└─────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────┐
│                PHASE 2: Transport                           │
│                                                             │
│  Standard DPM transport with modifications:                 │
│  ├─ Material lookup from CT-derived map                    │
│  ├─ Density variations per voxel                           │
│  ├─ Score absorbed dose per voxel                          │
│  └─ Track emission probability                             │
│                                                             │
└─────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────┐
│                PHASE 3: Dose Calculation                    │
│                                                             │
│  1. Convert deposited energy to dose:                      │
│     dose = energy_deposited / (mass * TIA)                 │
│  2. Create 3D dose map (Gy per decay)                      │
│  3. Save results in NIfTI format                           │
│  4. Calculate total dose: dose_total = dose * administered  │
│     activity (MBq) × time (s)                              │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

### 8.6 Implementation Details

#### 8.6.1 Source Sampling Algorithm

```cpp
class VoxelSourceSampler {
private:
    const SimTIAData* fTIAData;
    std::vector<double> fCumulativeTIA;
    double fTotalTIA;

public:
    VoxelSourceSampler(const SimTIAData* tiaData) : fTIAData(tiaData) {
        BuildCumulativeDistribution();
    }

    void BuildCumulativeDistribution() {
        int n = fTIAData->GetNumberOfVoxels();
        fCumulativeTIA.resize(n + 1, 0.0);
        fTotalTIA = 0.0;

        for (int i = 0; i < n; i++) {
            fTotalTIA += fTIAData->GetTIA(i);
            fCumulativeTIA[i + 1] = fTotalTIA;
        }

        // Normalize
        for (double& val : fCumulativeTIA) {
            val /= fTotalTIA;
        }
    }

    int SampleVoxel(double r) const {
        // Binary search in cumulative distribution
        auto it = std::upper_bound(fCumulativeTIA.begin(),
                                 fCumulativeTIA.end(), r);
        return std::distance(fCumulativeTIA.begin(), it) - 1;
    }
};
```

#### 8.6.2 Dose Scoring Enhancement

```cpp
class VoxelDoseScorer {
private:
    // 3D dose accumulation
    std::vector<double> fDoseMap;
    std::vector<double> fMassMap;
    const MedicalGeometry* fGeometry;

public:
    void ScoreEnergy(double energyDeposit, int ix, int iy, int iz) {
        int idx = fGeometry->GetVoxelIndex(ix, iy, iz);
        fDoseMap[idx] += energyDeposit;  // MeV

        // Convert to Gy per decay later
        // Gy = MeV × 1.602e-13 J/MeV / (mass in kg)
    }

    void NormalizeByTIAAndSave(const SimTIAData& tiaData,
                               const std::string& filename) {
        for (int i = 0; i < fDoseMap.size(); i++) {
            double tia = tiaData.GetTIA(i);
            double mass = fMassMap[i];  // kg

            // Gy per decay
            if (tia > 0 && mass > 0) {
                fDoseMap[i] *= 1.602e-13 / (mass * tia);
            }
        }

        SaveAsNIfTI(filename);
    }
};
```

### 8.7 Usage Example

```bash
# Generate data for tissue materials
./dpm_GenerateData -c 4 -d data_tissues -t config/tissues.json

# Run internal dosimetry simulation
./dpm_Simulate_Internal \
  --tia-data patient_Lu177_tia.nii \
  --ct-data patient_ct.nii \
  --radionuclide Lu-177 \
  --tissue-config config/tissues.json \
  --output dose_map.nii
```

### 8.8 Validation Strategy

#### 8.8.1 Test Cases

1. **Uniform TIA in water phantom** - Compare to analytical solution
2. **Point source in heterogeneous phantom** - Verify tissue interfaces
3. **Clinical SPECT/CT data** - Compare to commercial TPS
4. **Multi-radionuclide scenarios** - Verify mixed dose calculations

#### 8.8.2 Output Formats

```nifti
# Dose map output (NIfTI header info)
datatype = 64 (FLOAT64)          // 64-bit float for precision
pixdim = [dx, dy, dz, dt, ...]  // Voxel sizes
scl_slope = 1.0                  // Gy per unit
scl_inter = 0.0                  // No offset
descrip = "DPM absorbed dose Gy per decay"
```

---

## 9. Conclusions

The DPM-g4cpp framework represents a sophisticated balance between physical accuracy and computational efficiency:

### 9.1 Strengths
- **Deterministic execution** ideal for parallel hardware
- **Rejection-free sampling** eliminates thread divergence
- **Optimized for radiotherapy** energy ranges
- **Well-validated** against original DPM

### 9.2 Approximations
- **Energy-independent Moller** cross section
- **Simplified angular distributions** for some processes
- **Fixed reference material** scaling

### 9.3 Innovation
- **Two-step MSC** eliminates continuous scattering
- **Pre-computed tables** for all physics processes
- **Clean architecture** facilitating GPU/FPGA ports

This implementation successfully modernizes the original DPM algorithm while maintaining its performance advantages for treatment planning applications.

### 9.4 Extension Potential
The proposed voxel-wise support for TIA and CT data would significantly expand DPM-g4cpp's capabilities:
- **Internal dosimetry** for nuclear medicine applications
- **Patient-specific dose calculations** using actual clinical data
- **Multi-radionuclide therapy** planning and optimization
- **Integration with existing PACS** and treatment planning systems

The modular architecture ensures these extensions can be implemented without compromising the core performance optimizations that make DPM-g4cpp uniquely suited for real-time dose calculations.