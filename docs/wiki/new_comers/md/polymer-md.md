# Kremer–Grest Bead–Spring Polymer Model

>A minimal yet powerful coarse‑grained polymer model

!!! note "Compatibility"
    This tutorial is compatible with the 22Jul2025 LAMMPS version.

!!! info "Cite"
    This tutorial is originally developed by Dr. Zhenghao Wu.

## Introduction

The Kremer–Grest (KG) model is a widely used coarse-grained polymer model in which beads are connected by nonlinear springs. Non-bonded interactions are described by a standard Lennard-Jones potential, while chain connectivity is maintained using finitely extensible nonlinear elastic (FENE) bonds. This simple yet powerful model reproduces key features of polymer melts—including generic melt behavior, reptation, entanglement effects (for sufficiently long chains), and Rouse dynamics—without requiring atomistic detail.

In reduced Lennard–Jones units, the canonical KG parameters are:

- Non‑bonded: LJ with cutoff \(r_c = 2.5 \sigma, \epsilon = 1\), \(\sigma = 1\)
- Bonds: FENE with spring constant \(k = 30\,\epsilon/\sigma^2\) and maximum extension \(R_0 = 1.5\,\sigma\)
- Thermostat: Langevin or Nosé–Hoover targeting \(T = 1.0\)
- Time step: \(\Delta t = 0.005\)

> **Note**: [Here](https://hoomd-blue.readthedocs.io/en/stable/units.html) to learn more about unit conversion 

## My first polymer model

We will start from scratch from building single polymer chain, many polymer chains, and MD simulations of the system using KG model in LAMMPS. We will build a short melt of flexible chains using the KG model. The script is split into five categories: 1. Model construction; 2. Initilization; 3. System definition; 4. Settings; 5. Short equilibration

Create a new folder, create a new file called `kg_initial.lmp`, and copy the following lines into it.

```lammps
# PART A - KG melt setup & short equilibration
# 1) Initialization
# 2) System definition
# 3) Settings
# 4) Monitoring
# 5) Run
```

### Model construction

We start from constructing single short flexible chain model with 10 beads by creating a xyz file:

!!! hints "single_chain.xyz"
    ```
    10
    # type x y z
    1 0.0 0.0 0.0
    1 0.97 0.0 0.0
    1 1.94 0.0 0.0
    1 2.91 0.0 0.0
    1 3.88 0.0 0.0
    1 4.85 0.0 0.0
    1 5.82 0.0 0.0
    1 6.79 0.0 0.0
    1 7.76 0.0 0.0
    1 8.73 0.0 0.0
    ```

After creating single_chain.xyz, we need to edit the packmol input file to assemble the single chain to multiple chains:

!!! hints "many_chain_packmol.inp"
    ```
    tolerance 2.0
    filetype xyz
    output many_chain.xyz

    structure single_chain.xyz 
      number 100 
      inside box 0. 0. 0. 25. 25. 25. 
    end structure
    ```
Execute Packmol with the command: */path/packmol < many_chain_packmol.inp*. This will generate the file many_chain.xyz, a snapshot of which is shown below:

<div style="justify-content: center; width: 400px;">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="../../../../images/lammps/tutorial2/figures/many_chain.png">
    <img alt="Short bead–spring chains in a melt" src="../../../../images/lammps/tutorial2/figures/many_chain.png">
  </picture>
</div>

_Example snapshots of short bead–spring chains assembled into many chain by PACKMOL._

Now we have obtained the coordinates for the system, next we need to add bond topology for the 100 polymer chains, each with 9 bonds.

we can simply write a python code to generate the bond topology following the lammps data format:

!!! hints "generate_bond.py"
    ```
    all_bonds = []

    def generate_bond(number_chain,chain_length):
        bond_id=1
        for i_chain in range(number_chain):
            for i_bead in range(chain_length-1):
                # same bond type for this particular system, default 1
                _bond = [bond_id, 1, i_bead+i_chain*chain_length+1,i_bead+1+i_chain*chain_length+1]
                bond_id += 1
                all_bonds.append(_bond)
        return all_bonds
        
    all_bonds = generate_bond(100,10)
    with open('bonds.txt','w') as write_f:
        for _bond in all_bonds:
            write_f.write(f'{_bond[0]} {_bond[1]} {_bond[2]} {_bond[3]}\n')

    ```

After running the python script, you will get "bonds.txt". Assemble all these information into the lammps data file for this system:

!!! hints "bead_spring_polymer.data"
    ```
    LAMMPS Description           (1st line of file)

    1000 atoms         (this must be the 3rd line, 1st 2 lines are ignored)
    900 bonds                (# of bonds to be simulated)

    1 atom types           (# of nonbond atom types)
    1 bond types          (# of bond types = sets of bond coefficients)

    0 26.0 xlo xhi       (for periodic systems this is box size,
    0 26.0 ylo yhi        for non-periodic it is min/max extent of atoms)
    0 26.0 zlo zhi       (do not include this line for 2-d simulations)

    Masses

      1 1.0 
                          

    Atoms

    {atom information}


    Bonds

    {bond information}
    ```

Here the atom information is the definition of particles in the system, it corresponds to the *atom_style* command in lammps input script. We use *atom_style=bond* in this example, thus we follow the requirement to have atom information in this format

```text
# index molecule-tag atom-type x y z
1 molecule-tag atom-type x_atom1 y_atom1 y_atom1
2 molecule-tag atom-type x_atom2 y_atom2 y_atom2
...
```

The bond format has been implemented in the python script, thus we only need to copy paste the data in bonds.txt to the *bond information* session.

>Note: [Here](https://docs.lammps.org/2001/data_format.html) to learn more details about the format of LAMMPS data file

### Initialization

Set LJ units and choose an atom style that supports bonds. Update the first section to:

```lammps
# 1) Initialization
units           lj
dimension       3
atom_style      bond
boundary        p p p
neighbor        0.4 bin
neigh_modify    every 1 delay 0 check yes
```

### System definition

Different from the example of simple LJ liquid, we have defined the system in the lammps data file. Here we need to load the data file into the simulation. Add the following to the “System definition” category:

```lammps
# 2) System definition
read_data bead_spring_polymer.data
```

### Settings

Specify the KG non‑bonded and bonded interactions. Add the following to the “Settings” category:

```lammps
# 3) Settings

# LJ with rc = 2.5
pair_style      lj/cut 2.5
pair_modify     shift yes
pair_coeff      1 1 1.0 1.0

# FENE bonds: k = 30, R0 = 1.5, epsilon = 1.0, sigma = 1.0
bond_style      fene
bond_coeff      1 30.0 1.5 1.0 1.0

# Recommended for KG: exclude 1-2 LJ, include 1-3 and 1-4
special_bonds   fene
```

!!! note
    `special_bonds fene` sets the 1–2 LJ scaling to 0 and 1–3, 1–4 to 1, which is the standard choice for KG chains with WCA non‑bonded interactions.

### Short equilibration

Define thermo output, a simple image dump, and perform a short NVT run with a Langevin thermostat:

```lammps
# 4) Monitoring
thermo          1000
thermo_style    custom step temp pe ke etotal press
dump            viz all image 5000 kg-*.ppm type type size 800 800 zoom 1.6 fsaa yes
dump_modify     viz adiam 1 1 acolor 1 royalblue backcolor white box no

# 5) Run
velocity        all create 1.0 48291 mom yes dist gaussian
fix             mynve all nve
fix             mylgv all langevin 1.0 1.0 1.0 10917 zero yes
timestep        0.005
run             50000
```

You should observe fast relaxation toward a homogeneous melt structure at the target temperature.

!!! note
    The KG model is usually run without a prior energy minimization. The Langevin thermostat with modest damping (here 1.0) stabilizes the initial dynamics.

## Improving the script

We next save a minimized, relaxed configuration to restart longer simulations without re‑creating molecules.

### Save a melt and restart from data

Create `kg_save.lmp` that writes a data file after an initial thermalization:

```lammps
# 1) Initialization
units           lj
dimension       3
atom_style      bond
boundary        p p p

# 2) System definition
read_data bead_spring_polymer.data

# 3) Settings
pair_style      lj/cut 2.5
pair_modify     shift yes
pair_coeff      1 1 1.0 1.0
bond_style      fene
bond_coeff      1 30.0 1.5 1.0 1.0
special_bonds   fene

# 4) Monitoring
thermo          1000
thermo_style    custom step temp pe ke etotal press

# 5) Run
velocity        all create 1.0 98765 mom yes dist gaussian
fix             mynve all nve
fix             mylgv all langevin 1.0 1.0 1.0 7777 zero yes
timestep        0.005
run             50000

# 6) Save system
write_data      kg_melt.data
```

<div style="justify-content: center; width: 400px;">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="../../../../images/lammps/tutorial2/figures/many_chain_relax.png">
    <img alt="Short bead–spring chains in a melt" src="../../../../images/lammps/tutorial2/figures/many_chain_relax.png">
  </picture>
</div>

_Example snapshots of bead–spring polymer melts after quick relaxtaion._

Now prepare `kg_restart.lmp` that restarts from the saved data file and continues at constant temperature and pressure if desired.

```lammps
# 1) Initialization
units           lj
dimension       3
atom_style      bond
boundary        p p p

# 2) System definition
pair_style      lj/cut 2.5
pair_modify     shift yes
read_data       kg_melt.data

# 3) Settings
pair_coeff      * * 1.0 1.0
bond_style      fene
bond_coeff      1 30.0 1.5 1.0 1.0
special_bonds   fene

# 4) Monitoring
thermo          2000
thermo_style    custom step time temp pe ke etotal press density c_rg c_msd[4]
compute         rg molecule gyration
compute         msd all msd

# 5) Run (choose one)
fix             mynvt all nvt temp 1.0 1.0 1.0
timestep        0.005
run             200000
```

<div style="display: flex; justify-content: center; gap: 24px;">
  <div style="width: 400px;">
    <picture>
      <source media="(prefers-color-scheme: dark)" srcset="../../../../images/lammps/tutorial2/figures/many_chain_equ_unwrap.png">
      <img alt="Short bead–spring chains in a melt (unwrapped)" src="../../../../images/lammps/tutorial2/figures/many_chain_equ_unwrap.png">
    </picture>
  </div>
  <div style="display: flex; flex-direction: column; align-items: center; width: 300px;">
    <div style="height: 180px;"></div>
    <picture>
      <source media="(prefers-color-scheme: dark)" srcset="../../../../images/lammps/tutorial2/figures/many_chain_equ_wrap.png">
      <img alt="Short bead–spring chains in a melt (wrapped)" src="../../../../images/lammps/tutorial2/figures/many_chain_equ_wrap.png">
    </picture>
  </div>
</div>

_Snapshot of the pre-equilibrated polymer melt system. Left: unwrapped representation; Right: wrapped representation_

!!! note
    `compute gyration` gives the instantaneous radius of gyration `c_rg` of the whole group. `compute msd` reports mean‑squared displacement per atom; its fourth column `c_msd[4]` is the total MSD.

### Optional: density and pressure control

To relax the density, replace the thermostat with a barostat for a short time:

```lammps
unfix           mynvt
fix             mynpt all npt temp 1.0 1.0 1.0 iso 0.0 0.0 1.0 
run             50000
```

!!! note
    The KG melt has a typical reduced density around 0.85–0.9 for temperature near \(T=1\). Use a short NPT stage only to remove vacuum or strong pressure drifts

## Going further with exercises

### Chain stiffness (semiflexible KG)

- Add bending stiffness using `angle_style cosine/squared` or `angle_style harmonic`, define one angle type, and add angles to the molecule template (or data). Compare the scaling of the radius of gyration with and without stiffness.

### Pair structure and dynamical observables

- Compute a radial distribution function with `compute rdf` and average it over time with `fix ave/time`.
- Output the shear stress autocorrelation function with `compute pressure` and `fix ave/correlate` to estimate a Green–Kubo viscosity (challenge).

### Avoid “FENE bond blown” errors

- If you observe bond stretching beyond `R0` and runs abort, reduce the time step (e.g., to `0.003`), increase thermostat damping slightly, or equilibrate at lower temperature first. Ensure `special_bonds fene` is set and the non‑bonded cutoff stays at 2.5.

>Reference: Original KG model: K. Kremer and G. S. Grest, J. Chem. Phys. 92, 5057 (1990). See also the LAMMPS documentation for `bond_style fene` and `special_bonds`.

!!! info "Access the input files"
    You can access the input scripts and data files used in these tutorials (and full exercise solutions) from the dedicated [GitHub repository](https://github.com/lammpstutorials/lammpstutorials-inputs/).
