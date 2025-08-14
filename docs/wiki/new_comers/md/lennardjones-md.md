# Lennard-Jones fluid

>The very basics of LAMMPS through a simple example

!!! note "Compatibility"
    This tutorial is compatible with the 22Jul2025 LAMMPS version.

!!! info "Cite"
    If you find these tutorials useful, you can cite “A Set of Tutorials for the LAMMPS Simulation Package” by Simon Gravelle, Jacob R. Gissinger, and Axel Kohlmeyer (2025). See the full paper on [arXiv](https://doi.org/10.48550/arXiv.2503.14020).

## Introduction

The objective of this tutorial is to perform simple MD simulations using LAMMPS. The system consists of a Lennard-Jones fluid composed of neutral particles with two different effective diameters, contained within a cubic box with periodic boundary conditions. In this tutorial, basic MD simulations in the microcanonical (NVE) and canonical (NVT) ensembles are performed, and basic quantities are calculated, including the potential and kinetic energies.

!!! info "Access the input files"
    You can access the input scripts and data files used in these tutorials (and full exercise solutions) from the dedicated [GitHub repository](https://github.com/lammpstutorials/lammpstutorials-inputs/).

<div style="float: right; margin: 0 0 1rem 1rem; width: 300px;">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="../../../../images/lammps/tutorial1/figures/binary_LJ_fluid_dark.webp">
    <img alt="Binary Lennard-Jones mixture used in Tutorial 1" src="../../../../images/lammps/tutorial1/figures/binary_LJ_fluid_dark.webp" style="width: 100%; height: auto;">
  </picture>
</div>

## My first input

To run a simulation using LAMMPS, you need to write an input script containing a series of commands for LAMMPS to execute, similar to Python or Bash scripts. For clarity, the input scripts for this tutorial will be divided into five categories, which will be filled out step by step.

To set up this tutorial using LAMMPS graphical user interface (LAMMPS–GUI), select “Start LAMMPS Tutorial 1” from the “Tutorials” menu and follow the instructions. This will select (or create, if needed) a folder, place the initial input file `initial.lmp` in it, and open the file in the LAMMPS–GUI Editor window:

```lammps
# PART A - ENERGY MINIMIZATION
# 1) Initialization
# 2) System definition
# 3) Settings
# 4) Monitoring
# 5) Run
```

!!! info "If you are not using LAMMPS–GUI"
    Create a new folder and add a file named `initial.lmp` inside it. Open the file in a text editor of your choice and copy the previous lines into it.

Everything that appears after a hash symbol (#) is a comment and ignored by LAMMPS. These five categories are not required in every input script and do not necessarily need to be in that exact order. For instance, the “Settings” and the “Monitoring” categories could be inverted, or the “Monitoring” category could be omitted. However, LAMMPS reads input files from top to bottom and processes each command immediately. Therefore, the “Initialization” and “System definition” categories must appear at the top of the input, and the “Run” category must appear at the bottom. Also, the specifics of some commands can change after global settings are modified, so the order of commands in the input script is important.

### Initialization

In the first section, global parameters for the simulation are defined, such as units, boundary conditions (e.g., periodic or non‑periodic), and atom types (e.g., uncharged point particles or extended spheres with a radius and angular velocities). These commands must be executed before creating the simulation box or they will cause an error. Similarly, many LAMMPS commands may only be entered after the simulation box is defined. Only a limited number of commands may be used in both cases. Update `initial.lmp` so that the “Initialization” section appears as follows:

```lammps
# 1) Initialization
units lj
dimension 3
atom_style atomic
boundary p p p
```

!!! note
    Strictly speaking, none of the four commands specified in the “Initialization” section are mandatory, as they correspond to the default settings for their respective global properties. However, explicitly specifying these defaults is considered good practice to avoid confusion when sharing input files with other LAMMPS users.

The first line, `units lj`, specifies the use of reduced units, where all quantities are dimensionless. This unit system is a popular choice for simulations that explore general statistical mechanical principles, as it emphasizes relative differences between parameters rather than representing any specific material. The second line, `dimension 3`, specifies that the simulation is conducted in 3D space, as opposed to 2D, where atoms are confined to move only in the xy‑plane. The third line, `atom_style atomic`, designates the atomic style for representing simple, individual point particles. In this style, each particle is treated as a point with a mass, making it the most basic atom style. Other atom styles can incorporate additional attributes for atoms, such as charges, bonds, or molecule IDs, depending on the requirements of the simulated model. The last line, `boundary p p p`, indicates that periodic boundary conditions are applied along all three directions.

!!! note
    Each LAMMPS command has extensive online documentation that lists and discusses the different options for that command. Most LAMMPS commands also have default settings that are applied if no value is explicitly specified. From the LAMMPS–GUI editor, you can right‑click on a line containing a command (e.g., `units lj`) and select “View Documentation for `units`” to open the corresponding manual page.

### System definition

Next, create the simulation box and populate it with atoms. Modify the “System definition” category of `initial.lmp` as shown below:

```lammps
# 2) System definition
region simbox block -20 20 -20 20 -20 20
create_box 2 simbox
create_atoms 1 random 1500 34134 simbox overlap 0.3
create_atoms 2 random 100 12756 simbox overlap 0.3
```

The first line defines a region named `simbox` that is a block (rectangular cuboid) extending from −20 to 20 units along all three dimensions. The second line initializes a simulation box based on `simbox` and reserves space for two atom types. In LAMMPS, every atom is assigned an atom type. This property selects which force field parameters (here, the Lennard‑Jones parameters, \(\epsilon_{ij}\) and \(\sigma_{ij}\)) are applied to each pair of atoms via the `pair_coeff` command. The next two lines create atoms of type 1 and type 2 in the box at random positions, enforcing a minimum distance of 0.3 between randomly placed atoms to avoid close contacts.

!!! note
    From this point on, the number of atom types is locked in, and any command referencing an atom type larger than 2 will trigger an error. While it is possible to allocate more atom types than needed, you must assign a mass and provide force field parameters for each atom type.

Another way to define a system is to import an existing topology file containing atomic coordinates (and optionally other attributes) using `read_data` (see later tutorials).

### Settings

Specify the settings for the two atom types. Modify the “Settings” category of `initial.lmp` as follows:

```lammps
# 3) Settings
mass 1 1.0
mass 2 5.0
pair_style lj/cut 4.0
pair_coeff 1 1 1.0 1.0
pair_coeff 2 2 0.5 3.0
```

The two `mass` commands assign masses to types 1 and 2. The `pair_style lj/cut 4.0` specifies that atoms interact via a Lennard‑Jones (LJ) potential with cutoff \(r_c=4.0\):

\[
E_{ij}(r) = 4 \, \epsilon_{ij} \left[ \left( \dfrac{\sigma_{ij}}{r} \right)^{12} - \left( \dfrac{\sigma_{ij}}{r} \right)^{6} \right], \quad r < r_c.
\]

`pair_coeff 1 1 1.0 1.0` sets \(\epsilon_{11}=1.0\), \(\sigma_{11}=1.0\). `pair_coeff 2 2 0.5 3.0` sets \(\epsilon_{22}=0.5\), \(\sigma_{22}=3.0\).

!!! note
    By default, mixed coefficients use geometric mixing: \(\epsilon_{12} = \sqrt{\epsilon_{11}\epsilon_{22}}\) and \(\sigma_{12} = \sqrt{\sigma_{11}\sigma_{22}}\). Here, \(\epsilon_{12} \approx 0.707\), \(\sigma_{12} \approx 1.732\). Other rules can be selected with `pair_modify`.

### Single‑point energy

Complete “Monitoring” and “Run” in `initial.lmp`:

```lammps
# 4) Monitoring
thermo 10
thermo_style custom step etotal press
# 5) Run
run 0 post no
```

`thermo 10` prints thermo data every 10 steps. `thermo_style custom` selects outputs. `run 0 post no` initializes forces/energy without running.

!!! note
    Thermo outputs are instantaneous values, not cumulative averages. Many custom data from computes, fixes, and variables can be added for on‑the‑fly analysis.

You can now run LAMMPS. With LAMMPS–GUI, the Output window shows the total energy and pressure at step 0. Since no simulation steps were performed, the Charts window is empty. You can use the GUI’s Image Viewer to create a snapshot with `dump image`.

### Energy minimization

Replace `run 0 post no` with a minimization:

```lammps
# 5) Run
minimize 1.0e-6 1.0e-6 1000 10000
```

LAMMPS performs an iterative energy minimization (by default conjugate gradient). The run stops when one of the stopping criteria is met. Typically, the potential energy decreases from a positive value to a negative value as atoms move away from overlapping positions. Create and save a snapshot after minimization and compare to the initial image.

### Molecular dynamics

Append MD commands to continue from the minimized state:

```lammps
# PART B - MOLECULAR DYNAMICS
# 4) Monitoring
thermo 50
thermo_style custom step temp etotal pe ke press
```

Then add a second “Run” section for NVE:

```lammps
# 5) Run
fix mynve all nve
timestep 0.005
run 50000
```

This performs MD in the microcanonical (NVE) ensemble with velocity‑Verlet integration.

Now switch to a Langevin thermostat (NVT):

```lammps
# 5) Run
fix mynve all nve
fix mylgv all langevin 1.0 1.0 0.1 10917
timestep 0.005
run 15000
```

A Langevin thermostat drives the temperature to the target with damping 0.1 and random noise (seed 10917). In the presence of a thermostat, the simulation is in the canonical (NVT) ensemble.

 Across these runs, you should observe equilibration: potential energy drops during minimization, then stabilizes; kinetic energy grows during MD to a plateau near the target temperature.

<figure>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="../../../../images/lammps/tutorial1/figures/LJ-energy-dm.png">
    <img alt="Evolution of the Lennard-Jones fluid energy during minimization and MD" src="../../../../images/lammps/tutorial1/figures/LJ-energy.png">
  </picture>
  <figcaption>
    (a) Potential energy U versus step during energy minimization; (b) potential energy U versus time t during MD (NVT);
    (c) kinetic energy K during minimization; (d) kinetic energy K during MD.
  </figcaption>
  </figure>

!!! note
    All simulations here are deliberately short to run on a personal computer. They are not meant to be statistically converged.

### Trajectory visualization

To visualize trajectories, write atom positions periodically:

```lammps
dump mydmp all atom 100 dump.lammpstrj
```

Open `dump.lammpstrj` in VMD or OVITO. To render images during the run, use `dump image`. For example:

```lammps
dump viz all image 100 myimage-*.ppm type type size 800 800 zoom 1.452 shiny 0.7 fsaa yes &
    view 80 10 box yes 0.025 axes no 0.0 0.0 center s 0.483725 0.510373 0.510373
dump_modify viz pad 9 boxcolor royalblue backcolor white adiam 1 1.6 adiam 2 4.8
```

The ampersand `&` continues the command on a new line. These commands create NetPBM images every 100 steps; the two `type` keywords specify color and diameter.

## Improving the script

Let us improve the input and perform more advanced operations, such as specifying initial positions and restarting from a saved configuration.

### Control the initial atom positions

Open `improved.min.lmp`, which contains “Part A” of `initial.lmp` but without any commands in “System definition”:

```lammps
# 1) Initialization
units lj
dimension 3
atom_style atomic
boundary p p p
# 2) System definition
# 3) Settings
mass 1 1.0
mass 2 10.0
pair_style lj/cut 4.0
pair_coeff 1 1 1.0 1.0
pair_coeff 2 2 0.5 3.0
# 4) Monitoring
thermo 10
thermo_style custom step etotal press
# 5) Run
minimize 1.0e-6 1.0e-6 1000 10000
```

Create atoms of types 1 and 2 in two separate regions:

```lammps
# 2) System definition
region simbox block -20 20 -20 20 -20 20
create_box 2 simbox
# for creating atoms
region cyl_in cylinder z 0 0 10 INF INF side in
region cyl_out cylinder z 0 0 10 INF INF side out
create_atoms 1 random 1000 34134 cyl_out
create_atoms 2 random 150 12756 cyl_in
```

The `side in` and `side out` keywords define the inside and outside of the cylinder of radius 10.

Append a sixth section to save the system after minimization:

```lammps
# 6) Save system
write_data improved.min.data
```

!!! note
    `write_data` saves the system state to a text file. We will use it to restart without repeating system creation and minimization.

### Restarting from a saved configuration

Open `improved.md.lmp`, which contains the Initialization part only. Since we read most information from the data file, we don’t repeat the system definition and settings, except for `pair_style`, which must precede `read_data`:

```lammps
# 1) Initialization
units lj
dimension 3
atom_style atomic
boundary p p p
# 2) System definition
# 3) Settings
# 4) Monitoring
# 5) Run
```

Add:

```lammps
# 2) System definition
pair_style lj/cut 4.0
read_data improved.min.data
```

To restore clean initial regions, delete misplaced atoms:

```lammps
region cyl_in cylinder z 0 0 10 INF INF side in
region cyl_out cylinder z 0 0 10 INF INF side out
group grp_t1 type 1
group grp_t2 type 2
group grp_in region cyl_in
group grp_out region cyl_out
group grp_t1_in intersect grp_t1 grp_in
group grp_t2_out intersect grp_t2 grp_out
delete_atoms group grp_t1_in
delete_atoms group grp_t2_out
```

Since LAMMPS has a limited number of custom groups (30), delete groups that are no longer needed:

```lammps
# delete no longer needed groups
group grp_in delete
group grp_out delete
group grp_t1_in delete
group grp_t2_out delete
```

Monitor the number of atoms of each type inside the cylinder using equal‑style variables:

```lammps
variable n1_in equal count(grp_t1,cyl_in)
variable n2_in equal count(grp_t2,cyl_in)
```

Compute a coordination number of type‑2 atoms around type‑1 atoms and its average:

```lammps
compute coor12 grp_t1 coord/atom cutoff 2 group grp_t2
compute sumcoor12 grp_t1 reduce ave c_coor12
```

!!! note
    Compute commands produce data that must be consumed by other commands (e.g., referenced in `thermo_style`, written to dumps, used in variables, or consumed by fixes/other computes).

Complete monitoring and the run:

```lammps
# 4) Monitoring
thermo 1000
thermo_style custom step temp pe ke etotal press v_n1_in v_n2_in c_sumcoor12
dump viz all image 1000 myimage-*.ppm type type shiny 0.1 box no 0.01 view 0 0 zoom 1.8 fsaa yes size 800 800
dump_modify viz adiam 1 1 adiam 2 3 acolor 1 turquoise acolor 2 royalblue backcolor white

# 5) Run
velocity all create 1.0 49284 mom yes dist gaussian
fix mynve all nve
fix mylgv all langevin 1.0 1.0 0.1 10917 zero yes
timestep 0.005
run 300000
```

Here `velocity create` assigns initial velocities for a target temperature of 1.0, zeroing net momentum and using a Gaussian distribution. `zero yes` in the Langevin thermostat makes the total random force zero to prevent center‑of‑mass drift.

!!! note
    A bulk system with periodic boundary conditions is expected to remain in place. In temperature calculations, use \(3N-3\) degrees of freedom (no global translation). Otherwise, the “flying ice cube” syndrome may appear.

Over time, the two populations mix; variables `v_n1_in`, `v_n2_in` and compute `c_sumcoor12` help quantify the progression.

<figure>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="../../../../images/lammps/tutorial1/figures/mixing-vmd-dark.png">
    <img alt="Snapshots of the system during mixing (t=0, 75, 1500)" src="../../../../images/lammps/tutorial1/figures/mixing-vmd-light.png">
  </picture>
  <figcaption>
    Snapshots showing the system at t = 0 (left), t = 75 (middle), and t = 1500 (right). Type 1: small green; Type 2: large cyan.
  </figcaption>
</figure>

<figure>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="../../../../images/lammps/tutorial1/figures/LJ-mixing-dm.png">
    <img alt="Evolution of mixing metrics over time" src="../../../../images/lammps/tutorial1/figures/LJ-mixing.png">
  </picture>
  <figcaption>
    (a) Numbers N1,in and N2,in inside the cylinder over time; (b) average coordination number C1-2 between types 1 and 2.
  </figcaption>
</figure>

## Going further with exercises

### Experiments

Here are suggestions for further experiments that may lead to additional insights into how different systems are configured and how various features function:

- Use a Nosé–Hoover thermostat (`fix nvt`) instead of a Langevin thermostat (`fix nve` + `fix langevin`).
- Omit the energy minimization step before starting the MD simulation using either the Nosé–Hoover or the Langevin thermostat.
- Apply a thermostat to only one type of atoms and observe the temperature for each type separately.
- Append an NVE run (without any thermostat) and observe the energy levels.

!!! info "If you are using LAMMPS–GUI"
    A useful experiment is coloring the atoms in the Slide Show according to an observable, such as their respective coordination numbers. Replace the `dump` and `dump_modify` commands with the following lines, then run again:

    ```lammps
    variable coor12 atom (type==1)*(c_coor12)+(type==2)*-1
    dump viz all image 1000 myimage-*.ppm v_coor12 type &
    shiny 0.1 box no 0.01 view 0 0 zoom 1.8 fsaa yes size 800 800
    dump_modify viz adiam 1 1 adiam 2 3 backcolor white &
    amap -1 2 ca 0.0 4 min royalblue 0 turquoise 1 yellow max red
    ```

    Atoms of type 1 are colored based on the value of `c_coor12`, mapped continuously from turquoise to yellow and red for highest coordination. In the definition of `v_coor12`, atoms of type 2 are assigned a value of −1 and remain their default blue.

### Solve “Lost atoms” error

For this exercise, the following input script is provided:

```lammps
units lj
dimension 3
atom_style atomic
pair_style lj/cut 2.5
boundary p p p

region simulation_box block -20 20 -20 20 -20 20
create_box 1 simulation_box
create_atoms 1 random 1000 341841 simulation_box

mass 1 1
pair_coeff 1 1 1.0 1.0

dump mydmp all atom 100 dump.lammpstrj
thermo 100
thermo_style custom step temp pe ke etotal press

fix mynve all nve
fix mylgv all langevin 1.0 1.0 0.1 1530917
timestep 0.005

run 10000
```

As it is, this input returns a common error:

```bash
ERROR: Lost atoms: original 1000 current 984
```

The goal is to fix the error without using any other command than the ones already present. Only adjust parameter values and/or replicate commands as needed.

!!! info
    This script fails because particles are created randomly in space, some of them overlap, and no energy minimization is performed prior to starting MD.

### Create a demixed dense phase

Starting from one of the inputs created in this tutorial, fine‑tune parameters such as particle numbers and interactions to create a simulation with the following properties:

- The number density is high.
- Types 1 and 2 have the same size.
- Types 1 and 2 demix.

!!! hint
    An easy way to create a dense phase is to allow the box dimensions to relax until the vacuum disappears. Replace `fix nve` with `fix nph`.

<figure>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="../../../../images/lammps/tutorial1/figures/demixing-dark.png">
    <img alt="Demixing of Lennard-Jones binary fluid" src="../../../../images/lammps/tutorial1/figures/demixing-light.png">
  </picture>
  <figcaption>
    Snapshots at different times showing progressive demixing and formation of large demixed domains.
  </figcaption>
</figure>

### From atoms to molecules

Add a bond between particles of type 2 to create dumbbell molecules instead of single particles. Similarly, create a small polymer, i.e., a long chain of particles linked by bonds and angles.

<figure>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="../../../../images/lammps/tutorial1/figures/dumbell-dark.png">
    <img alt="Dumbbell Lennard-Jones molecules" src="../../../../images/lammps/tutorial1/figures/dumbell-light.png">
  </picture>
  <figcaption>
    Dumbbell molecules made of 2 large spheres mixed with smaller particles.
  </figcaption>
</figure>

See also: [video](https://youtu.be/R_oHonOQi68) for the dumbbell system, and [video](https://youtu.be/LfqcfP3ZQcY) for the polymer system.

<figure>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="../../../../images/lammps/tutorial1/figures/polymer-dark.png">
    <img alt="Small polymer made of Lennard-Jones spheres" src="../../../../images/lammps/tutorial1/figures/polymer-light.png">
  </picture>
  <figcaption>
    A single small polymer of 9 large spheres mixed with smaller particles.
  </figcaption>
</figure>

!!! hints
    Use a molecule template to easily insert as many atoms connected by bonds (i.e., molecules) as you want. A molecule template typically begins as follows:

    ```lammps
    2 atoms
    1 bonds

    Coords

    1 0.5 0 0
    2 -0.5 0 0

    (...)
    ```

    A bond section also needs to be added. See the LAMMPS documentation for `molecule` templates: <https://docs.lammps.org/molecule.html>.


>Reference: [lammpstutorials](https://lammpstutorials.github.io/sphinx/build/html/tutorial1/lennard-jones-fluid.html)
