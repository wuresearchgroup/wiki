---
title: HPC XJTLU Tutorial
authors: Zhangtao Yi, Zhenghao Wu
comments: false
---
Reference: [https://hpc.xjtlu.edu.cn/](https://hpc.xjtlu.edu.cn/)

## 1. Introduction of HPC

**What is HPC (high performance computing)?**
HPC is the technology that uses clusters of powerful processors, working in parallel, to process massive multi-dimensional datasets (big data) and solve complex problems at extremely high speeds.

### Hardware configuration of the HPC cluster

| System | Partition | Nodes | Processor | Cores/ Socket | Threads/ Core | Memory/ Node (GB) | GPUs/ Node |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **cpu6348** | cpu6348n[1–2] | 2 | Intel(R) Xeon(R) Gold 6348 | 28 | – | 512 | – |
| **gpu3090** | gpu3090n[1–6] | 6 | 2*Intel(R) Xeon(R) Gold 6336Y | 24 | 2 | 256 | 8*RTX 3090(24G) |

---

## 2. Log in

**Accessing the Login Node**
We can use SSH to access the login node of the HPC cluster (e.g., via Remote Explorer plug-in in Visual Studio Code).

**Local SSH Configuration File**
Default path: `C:\Users\username\.ssh\config`

```ssh
Host HPC
    HostName 10.7.88.95
    Port 22
    User (Type your user name, e.g. zhangtaoyi24)
    ForwardAgent yes
    ForwardX11 yes
```

> **Note:** Please note that the capacity limit of the `/gpfs/home/` personal directory is **10GB**. It is recommended to use the `/gpfs/work/` personal directory for storage and save data in a timely manner.

---

## 3. How to use HPC

**Specify the computational parameters in a job script:**
We use Slurm as the job scheduling and management system. QoS (Quality of Service) and Partition must be specified in the job script, which is written in a `.sh` file and begins with `#!/bin/bash`.

```bash
#!/bin/bash

#SBATCH --qos 8gpus
#SBATCH --partition gpu3090
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-12:00:00
#SBATCH --mem-per-gpu=20G
#SBATCH --mail-type=All
#SBATCH --job-name=GPU_TEST
#SBATCH --error=%j.err
#SBATCH --output=%j.out
```

*   We use `time` to specify the maximum time that you can use. Recall to complete your computation and save the results before endtime.
*   You can change the "error" and "output" parameters to specify your log directory.
*   A larger request for computational resources will downgrade the priority of the job.
*   QOS can be viewed by `sacctmgr` (see below "Available QOS").

**Write down your executing codes in the job script:**

```bash
module load openmpi/4.1.5-gcc-8.5.0-2wx4wru
cd /gpfs/work/che/zhangtaoyi24/MLFF_small/methanol/Lammps/
nvidia-smi -l 300 > gpu_usage_NVTprod.log &
nohup bash -c 'while true; do echo "$(date "+%F %T") CPU: $(top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk "{print 100 - \$1\"%\"}")" >> cpu_usage_NVTprod.log; sleep 300; done' &

mpirun -np 1 lmp \
-sf gpu \
-pk gpu 1 neigh yes binsize 6 split 1 \
-in system_NVTprod.inp > output_NVTprod.log 2>&1
```
*   *We use `module load` to load software (https://lmod.readthedocs.io/en/latest/010_user.html).*

**Load Python virtual environment (Alternatively):**
Sometimes we would like to employ packages in Python virtual environment, especially when analyzing data.

```bash
source /gpfs/work/che/zhangtaoyi24/miniconda3/etc/profile.d/conda.sh
conda activate mlff
cd /gpfs/work/che/zhangtaoyi24/MLFF_small/methanol/cgmap/

cgmap --dump ../5nsNVT_prod.dump --system system.yaml --output cg_system > output_cgmap.log 2>&1
```

**Submit your job script:**
Save your job script and use `sbatch` command to submit it in the login node.

```bash
(base) [zhangtaoyi24@xpszlogin1 ~]$ sbatch sh/mlff_small_methanol.sh
```

**View your job state:**
We can use `squeue` command to view the state of submitted job. The job can be at running (R) or pending (PD) and other states.

```bash
(base) [zhangtaoyi24@xpszlogin1 ~]$ squeue -u $USER
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
```

> **! Never** run computing jobs at the login node. Computing jobs should be submitted to the computing nodes via `sbatch` command. Installation of software and Python packages can be run at the login node.

---

## 4. Storage

**The storage we have:**
Parallel Storage System
*   GPFS

**Home Directory**
*   Default login directory
*   10GB
*   Example: `/gpfs/home/che/zhangtaoyi24/`

**Work Directory**
*   Large file storage
*   Default 500GB
*   Example: `/gpfs/work/che/zhangtaoyi24/`

---

## 5. Software Management

**Module Management:**
We can use `module avail` command to view all available software.

```text
(base) [zhangtaoyi24@xpszlogin1 ~]$ module avail

-------------------- /gpfs/spack/modules/linux-rocky8-icelake -----------------------
abaqus/2019-gcc-8.5.0-jjwvu43           py-contourpy/1.0.7-gcc-8.5.0-tzsepyi
adgpu/1.6-gcc-8.5.0-ymcaenf             py-contourpy/1.0.7-gcc-9.5.0-slsabzj(D)
adgpu/1.6-gcc-8.5.0-3yd4vk5(D)          py-cppy/1.1.0-gcc-8.5.0-isynyft
adios2/2.9.0-gcc-8.5.0-4nle5cs          py-cppy/1.2.1-gcc-8.5.0-ah5ii5k
alsa-lib/1.2.3.2-gcc-8.5.0-mq2xqit      py-cppy/1.2.1-gcc-8.5.0-h3dpfkz
amber/22-gcc-8.5.0-kzejyck              py-cppy/1.2.1-gcc-9.5.0-pfpumnr(D)
anaconda3/2022.10-gcc-8.5.0-4dp3trd     py-cutadapt/4.4-gcc-8.5.0-cc55eb5
anicalculator/1-gcc-8.5.0-nzcdo4c       py-cycler/0.11.0-gcc-8.5.0-uiyuo52
anicalculator/1-gcc-9.5.0-xbqsnem(D)    py-cycler/0.11.0-gcc-8.5.0-uq42bi4
ansys/2022R1-gcc-8.5.0-6opggez          py-cycler/0.11.0-gcc-8.5.0-76o2wsp
ant/1.10.13-gcc-8.5.0-owgwryr           py-cycler/0.11.0-gcc-9.5.0-4rqppp2(D)
...
--More--
```

**Conda environment Management:**
We can use `conda env list` command to view all available conda environment.

```text
(base) [zhangtaoyi24@xpszlogin1 ~]$ conda env list

# conda environments:
#
base                  * /gpfs/work/che/zhangtaoyi24/miniconda3
mlff                    /gpfs/work/che/zhangtaoyi24/miniconda3/envs/mlff
```

---

## 6. Available Cluster status and Job Scheduling

**Available QOS**

```bash
(mlff) [zhangtaoyi24@xpszlogin1 ~]$ sacctmgr show assoc user=$USER format=User,Account,Cluster,QOS

      User    Account    Cluster                  QOS
---------- ---------- ---------- --------------------
zhangtaoy+ zhenghaowu    xipuhpc 104cores,8gpus,cpud+
```

**View job status**
We can use `scontrol` to see the status of submitted job (e.g. start time, priority)

```bash
(base) [zhangtaoyi24@xpszlogin1 ~]$ scontrol show job (your jobid)
```

**Cancel job**

```bash
(base) [zhangtaoyi24@xpszlogin1 ~]$ scancel (your jobid)
```

For more information, please visit [https://hpc.xjtlu.edu.cn](https://hpc.xjtlu.edu.cn).