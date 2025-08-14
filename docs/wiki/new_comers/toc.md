---
authors: Zhenghao Wu
title: Guide
---

# New Member Tutorial

Welcome to the Wu Research Group! Everyone needs to familiarize themselves with the rules when joining a new environment. Please follow the checklist below in order to better prepare for joining the group.

## First Few Days
- Register Feishu and join the group (email [Zhenghao](mailto:zhenghao.wu@xjtlu.edu.cn) if you have not yet done it)
> **Note**: The commnucation and management of the group heavily relies on Feishu. We use Feishu for quick messaging, discussions, and online meetings. We use feishu calendar to manage time availability.

- Login in group cloud: [http://QuickConnect.to/wuresearchgroup](http://quickconnect.to/wuresearchgroup)
> **Note**: We use synology cloud to share data and documents.
    - Account: YOURNAME (e.g., zhenghaowu), password: Wu123456789
    - Change your password
    - Create a folder under YOURNAME (e.g., ZhenghaoWu)
    - [A quick guide about basics of synology cloud](https://kb.synology.com/en-global/DSM/tutorial/Quick_Start_Synology_Drive_users)
    - Feel free to ask [Zhenghao](mailto:zhenghao.wu@xjtlu.edu.cn) about the cloud settings
- Data Sharing: please send [Zhenghao](mailto:zhenghao.wu@xjtlu.edu.cn) the following items:
    - A selfie for our group website
    - Your email address for Zotero (where we share and manage literature)
    - If your timetable is available, we can arrange 1-on-1 meeting accordingly

## Personal Workspace

Each new PhD student will be assigned a workspace and a workstation (Windows) for daily research. For further information, please contact **[Zhenghao](mailto:zhenghao.wu@xjtlu.edu.cn)** by email.

> **Note**: Due to limited space, dedicated workspaces are currently not assigned to master's and undergraduate students.

## HPC and HPC account

XJTLU has university-wide [HPC](https://hpc.xjtlu.edu.cn/) for scientific research, which runs on the **Linux** operating system. Like Windows, Linux is another type of computer operating system, but it is mainly operated via the command line (keyboard). Therefore, if you are not familiar with the **Linux** system, please first do some basic self-study using [:fontawesome-brands-bilibili: Introduction to Linux](./linux/2024.md).

To log in to the cluster, you also need a HPC account. Please contact **[Zhenghao](mailto:zhenghao.wu@xjtlu.edu.cn)** to apply for account for you.

It is **recommended to use the terminal on Linux or MacOS** to log in to the cluster. Here, MacOS refers to Apple's operating system. Since both Apple's MacOS and Linux are derived from the Unix system, using an Apple computer to log in to the cluster is the most convenient. For computers running Windows, you will need to install additional software.

To log in to the cluster on MacOS, simply press ++command+space++ to bring up the search bar. Type `terminal` to open the Terminal application. You can then use the `SSH` command in the terminal. For details on using `SSH`, see below.

For Windows users, please refer to [Introduction to WSL on Windows](../cluster_usage/wsl_usage.md) and [:fontawesome-brands-bilibili: Introduction to Linux](./linux/2024.md) for instructions on configuring your system to achieve a seamless terminal experience.

To create an account, you need to generate an [SSH key pair](../cluster_usage/ssh_note.md#create-key-pair). Logging in to the cluster requires using [SSH](../cluster_usage/ssh_note.md).

Before using the cluster, please familiarize yourself with the [basic knowledge](../cluster_usage/cluster_usage.md) and operations. If you need to use GPU resources, please also learn [how to use GPUs on the cluster](../cluster_usage/gpu_usage.md).

If you find any of the above content difficult to understand, please report it immediately to **[Zhenghao](mailto:zhenghao.wu@xjtlu.edu.cn)**.

## Using Python on Linux/MacOS and HPC

Python is a very convenient programming language that can help us process computational data. However, installing pure Python and the required Python libraries can be quite troublesome. Therefore, a software called `Anaconda` can help us solve this problem.

On MacOS, to install `Anaconda`, simply search for "Anaconda" in your preferred search engine and download the appropriate installer from the official website.

On the cluster, `Anaconda` has already been installed for everyone. For instructions on how to use and set it up, please refer to [Anaconda on the Cluster](../cluster_usage/conda.md).

## Required Learning Items
> **Note**: some of these materials can be found on our cloud.

<div class="annotate" markdown>

- [:fontawesome-solid-book: Simulations: the dark side](https://arxiv.org/pdf/1211.4440v1)
- [:fontawesome-solid-book: Understanding Molecular Dynamics (First 3 chapters)](../book_recommendation.md)
- [:fontawesome-brands-bilibili: Introduction to Linux (including Windows environment setup + Git introduction)](./linux/2024.md) 
- [:fontawesome-brands-bilibili: How to Read a Paper](https://web.stanford.edu/class/ee384m/Handouts/HowtoReadPaper.pdf)
- [:fontawesome-brands-bilibili: Python Tutorial](./python/2024.md)
- [:fontawesome-brands-bilibili: ASE & Packmol Model Build Basics](./tools/2024-ase.md)

</div>

## Optional Learning (Project-specific)

### Statistical Mechanics

- [:material-file-multiple: Introduction to Statistical Mechanics](https://web.stanford.edu/~peastman/statmech/#contents) 
- [:material-file-multiple: David Tong at DAMTP, Cambridge: Lectures on Theoretical Physics](http://www.damtp.cam.ac.uk/user/tong/teaching.html) 
- [:material-file-multiple: Lectures on Statistical Physics](https://www.damtp.cam.ac.uk/user/tong/statphys.html) 


### Molecular Dynamics

- [:fontawesome-solid-book: Foundations in Molecular Simulations](https://livecomsjournal.org/index.php/livecoms/article/view/v1i1e5957)
- [:fontawesome-solid-book: Quantification of Uncertainty and Sampling Quality in Molecular Simulations](https://livecomsjournal.org/index.php/livecoms/article/view/v1i1e5067)
- [:fontawesome-solid-book: Understanding Molecular Dynamics by Smit and Frenkel (First 3 Chapters)](../book_recommendation.md)
- [:fontawesome-brands-bilibili: 分子动力学/机器学习分子动力学实践](./md/2024-md.md)

### Machine Learning

<div class="annotate" markdown>

- [:fontawesome-solid-book: Deep Learning for Molecules & Materials](https://dmol.pub/index.html) 
- [:fontawesome-solid-book: Pattern Recognition and Machine Learning](https://www.microsoft.com/en-us/research/uploads/prod/2006/01/Bishop-Pattern-Recognition-and-Machine-Learning-2006.pdf) 
- [:fontawesome-solid-book: Deep Learning by Goodfellow, Bengio, and Courville](http://alvarestech.com/temp/deep/Deep%20Learning%20by%20Ian%20Goodfellow,%20Yoshua%20Bengio,%20Aaron%20Courville%20(z-lib.org).pdf) 
- [:fontawesome-brands-youtube: Deep Learning Lecture by Frank Noe](https://www.youtube.com/playlist?list=PLqPI2gxxYgMKN5AVcTajQ79BTV4BiFN_0)
- [:fontawesome-brands-youtube: Machine Learning for Physics and the Physics of Learning 2019](https://www.youtube.com/playlist?list=PLHyI3Fbmv0SfQfS1rknFsr_UaaWpJ1EKA)
- [:fontawesome-brands-bilibili: 人工智能技术入门](./ai/2024-ai.md)
- [:fontawesome-brands-bilibili: AI 模型训练](./ai/2024-train.md)

</div>

### Machine Learning Force Fields

- [:material-file-multiple: A practical guide to machine learning interatomic potentials](https://pdf.sciencedirectassets.com/272101/1-s2.0-S1359028625X00024/1-s2.0-S1359028625000014/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEN3%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIQCytEsAip2qstC%2F%2BPOkn4sVQHYbe0%2BvWhmBhBjtilDtaAIgda4IWtcpn3HZdzyfGAPvAOPoApTkrghMnN1SbGOXmpEqsgUIJRAFGgwwNTkwMDM1NDY4NjUiDL4eYLAcGEqw2vmZVSqPBf25%2Ffl4Na86og7b0z%2FaIddSfVDb4konnUlexm1CPvLD4McTdHt14US95gLr4TM4JVrvRJdjZwlLCkAOul9k%2FEKV4ZD9V3z3F%2BkRWZOZCVtjaWPbOqnJt5zumfa5XiKCgPdnpQjtOoGA0LvybxXkVRPqakxAg4xJG%2FQ3XKu8XSGss66aFdcvtCyhoKrUKJ45vyXzCvmHdCuprywCW5ZI%2BoBolGA70XxF2q2KyMpRQ4Vmt%2Bmew5WU4XdP8KxmeCXseMDU%2FWKdhtnyEfS9nhh3UIYVXEnvXYr6MCcbRBXya0GpqY5AJxmp0C5QTW37qjNxdPSZJEmDGrRJ6j2d2jjbDxqEtdAO8tjhSqVzmOmDasVBZyyyfnLsdezvIAmqAph8Lh3tG9SA4CkvWzMuxho2a83aw5Qu%2BVN%2FY%2Fe%2F%2ByVMz%2B786Vy9rKSI2G0qtE0sdEuQq%2BqxZI1cwE962xzvMCNNF8uz3cr8rJ4i2nApdQb3h00NtidZA0b62NtL6bJnY3T32ZSuh0Bya0ykRiD4ag2iO2wUbU%2B9kFQT1YFJw7O9qUfgrlp300ys%2FhJPAKNESMxmj2WkOyU%2BHtrYIB6HtJNwKGUGxcpHq%2FJwQWexhqahbpfkFBvoizX1BBzJ0gq8Lkj2R3G%2BjdYd0wpUGA2k6X6%2BWnOLk8wEVQ%2FLPtdWY9k176H2zsxWXLJdSp04VSBKPdSoplMnv2FbjB%2FH2oQ8D4K2z30ryOSgZYrQzCWeCHJA8Ep5ozMJHPdDGXeLqEP1uijYLvYDTKZ7Px8gSvhu0mshcrglMv5MtZffIzR1sxpIjsCjs1DdOusbX2m4537EdzbzyoA3s7a7olmBa6MtdTPPnbfeL6A5vwThc8MEMLSLgRswk6HwxAY6sQGlfeK4KUZ4n%2BS%2FCFvhrvVlm9C4h1zetiX%2FlVdF1Y%2FEgDL2iVtIyrAIwUzSTKCzNSRmA%2FDAXV9iOUhsnM4%2Bdp%2F%2BuSrj4ES%2F35i387h%2BgEUGqq9vpnrSaAexO%2FakMFkiIKuFA%2FHs9V6wTG82KmW0ZcS%2B01eEIVm8qlh%2FHh%2BkxGDxFtGekolT%2FIVyZ9AyGqEdGdzIo2%2BLfw3gkkN7B13zrRJ9j3Qk8KW7bSw8dZf14JwKzrA%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20250813T050522Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYV3CX7WE6%2F20250813%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=e095ae987a4cc924c685ee52ce35fffb11c630786af8bdbdfb8a9ccfa57545ea&hash=060ecfbfc4a8135b7c963f314400e455ad00b581498c18351fe468bb1b7f02c4&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S1359028625000014&tid=spdf-b4f509c4-7431-41c8-9793-54f0f5d1aa84&sid=4f7e78e89550d04b7789cb74737aaba9ddbfgxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&rh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0b0c5856035a070150&rr=96e5ad516f4609b8&cc=hk&kca=eyJrZXkiOiIvem1NbFNVMk9nR3pZdGRjbWdwTjA1cURqc0JUckxDSHVwWG10QUlZQUpyN2d4dXh6akRTMlZ1RXFKdzlKbUt1cWFQRGF4eTNTYmF0MGg2SGZRSFZITWJ3K2tNVXdyckI5QVNVNHNwRkNoSXk4SHQza1V5RnpWek9RSENEODJQM1VGQkhaNFJEYi9lRlNpU1FNSVRtRU5SckxPZWRmOG9pR0JvcGg3RElmdG1xL0YzWVV3PT0iLCJpdiI6IjU2NDE1MDdhNTVlNjZkMTg3MTE2NDUxNzZmODMzZTE3In0=_1755061528013)
- [:fontawesome-brands-youtube: A General Introduction of MLFF by Gabor Csanyi at IPAM](https://www.youtube.com/watch?v=JxIFxM1f40U)
- [:fontawesome-brands-bilibili: 分子动力学/机器学习分子动力学实践](./md/2024-md.md)


### Polymer Physics and Chemistry
- [:fontawesome-solid-book: Polymer Physics by Rubinstein and Colby](../book_recommendation.md) 
- [:fontawesome-solid-book: The Theory of Polymer Dynamics by Edwards and Doi](../book_recommendation.md) 
- [:fontawesome-solid-book: Scaling concepts in polymer physics](../book_recommendation.md) 


!!! info "说明"
    - :fontawesome-solid-book: : 书籍或文献资料等
    - :fontawesome-brands-youtube: :fontawesome-brands-bilibili: : 视频教程
    - :material-file-multiple: : 博客或文档等
