![header](imgs/header.jpg)

# AlphaFold

<p align="center">
  <a href="README.md">English</a> · <b>简体中文</b>
</p>

本代码库提供了 AlphaFold v2 推理流水线的开源实现。为简便起见，在本文档的其余部分统称为 AlphaFold。

我们还提供：

1.  **AlphaFold-Multimer 实现**：目前仍在持续完善中，AlphaFold-Multimer 的稳定性可能暂不及我们的单体（monomer）AlphaFold 系统。[请阅读指南](#更新现有安装)了解如何升级和更新代码。
2.  **[技术备忘录](docs/technical_note_v2.3.0.md)**：包含升级版 AlphaFold v2.3.0 的模型架构与推理流程。
3.  **[CASP15 基线预测集](docs/casp15_predictions.zip)**：包含预测结果及所有人工干预步骤的文档记录。

任何发表基于使用本源代码或模型参数的研究成果，均应[引用](#引用本研究) [AlphaFold 论文](https://doi.org/10.1038/s41586-021-03819-2)；如适用，还应引用 [AlphaFold-Multimer 论文](https://www.biorxiv.org/content/10.1101/2021.10.04.463034v1)。

关于该方法的详细技术描述，请参阅[补充材料（Supplementary Information）](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03819-2/MediaObjects/41586_2021_3819_MOESM1_ESM.pdf)。

**您也可以通过社区支持的版本使用略微简化的 AlphaFold（见下文）。**

如果您有任何疑问，请通过 [alphafold@deepmind.com](mailto:alphafold@deepmind.com) 联系 AlphaFold 团队。

![CASP14 predictions](imgs/casp14_predictions.gif)

## 安装与运行您的首次预测

您需要一台运行 Linux 的机器，AlphaFold 不支持其他操作系统。完整安装最多需要 3 TB 的磁盘空间来存储基因数据库（推荐使用固态硬盘 SSD），以及一块现代 NVIDIA GPU（显存更大的 GPU 可以预测更大的蛋白质结构）。

请按照以下步骤操作：

1.  安装 [Docker](https://www.docker.com/)。

    *   安装 [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) 以获取 GPU 支持。
    *   配置[以非 root 用户身份运行 Docker](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user)。

2.  克隆本仓库并进入目录：

    ```bash
    git clone https://github.com/deepmind/alphafold.git
    cd ./alphafold
    ```

3.  下载基因序列数据库与模型参数：

    *   安装 `aria2c`。在大多数 Linux 发行版中均可通过包管理器获取 `aria2` 软件包（在基于 Debian 的系统上，运行 `sudo apt install aria2` 即可安装）。`rsync` 同理。

    *   请使用脚本 `scripts/download_all_data.sh` 下载并配置完整数据库。此过程可能耗时较长（下载数据量约为 556 GB），建议在后台运行该脚本：

    ```bash
    scripts/download_all_data.sh <DOWNLOAD_DIR> > download.log 2> download_all.log &
    ```

    *   **注意：下载目录 `<DOWNLOAD_DIR>` *切勿* 设为 AlphaFold 仓库目录下的子目录。** 否则在构建 Docker 镜像时，庞大的数据库会被复制进 Docker 构建上下文，导致构建极其缓慢。

    *   也可以使用精简版数据库运行 AlphaFold；详情请参阅[完整说明](#基因序列数据库)。

4.  运行以下命令，验证 AlphaFold 是否能够正常调用 GPU：

    ```bash
    docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi
    ```

    该命令的输出应列出您的 GPU 设备列表。如果未显示，请检查设置 [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) 时的所有步骤是否正确，或参考此 [NVIDIA Docker Issue](https://github.com/NVIDIA/nvidia-docker/issues/1447#issuecomment-801479573)。

    如果您希望使用 Singularity（超算集群 HPC 系统上常用的容器化平台）运行 AlphaFold，我们建议参考第三方 Singularity 配置方案，相关链接见 https://github.com/deepmind/alphafold/issues/10 或 https://github.com/deepmind/alphafold/issues/24。

5.  构建 Docker 镜像：

    ```bash
    docker build -f docker/Dockerfile -t alphafold .
    ```

    如果遇到如下错误：

    ```
    W: GPG error: https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64 InRelease: The following signatures couldn't be verified because the public key is not available: NO_PUBKEY A4B469963BF863CC
    E: The repository 'https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64 InRelease' is not signed.
    ```

    请参考 https://github.com/deepmind/alphafold/issues/463#issuecomment-1124881779 中描述的权宜方案。

6.  安装 `run_docker.py` 的依赖项。注意：您也可以选择创建 [Python 虚拟环境（Virtual Environment）](https://docs.python.org/3/tutorial/venv.html) 以避免与系统 Python 环境冲突。

    ```bash
    pip3 install -r docker/requirements.txt
    ```

7.  确保输出目录已存在（默认路径为 `/tmp/alphafold`），且您对其拥有足够的写入权限。

8.  运行 `run_docker.py`，指向包含目标预测蛋白质序列的 FASTA 文件（`--fasta_paths` 参数）。AlphaFold 会检索由 `--max_template_date` 参数指定日期之前可用的模板；此功能可用于在建模过程中排除某些模板。`--data_dir` 为已下载基因数据库的所在目录，`--output_dir` 为输出目录的绝对路径。

    ```bash
    python3 docker/run_docker.py \
      --fasta_paths=your_protein.fasta \
      --max_template_date=2022-01-01 \
      --data_dir=$DOWNLOAD_DIR \
      --output_dir=/home/user/absolute_path_to_the_output_dir
    ```

9.  运行结束后，输出目录中将包含目标蛋白质的预测结构。更多详细选项与故障排查建议，请查阅下文文档。

### 基因序列数据库

此步骤需要在您的机器上安装 `aria2c`。

AlphaFold 运行需要多个基因（序列）数据库：

*   [BFD](https://bfd.mmseqs.com/)，
*   [MGnify](https://www.ebi.ac.uk/metagenomics/)，
*   [PDB70](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/)，
*   [PDB](https://www.rcsb.org/)（mmCIF 格式结构），
*   [PDB seqres](https://www.rcsb.org/) – 仅用于 AlphaFold-Multimer，
*   [UniRef30 (原 UniClust30)](https://uniclust.mmseqs.com/)，
*   [UniProt](https://www.uniprot.org/uniprot/) – 仅用于 AlphaFold-Multimer，
*   [UniRef90](https://www.uniprot.org/help/uniref)。

我们提供了脚本 `scripts/download_all_data.sh`，可用于下载并配置所有这些数据库：

*   推荐默认用法：

    ```bash
    scripts/download_all_data.sh <DOWNLOAD_DIR>
    ```

    将下载全量完整数据库。

*   携带 `reduced_dbs` 参数：

    ```bash
    scripts/download_all_data.sh <DOWNLOAD_DIR> reduced_dbs
    ```

    将下载用于 `reduced_dbs` 数据库预设的精简版数据库。在后续 AlphaFold 运行中，应配合对应的参数 `--db_preset=reduced_dbs` 使用（详见[运行 AlphaFold](#运行-alphafold) 章节）。

:ledger: **注意：下载目录 `<DOWNLOAD_DIR>` *切勿* 设为 AlphaFold 仓库目录下的子目录。** 否则在镜像创建过程中，庞大的数据库会被复制，导致 Docker 构建异常缓慢。

我们未完全提供 CASP14 中使用的数据库版本 – 详情请参阅[关于 CASP14 可复现性的说明](#关于-casp14-可复现性的说明)。为了提高速度，部分数据库已建立了镜像分发，详见[镜像数据库](#镜像数据库)。

:ledger: **注意：完整数据库的下载总数据量约为 556 GB，解压后的总大小约为 2.62 TB。请确保拥有足够的硬盘空间、网络带宽和下载时间。建议使用 SSD 以获得更好的基因比对搜索性能。**

:ledger: **注意：若下载目录和数据集缺少完整的读写权限，可能导致 MSA 比对工具报错，并出现难以诊断的（外部）错误信息。请确保配置了相应权限，例如执行 `sudo chmod 755 --recursive "$DOWNLOAD_DIR"` 命令。**

`download_all_data.sh` 脚本还会一并下载模型参数文件。脚本执行完成后，目录结构应如下所示：

```
$DOWNLOAD_DIR/                             # 总计：~ 2.62 TB (下载大小：556 GB)
    bfd/                                   # ~ 1.8 TB (下载大小：271.6 GB)
        # 6 个文件。
    mgnify/                                # ~ 120 GB (下载大小：67 GB)
        mgy_clusters_2022_05.fa
    params/                                # ~ 5.3 GB (下载大小：5.3 GB)
        # 5 个 CASP14 模型，
        # 5 个 pTM 模型，
        # 5 个 AlphaFold-Multimer 模型，
        # LICENSE，
        # = 共计 16 个文件。
    pdb70/                                 # ~ 56 GB (下载大小：19.5 GB)
        # 9 个文件。
    pdb_mmcif/                             # ~ 238 GB (下载大小：43 GB)
        mmcif_files/
            # 约 199,000 个 .cif 文件。
        obsolete.dat
    pdb_seqres/                            # ~ 0.2 GB (下载大小：0.2 GB)
        pdb_seqres.txt
    small_bfd/                             # ~ 17 GB (下载大小：9.6 GB)
        bfd-first_non_consensus_sequences.fasta
    uniref30/                              # ~ 206 GB (下载大小：52.5 GB)
        # 7 个文件。
    uniprot/                               # ~ 105 GB (下载大小：53 GB)
        uniprot.fasta
    uniref90/                              # ~ 67 GB (下载大小：34 GB)
        uniref90.fasta
```

只有在下载完整数据库时才会下载 `bfd/`，而在下载精简数据库时只会下载 `small_bfd/`。

### 模型参数

虽然 AlphaFold 代码采用 Apache 2.0 许可证开源，但 AlphaFold 模型参数和 CASP15 预测数据是依据 CC BY 4.0 许可证条款提供的。有关详细信息，请参阅下文的[免责声明](#许可证与免责声明)。

AlphaFold 参数可从 https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar 获取，并在运行 `scripts/download_all_data.sh` 脚本时自动下载。该脚本将下载以下模型参数：

*   **5 个 CASP14 模型**：在 CASP14 期间使用，经过广泛的结构预测质量验证（详情请参阅 Jumper 等人 2021，Suppl. Methods 1.12）。
*   **5 个 pTM 模型**：经过微调，在输出结构预测的同时，生成 pTM（预测 TM-score）和 PAE（预测对齐误差）值（详情请参阅 Jumper 等人 2021，Suppl. Methods 1.9.7）。
*   **5 个 AlphaFold-Multimer 模型**：在输出复合物结构预测的同时生成 pTM 和 PAE 值。

### 更新现有安装

如果您已安装了旧版本，既可以从头完全重装（删除所有内容并重新执行配置），也可以进行增量更新。增量更新速度显著更快，但需要稍多的操作。请务必严格按照以下顺序执行：

1.  **更新代码。**
    *   进入已克隆的 AlphaFold 仓库目录，运行 `git fetch origin main` 获取所有代码更新。
2.  **更新 UniProt、UniRef、MGnify 和 PDB seqres 数据库。**
    *   删除 `<DOWNLOAD_DIR>/uniprot`。
    *   运行 `scripts/download_uniprot.sh <DOWNLOAD_DIR>`。
    *   删除 `<DOWNLOAD_DIR>/uniclust30`。
    *   运行 `scripts/download_uniref30.sh <DOWNLOAD_DIR>`。
    *   删除 `<DOWNLOAD_DIR>/uniref90`。
    *   运行 `scripts/download_uniref90.sh <DOWNLOAD_DIR>`。
    *   删除 `<DOWNLOAD_DIR>/mgnify`。
    *   运行 `scripts/download_mgnify.sh <DOWNLOAD_DIR>`。
    *   删除 `<DOWNLOAD_DIR>/pdb_mmcif`。PDB SeqRes 和 PDB 必须来自完全相同的日期。若忽略此步骤，在运行 AlphaFold-Multimer 进行模板搜索时可能会产生错误。
    *   运行 `scripts/download_pdb_mmcif.sh <DOWNLOAD_DIR>`。
    *   运行 `scripts/download_pdb_seqres.sh <DOWNLOAD_DIR>`。
3.  **更新模型参数。**
    *   删除 `<DOWNLOAD_DIR>/params` 中的旧模型参数。
    *   使用 `scripts/download_alphafold_params.sh <DOWNLOAD_DIR>` 下载最新模型参数。
4.  **遵循[运行 AlphaFold](#运行-alphafold)步骤。**

#### 使用已弃用的旧版模型权重

如需使用已弃用的 v2.2.0 AlphaFold-Multimer 模型权重：

1.  将 `scripts/download_alphafold_params.sh` 中的 `SOURCE_URL` 修改为 `https://storage.googleapis.com/alphafold/alphafold_params_2022-03-02.tar`，并下载旧版参数。
2.  将 `config.py` 中多聚体 `MODEL_PRESETS` 里的 `_v3` 修改为 `_v2`。

如需使用已弃用的 v2.1.0 AlphaFold-Multimer 模型权重：

1.  将 `scripts/download_alphafold_params.sh` 中的 `SOURCE_URL` 修改为 `https://storage.googleapis.com/alphafold/alphafold_params_2022-01-19.tar`，并下载旧版参数。
2.  移除 `config.py` 中多聚体 `MODEL_PRESETS` 里的 `_v3`。

## 运行 AlphaFold

**运行 AlphaFold 最简单的方式是使用提供的 Docker 脚本。** 该方案已在 Google Cloud 虚拟机上进行了充分验证，测试环境使用 `nvidia-gpu-cloud-image`、12 个 vCPU、85 GB RAM、100 GB 系统盘、挂载在额外 3 TB 磁盘上的数据库，以及一块 A100 GPU。首次运行请遵循[安装与运行您的首次预测](#安装与运行您的首次预测)部分的指引。

1.  默认情况下，AlphaFold 将尝试使用所有可见的 GPU 设备。若要仅使用指定的部分 GPU，可通过 `--gpu_devices` 标志指定逗号分隔的 GPU UUID 或索引列表。详情请参阅 [GPU 枚举说明](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/user-guide.html#gpu-enumeration)。

2.  您可以通过添加 `--model_preset=` 标志来指定运行哪种 AlphaFold 模型。我们提供以下模型预设：

    *   **monomer**：CASP14 中使用的原始模型，不采用集成（ensembling）。
    *   **monomer\_casp14**：CASP14 中使用的原始模型，配置 `num_ensemble=8`，与我们在 CASP14 中的设置完全一致。主要用于结果复现；其计算开销增加 8 倍，但精度提升有限（CASP14 结构域上平均 GDT 增加约 0.1）。
    *   **monomer\_ptm**：使用 pTM 头部微调的原始 CASP14 模型，可输出残基对级别的置信度评估。其精度略低于标准的 monomer 模型。
    *   **multimer**：[AlphaFold-Multimer](#引用本研究) 模型。使用此模型需提供包含多条序列的 FASTA 文件，且必须已下载 UniProt 数据库。

3.  您可以通过在运行命令中添加 `--db_preset=reduced_dbs` 或 `--db_preset=full_dbs` 来平衡多序列比对（MSA）的速度与质量。我们提供以下预设：

    *   **reduced\_dbs**：针对速度和较低硬件要求进行了优化。使用精简版 BFD 数据库运行。需要 8 个 CPU 核心（vCPU）、8 GB RAM 和 600 GB 磁盘空间。
    *   **full\_dbs**：使用 CASP14 中所用的全量基因数据库运行。

    使用 `monomer` 模型预设和 `reduced_dbs` 数据预设运行上述命令示例如下：

    ```bash
    python3 docker/run_docker.py \
      --fasta_paths=T1050.fasta \
      --max_template_date=2020-05-14 \
      --model_preset=monomer \
      --db_preset=reduced_dbs \
      --data_dir=$DOWNLOAD_DIR \
      --output_dir=/home/user/absolute_path_to_the_output_dir
    ```

4.  在生成预测模型后，AlphaFold 会执行松弛（relaxation）步骤以改善局部几何结构。默认情况下，仅对（基于 pLDDT 得分的）最佳模型进行松弛处理（`--models_to_relax=best`），也可以对所有模型（`--models_to_relax=all`）或不对任何模型（`--models_to_relax=none`）进行松弛。

5.  松弛步骤可在 GPU（速度更快，但稳定性稍弱）或 CPU（较慢，但更稳定）上执行。这可以通过 `--enable_gpu_relax=true`（默认）或 `--enable_gpu_relax=false` 来控制。

6.  AlphaFold 支持通过 `--use_precomputed_msas=true` 选项复用相同序列已计算好的多序列比对（MSA）；这在测试不同的 AlphaFold 运行参数时非常有用。该选项假定输出目录中已存在首次 AlphaFold 运行生成的目录结构，且蛋白质序列保持一致。

### 运行 AlphaFold-Multimer

所有操作步骤与运行单体系统相同，但您需要：

*   提供包含多条序列的输入 FASTA 文件；
*   设置 `--model_preset=multimer`。

折叠蛋白质复合物 `multimer.fasta` 的命令示例：

```bash
python3 docker/run_docker.py \
  --fasta_paths=multimer.fasta \
  --max_template_date=2020-05-14 \
  --model_preset=multimer \
  --data_dir=$DOWNLOAD_DIR \
  --output_dir=/home/user/absolute_path_to_the_output_dir
```

默认情况下，多聚体系统将为每个模型运行 5 个随机种子（共产生 25 个预测结构）。如果允许精度略微降低，您可以选择每个模型仅运行单个种子。这可以通过 `--num_multimer_predictions_per_model` 标志来设置，例如将其设为 `--num_multimer_predictions_per_model=1`。

### AlphaFold 预测速度

下表列出了不同长度蛋白质的预测运行耗时。我们仅测量带有 3 次循环回收迭代（recycles）的未松弛（unrelaxed）结构预测时间，不包含 MSA 构建与模板搜索耗时。当携带 `--benchmark=true` 运行 `docker/run_docker.py` 时，该耗时将记录在 `timings.json` 中。所有测试耗时均基于单块 NVIDIA A100 GPU。在 A100 上预测较小结构时，可以通过增大 `alphafold/model/config.py` 中的 `global_config.subbatch_size` 来提升预测速度。

氨基酸残基数 | 预测耗时 (秒)
-----------: | ------------------:
100          | 4.9
200          | 7.7
300          | 13
400          | 18
500          | 29
600          | 36
700          | 53
800          | 60
900          | 91
1,000        | 96
1,100        | 140
1,500        | 280
2,000        | 450
2,500        | 969
3,000        | 1,240
3,500        | 2,465
4,000        | 5,660
4,500        | 12,475
5,000        | 18,824

### 示例

以下是 AlphaFold 在不同场景下的使用示例。

#### 折叠单体蛋白质 (Monomer)

假设我们有一个序列为 `<SEQUENCE>` 的单体。输入 FASTA 文件应为：

```fasta
>sequence_name
<SEQUENCE>
```

然后运行以下命令：

```bash
python3 docker/run_docker.py \
  --fasta_paths=monomer.fasta \
  --max_template_date=2021-11-01 \
  --model_preset=monomer \
  --data_dir=$DOWNLOAD_DIR \
  --output_dir=/home/user/absolute_path_to_the_output_dir
```

#### 折叠同源多聚体 (Homomer)

假设我们要折叠具有相同序列 `<SEQUENCE>` 的 3 个拷贝的同源三聚体。输入 FASTA 文件应为：

```fasta
>sequence_1
<SEQUENCE>
>sequence_2
<SEQUENCE>
>sequence_3
<SEQUENCE>
```

然后运行以下命令：

```bash
python3 docker/run_docker.py \
  --fasta_paths=homomer.fasta \
  --max_template_date=2021-11-01 \
  --model_preset=multimer \
  --data_dir=$DOWNLOAD_DIR \
  --output_dir=/home/user/absolute_path_to_the_output_dir
```

#### 折叠异源多聚体 (Heteromer)

假设我们要折叠一个 A2B3 异源多聚体，即含有 2 个拷贝的 `<SEQUENCE A>` 和 3 个拷贝的 `<SEQUENCE B>`。输入 FASTA 文件应为：

```fasta
>sequence_1
<SEQUENCE A>
>sequence_2
<SEQUENCE A>
>sequence_3
<SEQUENCE B>
>sequence_4
<SEQUENCE B>
>sequence_5
<SEQUENCE B>
```

然后运行以下命令：

```bash
python3 docker/run_docker.py \
  --fasta_paths=heteromer.fasta \
  --max_template_date=2021-11-01 \
  --model_preset=multimer \
  --data_dir=$DOWNLOAD_DIR \
  --output_dir=/home/user/absolute_path_to_the_output_dir
```

#### 依次预测多个单体

假设我们有两个单体文件：`monomer1.fasta` 和 `monomer2.fasta`。

可以通过以下命令依序完成折叠预测：

```bash
python3 docker/run_docker.py \
  --fasta_paths=monomer1.fasta,monomer2.fasta \
  --max_template_date=2021-11-01 \
  --model_preset=monomer \
  --data_dir=$DOWNLOAD_DIR \
  --output_dir=/home/user/absolute_path_to_the_output_dir
```

#### 依次预测多个多聚体

假设我们有两个多聚体文件：`multimer1.fasta` 和 `multimer2.fasta`。

可以通过以下命令依序完成折叠预测：

```bash
python3 docker/run_docker.py \
  --fasta_paths=multimer1.fasta,multimer2.fasta \
  --max_template_date=2021-11-01 \
  --model_preset=multimer \
  --data_dir=$DOWNLOAD_DIR \
  --output_dir=/home/user/absolute_path_to_the_output_dir
```

### AlphaFold 输出文件

输出文件将保存在通过 `run_docker.py` 的 `--output_dir` 标志指定的子目录中（默认路径为 `/tmp/alphafold/`）。输出内容包括计算生成的 MSA、未松弛结构、松弛后结构、按置信度排序后的结构、原始模型输出、预测元数据以及各阶段运行耗时。`--output_dir` 目录的结构如下：

```
<target_name>/
    features.pkl
    ranked_{0,1,2,3,4}.pdb
    ranking_debug.json
    relax_metrics.json
    relaxed_model_{1,2,3,4,5}.pdb
    result_model_{1,2,3,4,5}.pkl
    timings.json
    unrelaxed_model_{1,2,3,4,5}.pdb
    msas/
        bfd_uniref_hits.a3m
        mgnify_hits.sto
        uniref90_hits.sto
```

各输出文件的具体内容如下：

*   `features.pkl` – 一个 `pickle` 文件，包含模型用于生成结构预测的输入特征 NumPy 数组。
*   `unrelaxed_model_*.pdb` – PDB 格式文本文件，包含模型直接输出的原始预测结构。
*   `relaxed_model_*.pdb` – PDB 格式文本文件，包含对未松弛预测结构执行 Amber 松弛处理后的结构（详情参阅 Jumper 等人 2021，Suppl. Methods 1.8.6）。
*   `ranked_*.pdb` – PDB 格式文本文件，包含按模型置信度重新排序后的预测结构。其中 `ranked_i.pdb` 表示置信度排名第 (`i + 1`) 位的预测（因此 `ranked_0.pdb` 置信度最高）。模型置信度采用预测 LDDT（pLDDT）得分进行排名（详情参阅 Jumper 等人 2021，Suppl. Methods 1.9.6）。如果设置了 `--models_to_relax=all`，则所有排名的结构均经过松弛；若设置为 `--models_to_relax=best`，则仅 `ranked_0.pdb` 执行了松弛（其余为未松弛）；若设置为 `--models_to_relax=none`，则所有排名的结构均为未松弛。
*   `ranking_debug.json` – JSON 格式文本文件，包含用于执行模型排名的 pLDDT 数值，以及对应回原始模型名称的映射关系。
*   `relax_metrics.json` – JSON 格式文本文件，包含结构松弛指标（例如剩余的几何冲突违规数量）。
*   `timings.json` – JSON 格式文本文件，记录运行 AlphaFold 流水线各部分所消耗的时间。
*   `msas/` – 包含用于构建输入 MSA 的各项基因比对工具命中结果文件的目录。
*   `result_model_*.pkl` – 一个 `pickle` 文件，包含模型直接生成的各项 NumPy 数组嵌套字典。除结构模块的输出外，还包含辅助输出信息，例如：
    *   **距离图（Distograms）**：`distogram/logits` 包含形状为 `[N_res, N_res, N_bins]` 的 NumPy 数组，`distogram/bin_edges` 包含区间区间的定义。
    *   **各残基 pLDDT 得分**：`plddt` 包含形状为 `[N_res]` 的 NumPy 数组，取值范围为 `0` 到 `100`（`100` 表示置信度最高）。可用于识别高置信度预测的序列区域，或在残基间取均值作为靶标蛋白的整体置信度评分。
    *   **预测 TM-score（仅在使用 pTM 模型时提供）**：`ptm` 字段包含一个标量数值。作为全局结构叠合指标的预测值，此分数旨在评估模型对整体结构域堆叠排布的把握程度。
    *   **预测残基对对齐误差（仅在使用 pTM 模型时提供）**：`predicted_aligned_error` 包含形状为 `[N_res, N_res]` 的 NumPy 数组，取值范围从 `0` 到 `max_predicted_aligned_error`（`0` 表示置信度最高）。可用于可视化结构内部结构域排布的相对位置置信度。

pLDDT 置信度度量值存储在输出 PDB 文件的 B-factor（温度因子）字段中（注意：与真实的 B-factor 不同，pLDDT 数值越大代表质量越好，在用于分子置换等下游任务时需特别注意）。

本代码库已经过测试，能够在 CASP14 测试集上通过 5 个模型预测的 pLDDT 排序达到平均 Top-1 精度（部分 CASP 目标是在早期版本的 AlphaFold 上运行的，部分经历了人工干预；详情请参阅我们后续的研究发表）。某些特定靶标（例如 T1064）在不同随机种子下的单次运行可能具有较大的方差。

## 批量蛋白质推理

提供的推理脚本针对单个蛋白质的结构预测进行了高度优化，它会将神经网络针对具体的序列长度、MSA 和模板尺寸进行专用编译。对于大型蛋白质，编译时间占总运行时间的比例微乎其微；但对于较短的蛋白质或已预先计算好 MSA 的情况，编译耗时占比可能会变得较为明显。在需要批量推理的大规模场景下，可以使用我们的 `make_fixed_size` 函数将输入填充至固定大小，从而减少所需的编译次数。

我们暂未提供开箱即用的批量推理脚本，但基于 `RunModel.predict` 方法并配合并行预计算多序列比对系统，自行构建批量推理流水线相当直观。此外，也可以仅以适度的系统开销重复调用现有脚本。

## 关于 CASP14 可复现性的说明

AlphaFold 对少数蛋白质的预测结果具有较高的运行间方差，并可能受输入数据变动的影响。CASP14 目标 T1064 便是一个典型示例：近期存入数据库的海量 SARS-CoV-2 相关序列显著改变了其 MSA 比对结果。通过模型集成与选择机制（运行 5 个模型并选取置信度最高的一个），这种变异性在一定程度上得到了缓解。

若要尽可能严谨地复现我们在 CASP14 系统中的预测结果，您必须使用我们在 CASP 期间所使用的相同版本数据库。这些版本可能与我们脚本默认下载的当前最新版本不同。

基因序列数据库版本：

*   UniRef90：[v2020_01](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2020_01/uniref/)
*   MGnify：[v2018_12](http://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2018_12/)
*   Uniclust30：[v2018_08](http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/)
*   BFD：[唯一可用版本](https://bfd.mmseqs.com/)

模板数据库版本：

*   PDB：（下载于 2020-05-14）
*   PDB70：[2020-05-13](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/old-releases/pdb70_from_mmcif_200513.tar.gz)

对于模板库，另一种替代方案是使用最新的 PDB 和 PDB70，但附加标志 `--max_template_date=2020-05-14`，以将模板严格限制在 CASP14 开始前已公开发布的结构中。

## 引用本研究

如果您在学术研究中使用了本软件包的代码或数据，请引用：

```bibtex
@Article{AlphaFold2021,
  author  = {Jumper, John and Evans, Richard and Pritzel, Alexander and Green, Tim and Figurnov, Michael and Ronneberger, Olaf and Tunyasuvunakool, Kathryn and Bates, Russ and {\v{Z}}{\'\i}dek, Augustin and Potapenko, Anna and Bridgland, Alex and Meyer, Clemens and Kohl, Simon A A and Ballard, Andrew J and Cowie, Andrew and Romera-Paredes, Bernardino and Nikolov, Stanislav and Jain, Rishub and Adler, Jonas and Back, Trevor and Petersen, Stig and Reiman, David and Clancy, Ellen and Zielinski, Michal and Steinegger, Martin and Pacholska, Michalina and Berghammer, Tamas and Bodenstein, Sebastian and Silver, David and Vinyals, Oriol and Senior, Andrew W and Kavukcuoglu, Koray and Kohli, Pushmeet and Hassabis, Demis},
  journal = {Nature},
  title   = {Highly accurate protein structure prediction with {AlphaFold}},
  year    = {2021},
  volume  = {596},
  number  = {7873},
  pages   = {583--589},
  doi     = {10.1038/s41586-021-03819-2}
}
```

此外，如果您使用了 AlphaFold-Multimer 模式，请同时引用：

```bibtex
@article {AlphaFold-Multimer2021,
  author       = {Evans, Richard and O{\textquoteright}Neill, Michael and Pritzel, Alexander and Antropova, Natasha and Senior, Andrew and Green, Tim and {\v{Z}}{\'\i}dek, Augustin and Bates, Russ and Blackwell, Sam and Yim, Jason and Ronneberger, Olaf and Bodenstein, Sebastian and Zielinski, Michal and Bridgland, Alex and Potapenko, Anna and Cowie, Andrew and Tunyasuvunakool, Kathryn and Jain, Rishub and Clancy, Ellen and Kohli, Pushmeet and Jumper, John and Hassabis, Demis},
  journal      = {bioRxiv},
  title        = {Protein complex prediction with AlphaFold-Multimer},
  year         = {2021},
  elocation-id = {2021.10.04.463034},
  doi          = {10.1101/2021.10.04.463034},
  URL          = {https://www.biorxiv.org/content/early/2021/10/04/2021.10.04.463034},
  eprint       = {https://www.biorxiv.org/content/early/2021/10/04/2021.10.04.463034.full.pdf},
}
```

## 社区贡献

由社区提供的 Colab Notebook（请注意，这些 Notebook 可能与我们完整的 AlphaFold 系统存在差异，我们并未验证其准确性）：

*   由 Martin Steinegger、Sergey Ovchinnikov 与 Milot Mirdita 开发的 [ColabFold AlphaFold2 notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)，使用托管于 Södinglab 基于 MMseqs2 服务器的 API [(Mirdita 等人 2019, Bioinformatics)](https://academic.oup.com/bioinformatics/article/35/16/2856/5280135) 创建多序列比对。

## 致谢

AlphaFold 调用或引用了以下独立的外部程序库与软件包：

*   [Abseil](https://github.com/abseil/abseil-py)
*   [Biopython](https://biopython.org)
*   [Colab](https://research.google.com/colaboratory/)
*   [Docker](https://www.docker.com)
*   [HH Suite](https://github.com/soedinglab/hh-suite)
*   [HMMER Suite](http://eddylab.org/software/hmmer)
*   [Haiku](https://github.com/deepmind/dm-haiku)
*   [JAX](https://github.com/google/jax/)
*   [Kalign](https://msa.sbc.su.se/cgi-bin/msa.cgi)
*   [matplotlib](https://matplotlib.org/)
*   [ML Collections](https://github.com/google/ml_collections)
*   [NumPy](https://numpy.org)
*   [OpenMM](https://github.com/openmm/openmm)
*   [OpenStructure](https://openstructure.org)
*   [pymol3d](https://github.com/avirshup/py3dmol)
*   [Sonnet](https://github.com/deepmind/sonnet)
*   [TensorFlow](https://github.com/tensorflow/tensorflow)
*   [Tree](https://github.com/deepmind/tree)
*   [tqdm](https://github.com/tqdm/tqdm)

我们衷心感谢所有上述项目的贡献者与维护者！

## 联系我们

如果您有任何本文档未涵盖的问题，请通过 [alphafold@deepmind.com](mailto:alphafold@deepmind.com) 联系 AlphaFold 团队。

我们非常期待收到您的使用反馈，并深入了解 AlphaFold 是如何在您的科研工作中发挥作用的。欢迎通过 [alphafold@deepmind.com](mailto:alphafold@deepmind.com) 与我们分享您的研究成果与使用故事。

## 许可证与免责声明

这不是 Google 官方支持的产品。

Copyright 2022 DeepMind Technologies Limited.

AlphaFold 2 及其预测结果仅用于理论建模研究。它们未经过临床使用验证或批准。您不应将 AlphaFold 2 或其输出用于临床目的，或将其作为医疗及其他专业建议的依据。与这些主题相关的任何内容仅供参考，不能替代合格专业人士的诊断或建议。

AlphaFold 2 的输出为具有不同置信度水平的预测结果，应审慎进行生物学解释。在依赖、发表、下载或以其他方式使用 AlphaFold 2 及其输出时，请审慎行事。

### AlphaFold 代码许可证

根据 Apache License 2.0 版本（“许可证”）获得许可；除非遵守许可证，否则您不得使用此文件。您可以在 https://www.apache.org/licenses/LICENSE-2.0 获取许可证副本。

除非适用法律要求或书面同意，根据许可证分发的软件均按“原样”分发，不附带任何明示或暗示的担保或条件。有关许可证管辖的特定语言权限和限制，请参阅许可证。

### 模型参数许可证

AlphaFold 模型参数基于知识共享署名 4.0 国际许可证（Creative Commons Attribution 4.0 International, CC BY 4.0）提供。您可以在以下网址找到详细条款：https://creativecommons.org/licenses/by/4.0/legalcode

### 第三方软件

使用上述[致谢](#致谢)部分中提及的第三方软件、库或代码可能受单独的条款和条件或许可证规定的约束。您对第三方软件、库或代码的使用须遵守任何此类条款，在使用前您应检查自己是否符合任何适用的限制或条款条件。

### 镜像数据库

DeepMind 对以下数据库建立了镜像，供参考使用：

*   [BFD](https://bfd.mmseqs.com/)（未修改版本），由 Steinegger M. 和 Söding J. 提供，依据 [知识共享署名-相同方式共享 4.0 国际许可证 (CC BY-SA 4.0)](http://creativecommons.org/licenses/by-sa/4.0/) 提供。
*   [BFD](https://bfd.mmseqs.com/)（修改版本），由 Steinegger M. 和 Söding J. 创建，经 DeepMind 修改，依据 [知识共享署名-相同方式共享 4.0 国际许可证 (CC BY-SA 4.0)](http://creativecommons.org/licenses/by-sa/4.0/) 提供。详见 [AlphaFold 蛋白质组论文](https://www.nature.com/articles/s41586-021-03828-1) 的“方法（Methods）”章节。
*   [Uniref30: v2021_03](http://wwwuser.gwdg.de/~compbiol/uniclust/2021_03/)（未修改版本），由 Mirdita M. 等人提供，依据 [知识共享署名-相同方式共享 4.0 国际许可证 (CC BY-SA 4.0)](http://creativecommons.org/licenses/by-sa/4.0/) 提供。
*   [MGnify: v2022_05](http://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2022_05/README.txt)（未修改版本），由 Mitchell AL 等人提供，免除所有版权限制，在 [CC0 1.0 通用公有领域贡献 (CC0 1.0 Universal)](https://creativecommons.org/publicdomain/zero/1.0/) 下完全免费供非商业和商业用途使用。

---

> 💡 **文档维护说明**：本中文文档由社区志愿者（@JasonYeYuhe）翻译维护，最后同步更新于 2026年10月4日。如发现内容与官方英文原版存在差异或新特性滞后，欢迎提交 PR 共同完善！
