# amplicon_qiime2
qiime2 pipline



## 方法

- 特征序列分析

  ```shell
  python Amplicon_run.py -n SRR18505770 -i SRR18505770_1.fastq SRR18505770_2.fastq -o $PWD/SRR18505770
  ```
  
  

- 多样本分析

  ```shell
  python ../Amplicon_run.py -type merge -m metadata.csv -o outdir -i $PWD/SRR18505770 $PWD/SRR18505774 $PWD/SRR18505775 $PWD/SRR18505775-1 $PWD/SRR18505775_trime $PWD/SRR18505786 $PWD/SRR18505800 $PWD/SRR18505803 $PWD/SRR18505854 $PWD/SRR18505864 $PWD/SRR18505865 $PWD/SRR18505869
  
  # $PWD/SRR18505770 为运行过特征序列分析的结果目录
  ```

  
  
  - metadata.csv格式如下（`group`与`11`对应组名）
  
  ```
  sampleid        group   11
  SRR18505770     1_1     测试1
  SRR18505774     1_1     1_2
  SRR18505775     1_1     1_2
  SRR18505775-1   1_1     1_2
  SRR18505775_trime       1_1     1_2
  SRR18505786     1_2     1_2
  SRR18505800     1_2     1_3
  SRR18505803     1_2     1_3
  SRR18505854     1_2     1_3
  SRR18505864     1_2     1_3
  SRR18505865     1_2     1_3
  SRR18505869     1_3     1_3
  ```
  
  
  
  

## 结果

- 特征序列分析

  - status_report

    - Fastqc
    - Asv

  - Fastqc.json
  - ASV.json

  | 参数                              | 类型  | 意义                   | 说明        |
  | --------------------------------- | ----- | ---------------------- | ----------- |
  | **seq_stats**                     | dict  | 序列统计               |             |
  | input                             | int   | 输入序列数量           |             |
  | filtered                          | int   | 质量过滤后序列数量     |             |
  | percentage_of_input_passed_filter | float | 质量过滤后序列百分比   | 百分制      |
  | denoised                          | int   | 噪音变异过滤后序列数量 |             |
  | non_chimeric                      | int   | 嵌合过滤后序列数量     |             |
  | percentage_of_input_non_chimeric  | float | 嵌合过滤后序列百分比   | 百分制      |
  |                                   |       |                        |             |
  | **asv_fre_sta**                   | dict  | 特征序列统计           |             |
  | asv_frequency_count               | int   | 特征序列数量           |             |
  | asv_frequency_mean                | float | 频数平均值             |             |
  | asv_frequency_std                 | float | 频数标准差             |             |
  | asv_frequency_min                 | int   | 频数最小值             |             |
  | asv_frequency_max                 | int   | 频数最大值             |             |
  | **asv_len_sta**                   | dict  |                        |             |
  | asv_seq_len_mean                  | float | 长度平均值             |             |
  | asv_seq_len_std                   | float | 长度标准差             |             |
  | asv_seq_len_min                   | int   | 长度最小值             |             |
  | asv_seq_len_max                   | int   | 长度最大值             |             |
  |                                   |       |                        |             |
  | asv_tax_csv                       | str   | 特征序列详情           | csv文件路径 |
  | `无标题`                          | str   | 特征序列ID             |             |
  | frequency                         | int   | 序列频数               |             |
  | Sequence                          | str   | 序列                   |             |
  | Taxon                             | str   | 物种分类               |             |
  | Confidence                        | float | 物种分类可信度         |             |

  

- 多样本分析

  - Merge.json
  
  | 参数                  | 类型  | 意义                                           | 说明                                                         |
  | --------------------- | ----- | ---------------------------------------------- | ------------------------------------------------------------ |
  | **merged_tab_csv**    | dict  | 物种丰度合并表                                 |                                                              |
  | `无标题`              | str   | 特征序列ID                                     |                                                              |
  | Sequence              | str   | 序列                                           |                                                              |
  | Taxon                 | str   | 物种分类                                       |                                                              |
  | Confidence            | float | 物种分类可信度                                 |                                                              |
  | `样本编号`            | int   | 序列频数                                       |                                                              |
  |                       |       |                                                |                                                              |
  | **alpha_rarefaction** | list  | 样本内多样性（Alpha多样性）                    | 列标题depth-X_iterY表示样本序列数量为X的第Y次重复结果（即：depth-1_iter10代表示样本序列数量为1的第10次重复结果）。每个样本序列数量有10次重复，共10种样本数量。每次任务的样本序列数量以1起始，但之后不固定的，但是总共为10种样本数量。 |
  | faith_pd              | str   | faith_pd（Alpha多样性指数）的结果文件          |                                                              |
  | observed_features     | str   | observed_features（Alpha多样性指数）的结果文件 |                                                              |
  | shannon               | str   | shannon（Alpha多样性指数）的结果文件           |                                                              |
  | **beta_diversity**    | list  | 样本间多样性（Beta多样性）                     | 四种分析方法（jaccard, bray_curtis, unweighted_unifrac, weighted_unifrac） |
  | sample_csv            | str   | PCoA分析                                       | 列标题为坐标轴（结果展示为前3个坐标轴），行标题为样本，值为样本在坐标轴上的位置 |
  | exp_csv               | str   | PCoA分析                                       | 各坐标轴的可解释方差，即括号内的百分比（未百分制）           |
  | distance_matrix_csv   | str   | 样本距离分析                                   | 行标题与列标题分别为样本，值为样本与样本间的距离关系值       |
  |                       |       |                                                |                                                              |
  | **alpha_res**         | list  | Alpha多样性差异                                |                                                              |
  | group                 |       | 分组方式                                       |                                                              |
  | metrics               | list  |                                                | 只选择metrics为faith_pd的                                    |
  | info                  | list  | 分组                                           |                                                              |
  | data                  |       | 分组的值                                       |                                                              |
  | H                     | float | H全部                                          |                                                              |
  | p                     | float | p全部                                          |                                                              |
  | pairwise              | list  | 各分组间比较的H, p-value和q-value              |                                                              |
  |                       |       |                                                |                                                              |
  | **ancom_res**         | list  | 物种差异                                       |                                                              |
  | group                 | str   | 分组方式                                       |                                                              |
  | data                  | list  | 差异物种分析                                   | ANCOM差异异分析火山图中所有物种分类，W和clr值                |
  | ancom                 | list  | 差异物种分析（标色的点）                       | ANCOM差异异分析火山图中，有差异的的物种分类                  |
  | abund_percent         | list  | 差异物种丰度                                   | 差异物种丰度箱线图的标记                                     |
  | abund_group           | list  | 差异物种丰度                                   | 差异物种丰度箱线图的分组                                     |
  | abund                 | dict  | 差异物种丰度                                   | 差异物种丰度箱线图,每个物种的值（按照abund_percent及abund_group顺序） |
  |                       |       |                                                |                                                              |

 
