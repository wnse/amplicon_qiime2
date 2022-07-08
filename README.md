# amplicon_qiime2
qiime2 pipline



## 方法

- 特征序列分析

  ```shell
  python /Bio/Amplicon/Amplicon_run.py \
  	-n SRR18505770 \
  	-i SRR18505770_1.fastq SRR18505770_2.fastq \
  	-o $PWD/SRR18505770
  	-tID testamplicon
  	
  # -n 为样本名称
  # -i 为样本的输入数据
  # -o 为结果目录
  # -tID taskID
  ```

  

- 多样本分析

  ```shell
  python /Bio/Amplicon/Amplicon_run.py \
  	-type merge \
  	-tag DiversityAmplicon \
  	-i $PWD/SRR18505770 $PWD/SRR18505774 $PWD/SRR18505775 $PWD/SRR18505775_trime $PWD/SRR18505786 $PWD/SRR18505800 $PWD/SRR18505803 $PWD/SRR18505854 $PWD/SRR18505864 $PWD/SRR18505865 $PWD/SRR18505869 $PWD/SRR18505869_1 \
  	-m $PWD/metadata.csv \
  	-o $PWD/test_docker_multi \
  	-tID testamplicon
  
  # -type 固定为merge，为了与特征序列分析进行区分
  # -tag 为应用名称，生成结果中的json文件及status_report.txt中节点名称
  # -i 为各样本运行过特征序列分析的结果目录
  # -m 样本信息分组csv表，具体格式见下
  # -o 结果目录
  # -tID taskID
  ```
  
  
  
  - metadata.csv格式如下（`group`与`sampletype`分类组名）
  
    **注意：metadata.csv中的列标题不能为中文，首列为样本名称，标题必须为sampleid**
  
  ```
  sampleid	group	sampletype
  SRR18505770	G1	血液
  SRR18505774	G1	血液
  SRR18505775	G1	血液
  SRR18505775_trime	G1	血液
  SRR18505786	G1	血液
  SRR18505800	G2	血液
  SRR18505803	G3	血液
  SRR18505854	G3	咽拭子
  SRR18505864	G3	咽拭子
  SRR18505865	G3	咽拭子
  SRR18505869	G3	咽拭子
  SRR18505869_1	G3	咽拭子
  ```
  
  
  
  

## 结果

- 特征序列分析

  - status_report

    - FastqcAmplicon
    - AsvAmplicon

  - FastqcAmplicon.json
  - AsvAmplicon.json

  | 参数                              | 类型  | 意义                   | 说明         |
  | --------------------------------- | ----- | ---------------------- | ------------ |
  | **seq_stats**                     | dict  | 序列统计               |              |
  | input                             | int   | 输入序列数量           |              |
  | filtered                          | int   | 质量过滤后序列数量     |              |
  | percentage_of_input_passed_filter | float | 质量过滤后序列百分比   | 百分制       |
  | denoised                          | int   | 噪音变异过滤后序列数量 |              |
  | non_chimeric                      | int   | 嵌合过滤后序列数量     |              |
  | percentage_of_input_non_chimeric  | float | 嵌合过滤后序列百分比   | 百分制       |
  |                                   |       |                        |              |
  | **asv_info**                      | list  | 特征序列详情           |              |
  | sample_id                         | str   | 特征序列ID             |              |
  | frequency                         | str   | 序列频数               |              |
  | Sequence                          | str   | 序列                   |              |
  | Taxon                             | str   | 物种分类               |              |
  | Confidence                        | float | 物种分类可信度         | 需转百分制   |
  |                                   |       |                        |              |
  | **asv_fre_sta**                   | dict  | 特征序列统计           |              |
  | asv_frequency_count               | int   | 特征序列数量           |              |
  | asv_frequency_mean                | float | 频数平均值             |              |
  | asv_frequency_std                 | float | 频数标准差             |              |
  | asv_frequency_min                 | int   | 频数最小值             |              |
  | asv_frequency_max                 | int   | 频数最大值             |              |
  |                                   |       |                        |              |
  | **asv_len_sta**                   | dict  |                        |              |
  | asv_seq_len_mean                  | float | 长度平均值             |              |
  | asv_seq_len_std                   | float | 长度标准差             |              |
  | asv_seq_len_min                   | int   | 长度最小值             |              |
  | asv_seq_len_max                   | int   | 长度最大值             |              |
  |                                   |       |                        |              |
  | **tax_krona_html**                | html  | 物种分类统计           | html文件路径 |
  | asv_tax_csv                       | str   | 特征序列详情           | csv文件路径  |
  | `无标题`                          | str   | 特征序列ID             |              |
  | frequency                         | int   | 序列频数               |              |
  | Sequence                          | str   | 序列                   |              |
  | Taxon                             | str   | 物种分类               |              |
  | Confidence                        | float | 物种分类可信度         |              |

  

- 多样本分析

  - DiversityAmplicon.json
  
    *文件名称随运行脚本时参数`-tag`改变*
  
  | 参数                  | 类型  | 意义                              | 说明                                                         |
  | --------------------- | ----- | --------------------------------- | ------------------------------------------------------------ |
  | merged_tab_csv        | dict  | 物种丰度合并表                    |                                                              |
  | `无标题`              | str   | 特征序列ID                        |                                                              |
  | Sequence              | str   | 序列                              |                                                              |
  | Taxon                 | str   | 物种分类                          |                                                              |
  | Confidence            | float | 物种分类可信度                    |                                                              |
  | `样本编号`            | int   | 序列频数                          |                                                              |
  |                       |       |                                   |                                                              |
  | **taxa_info**         | list  | 物种丰度分布                      |                                                              |
  | level                 | int   | 分类级别                          |                                                              |
  | legend                | str   | 图例                              |                                                              |
  | sortby                | str   | 排序方式                          |                                                              |
  | info                  | list  | 各样本各物种的序列数量            |                                                              |
  |                       |       |                                   |                                                              |
  | **alpha_rarefaction** | list  | 样本内多样性（Alpha多样性）       | 三种分析方法（shannon, faith_pd, observed_features）         |
  | group                 | str   | 分组方式                          |                                                              |
  | *metrics*             | list  | 各分析方法的详情                  |                                                              |
  | metric                | str   | 分析方法名称                      |                                                              |
  | info                  | dict  | 分析结果                          |                                                              |
  | columns               | list  | data中各元素意义                  | 【分组信息，横坐标，柱子高度最小值-最大值，样本数量】        |
  | data                  | list  | 对应坐标的值                      |                                                              |
  |                       |       |                                   |                                                              |
  | **beta_diversity**    | list  | 样本间多样性（Beta多样性）        | 四种分析方法（jaccard, bray_curtis, unweighted_unifrac, weighted_unifrac） |
  |                       |       |                                   |                                                              |
  | *sample_pcoa*         | list  | PCoA分析                          |                                                              |
  | sample                | str   | 样本名称                          |                                                              |
  | info                  | dict  | 各轴的位置                        | 结果只能展示前3个坐标轴                                      |
  | pcoa_exp              | dict  | PCoA分析                          | 各坐标轴的可解释方差，即坐标轴括号内的百分比（未百分制）     |
  | distance_matrix       | dict  | 样本距离分析                      | 每个样本与其他各样本之间的距离值                             |
  |                       |       |                                   |                                                              |
  | **alpha_res**         | list  | Alpha多样性差异                   | 四种分析方法：faith_pd，observed_features，pielou_evenness，shannon_entropy |
  | group                 |       | 分组信息                          |                                                              |
  | *metrics*             | list  | 各方法的差异分析结果              |                                                              |
  | metric                | str   | 分析方法                          |                                                              |
  | info                  | list  | 差异结果                          |                                                              |
  | index                 | list  | 各组横坐标                        |                                                              |
  | data                  | list  | 各组具体值                        |                                                              |
  | H                     | float | H全部                             |                                                              |
  | p                     | float | p全部                             |                                                              |
  | pairwise              | list  | 各分组间比较的H, p-value和q-value |                                                              |
  |                       |       |                                   |                                                              |
  | **ancom_res**         | list  | 物种差异                          |                                                              |
  | group                 | str   | 分组信息                          |                                                              |
  | data                  | list  | 差异物种分析                      | ANCOM差异异分析火山图中所有物种分类，W和clr值                |
  | ancom                 | list  | 差异物种分析（标色的点）          | ANCOM差异异分析火山图中，有差异的的物种分类                  |
  | abund_percent         | list  | 差异物种丰度                      | 差异物种丰度箱线图的标记                                     |
  | abund_group           | list  | 差异物种丰度                      | 差异物种丰度箱线图的分组                                     |
  | abund                 | dict  | 差异物种丰度                      | 差异物种丰度箱线图,每个物种的值（按照abund_percent及abund_group顺序） |
  |                       |       |                                   |                                                              |

 
