# 结构用于定位，NMR用于描述运动

四套文件采用公开PDB沉积的第一MODEL原始坐标，CaM原始A链及其对应肽链分别输出为A链和P链，保留对应CaM的四个原始Ca²⁺。未补齐任何缺失原子或残基，也未优化或重新采样坐标。CaM通过完整序列比对映射到四个BMRB记录共同的成熟CaM序列，即Ala1起始的148个残基，输出PDB采用这一成熟编号。所有原始作者编号都在CSV中保留；肽沿用原始作者编号。

CaMKKα使用1CKK（原始A/B链），eNOS使用1NIW（原始A/B链），nNOS使用2O60（原始A/B链），smMLCK使用1CDL（原始A/E链）。CaMKKα和nNOS的CaM没有缺失残基；1NIW缺少1、2、148；1CDL缺少1–4、147、148。强调的Leu18、Leu39、Leu105、Met109、Met144均存在，全部重原子完整，其中甲基CD1/CD2或CE完整。其他位置的缺失重原子逐个保存在审计JSON。

1NIW的9个原始Met位置是硒代甲硫氨酸（MSE），包含109与144。序列映射将MSE视作M的序列等价物，但输出文件保留原始MSE残基名、Se元素及原坐标，没有把Se改成S。因此该结构可作eNOS复合物空间示意，不代表NMR样品完全相同的化学构成；图注应明确这一点。1CDL的COMPND文本将肽错误描述为CaMKII，而DBREF指向MYLK_CHICK，其肽序列与smMLCK NMR结合肽匹配；采用DBREF及实际序列核对，不采用有冲突的COMPND名称。

图中应说明“结构作为匹配伙伴复合物的空间定位参考，实验S²来自独立沉积的NMR数据；结构示意不是NMR运动轨迹”。具体距离为单个静态结构中的最近重原子距离，仅帮助定位，不是动态机制或因果关系证据。nNOS与smMLCK的Leu39和Leu105均靠近肽（侧链至肽最近重原子距离依次为3.702/3.958 Å与3.508/3.689 Å），不能画成一组接触、另一组无接触。

cam_all_residue_mapping.csv每套结构包含148个成熟CaM位置，含缺失标记，共592行；cam_key_site_structure_context.csv包含20个关键位点及其CA坐标、完整性检查和最近肽原子。structure_mapping_audit.json包含原始及输出SHA256、比对序列、缺失段、修改残基、第一MODEL选择和钙离子审计。
