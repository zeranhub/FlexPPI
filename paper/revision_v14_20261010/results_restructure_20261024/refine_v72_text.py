"""Keep the source evidence while shortening repeated interpretation."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent

def read(name):
    return json.loads((ROOT / name).read_text(encoding='utf-8'))

def write(name, value):
    (ROOT / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

first = read('sections_31_32.json')
first['3.1']['paragraphs'][2] = first['3.1']['paragraphs'][2].replace(
    '完整模型的三项优势分别为0.046、0.034和0.014。',
    '按未舍入值计算，完整模型的三项优势分别为0.046、0.034和0.015。')
first['3.1']['paragraphs'][3] = first['3.1']['paragraphs'][3].replace(
    '相对BackFlip-2.1（assembly）的优势分别为0.068和0.158。',
    '按未舍入值计算，相对BackFlip-2.1（assembly）的优势分别为0.069和0.158。')
write('sections_31_32.json', first)

last = read('sections_35_discussion.json')
last['4']['paragraphs'] = [
    'FlexPPI的核心价值是在同一静态结构上同时定位骨架活动与侧链状态差异。RMSF组织主链活动区域的相对排序，Hχ区分相同化学种类在不同局部环境中的联合状态占据；共享几何表示使两个动态维度对应到相同残基位置。类型均值与静态特征对照、完整与链内编码比较，以及同类型和家族隔离评价，共同支持主链朝向、邻域排列和完整空间环境的预测价值，为单链与复合物提供统一的局部动力学表征。',
    '完整复合物环境使柔性预测能够解析界面两侧的接触组织。KIX活动环区与泛素同类型侧链的实例分别展示了骨架定位和取向差异；MEL5亲和成熟实例将HLA Q72的侧链状态集中与未替换TCR βGly54的主链波动降低联系到同一接触片区。起始结构中Q72酰胺两端与β53/β54骨架的接近，以及重复轨迹中接触和氢键几何占据的增加，共同支持该片区的局部稳定化。两类预测将受体改造所伴随的动态响应延伸到邻近未替换位点和对侧侧链，为理解已有结合实验对应的局部接触变化提供机制线索。',
    '侧链输出通过构象占据与经验投影自由能建立物理对应。Hχ概括固定χ角区间下多角联合状态的占据多样性，同类型控制使比较建立在相同侧链自由度上。83条记录的配对相关与预测分组，以及泛素Q41/Q2的逐角曲线，均支持侧链输出优先定位取向范围更广、角度占据更多样的位置。经验投影自由能F*=−ln(p̂/p̂max)以所分析轨迹中占据最高的角区为零，描述该投影坐标下的相对占据偏好[3,24,25]；F*≤1的角区覆盖与细箱熵分别刻画主占据范围和整体占据分布。这些对应将静态结构预测落实到可由轨迹解析的构象特征，可用于确定后续局部模拟、接触分析和功能实验关注的残基位置。'
]
last['5']['paragraphs'] = [
    'FlexPPI以共享主链几何表示联合预测逐残基RMSF与侧链联合状态熵，在同一位置提供骨架波动和侧链状态多样性两个互补读数。五折体系留出、外部方法比较与同类型评价支持其预测优势，家族隔离对照支持完整结构环境的贡献。KIX、泛素及MEL5实例将两类读数对应到活动骨架、侧链取向差异和跨界面接触网络；83条记录的构象占据验证进一步支持Hχ对实际角度覆盖与占据多样性的表征，使联合柔性预测成为从静态结构定位局部动态特征的工具。'
]
write('sections_35_discussion.json', last)
print(json.dumps({key: sum(map(len, last[key]['paragraphs'])) for key in ['4', '5']}, ensure_ascii=False))
