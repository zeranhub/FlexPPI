# Direct panel-f entry with the complete complex represented.
# File > Run Script; all paths are derived from this file.
python
from pathlib import Path
from pymol import cmd
_f_bundle = Path(__script__).resolve().parent
cmd.load(str(_f_bundle / "ubiquitin_2QHO.pse"))
cmd.scene("02_Q62_MD_observed", "recall", animate=0)
python end
