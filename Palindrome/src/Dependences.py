import os
from Config import ScrollOfEternalTruths

class ArsenalSummoner:
    """Weapons master responsible for equipping external bioinformatics blades like BLAST and DIAMOND."""
    
    def __init__(self):
        self.scroll = ScrollOfEternalTruths()

    def inspect_blade_durability(self, forged_path):
        """Verifies if the weapon (executable) truly exists and can be wielded."""
        if not forged_path:
            return False
        return os.path.isfile(forged_path) and os.access(forged_path, os.X_OK)

    def tame_blast_beast(self, new_sword_path):
        """Saves the path of the BLAST executable in the grimoire."""
        self.scroll.ancestral_memories["blast_path"] = new_sword_path
        self.scroll.carve_runes_into_stone()

    def find_lost_diamond(self, new_staff_path):
        """Saves the path of the DIAMOND executable in the grimoire."""
        self.scroll.ancestral_memories["diamond_path"] = new_staff_path
        self.scroll.carve_runes_into_stone()

    def reveal_current_armory(self):
        """Returns the saved executable paths."""
        return {
            "blast": self.scroll.ancestral_memories.get("blast_path", ""),
            "diamond": self.scroll.ancestral_memories.get("diamond_path", "")
        }
