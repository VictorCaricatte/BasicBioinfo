import json
import os

class ScrollOfEternalTruths:
    """Safeguards the profane secrets, paths, and configurations of the journey."""
    
    def __init__(self):
        self.path_of_runes = "secret_configurations.json"
        self.ancestral_memories = {
            "blast_path": "",
            "diamond_path": "",
            "database_cavern": "local_db_dungeon",
            "batch_forge_cavern": "batch_results_dungeon"
        }
        self.unearth_ancient_secrets()

    def unearth_ancient_secrets(self):
        """Reads the JSON scroll if it exists on the physical plane."""
        if os.path.exists(self.path_of_runes):
            try:
                with open(self.path_of_runes, 'r', encoding='utf-8') as grimoire:
                    data = json.load(grimoire)
                    self.ancestral_memories.update(data)
            except Exception:
                pass # If the grimoire is corrupted, accept the default memories

    def carve_runes_into_stone(self):
        """Seals the configurations in JSON format for eternity."""
        with open(self.path_of_runes, 'w', encoding='utf-8') as grimoire:
            json.dump(self.ancestral_memories, grimoire, indent=4)
