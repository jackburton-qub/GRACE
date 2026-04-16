"""
blast_worker.py — QThread worker for BLAST specificity checks.
BLAST binaries are resolved automatically from the bundled blast/ folder.
"""
from PyQt6.QtCore import QThread, pyqtSignal


class BLASTWorker(QThread):
    progress = pyqtSignal(int, int)   # done, total
    finished = pyqtSignal(list)       # results list
    error    = pyqtSignal(str)        # error message

    def __init__(self, genome_path, primers, blast_params, specificity_params, parent=None):
        super().__init__(parent)
        self.genome_path        = genome_path
        self.primers            = primers
        self.blast_params       = blast_params
        self.specificity_params = specificity_params

    def run(self):
        try:
            from core.primer_specificity_blast import check_specificity_blast
            results = check_specificity_blast(
                genome_fasta=self.genome_path,
                primer_results=self.primers,
                blast_params=self.blast_params,
                specificity_params=self.specificity_params,
                # blast_bin_dir omitted — resolved automatically from bundled blast/
            )
            self.finished.emit(results)
        except Exception as e:
            self.error.emit(str(e))