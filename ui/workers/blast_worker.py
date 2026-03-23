"""
blast_worker.py — QThread worker stub
"""
from PyQt6.QtCore import QThread, pyqtSignal


class BLASTWorker(QThread):
    progress = pyqtSignal(int, int)   # done, total
    finished = pyqtSignal(object)     # result
    error    = pyqtSignal(str)        # error message

    def __init__(self, parent=None):
        super().__init__(parent)

    def run(self):
        pass  # implement in full version
