
import os
import subprocess
from abc import ABC, abstractmethod

from stage.stage_config import StageConfig


class Runner(ABC):

    MSG_LOG = "==LOG=="
    MSG_DEBUG = "===DEBUG==="

    def __init__(self, prog_name, config: StageConfig):
        self.prog_name = prog_name  # for debug/naming only
        self.message = f"{Runner.MSG_LOG} Program: {self.prog_name}."
        self.capture_output = None
        try:
            self.verbose = config.verbose
            self.dry = config.dry
            self.logger = config.logger
        except AttributeError as e:
            print(e)


    @abstractmethod
    def run(self):
        """
        Return True/False for success/failure.
        """
        pass
