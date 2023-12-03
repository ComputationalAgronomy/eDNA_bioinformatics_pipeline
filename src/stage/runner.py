
from abc import ABC, abstractmethod
import os
import subprocess


class Runner(ABC):

    MSG_LOG = "==LOG=="
    MSG_DEBUG = "===DEBUG==="

    def __init__(self, prog_name, config):
        self.prog_name = prog_name  # for debug/naming only
        self.message = f"{Runner.MSG_LOG} Program: {self.prog_name}."
        self.comp_process = None
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
