import sys


class StageConfig:
    """
    Configuration shared between all stages
    """
    def __init__(self, verbose=False, dry=False,
                 logger=sys.stdout,
                 n_cpu=1, memory=8):
        self.verbose = verbose
        self.dry = dry
        self.logger = logger
        self.n_cpu = n_cpu
        self.memory = memory

    def get_machine_info(self):
        return {"n_cpu": self.n_cpu,
                "memory": self.memory
                }

    def get_basic_configuration(self):
        return {"verbose": self.verbose,
                "dry_run": self.dry,
                "logger": self.logger
                }
