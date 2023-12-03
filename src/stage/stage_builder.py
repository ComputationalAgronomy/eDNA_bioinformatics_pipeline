from abc import ABC, abstractmethod

from stage.runner import Runner
from stage.stage_config import StageConfig
from stage.subproces_runner import RedirectOutputRunner, SubprocessRunner


class StageBuilder(ABC):
    """
    config : StageConfig
    runners : list of Runner
    output : list out stdout (or others) return from Runner
    """

    def __init__(self, heading, config: StageConfig):
        self.heading = heading
        self.config = config
        self.runners = []
        self.output = []

    def add_stage(self, prog_name, command, shell=False):
        stage = SubprocessRunner(prog_name, command, self.config, shell=shell)
        self.runners.append(stage)

    def add_stage_output_to_file(self, prog_name, stage, outfile_name):
        rd_stage = RedirectOutputRunner(prog_name, self.runners[stage], outfile_name, self.config)
        self.runners.append(rd_stage)

    def summary(self):
        messages = []
        for i, runner in enumerate(self.runners):
            messages.append(f"Step {i}: {runner.message}")
        return messages

    def run(self):
        self.config.logger.write(f"{Runner.MSG_LOG} Running: {self.heading}\n")
        self.output = []
        for i, runner in enumerate(self.runners):
            out = runner.run()
            self.output.append(out)
        self.config.logger.flush()
        return all(self.output)

    @abstractmethod
    def setup(self):
        pass
