from abc import ABC, abstractmethod
import os

from fastq_processor.step_build.runner import Runner
from fastq_processor.step_build.stage_config import StageConfig
from fastq_processor.step_build.subproces_runner import RedirectOutputRunner, SubprocessRunner


class StageBuilder(ABC):
    """
    Abstract base class for building and managing stages of runners.

    Attributes:
        heading (str): The heading for the stage builder.
        config (StageConfig): The configuration for the stage.
        runners (list of Runner): List of runners to execute.
        output (list): List of outputs from the executed runners.
    """

    def __init__(self, heading: str, config: StageConfig):
        self.heading = heading
        self.config = config
        self.runners = []
        self.output = []

    def add_stage(self, prog_name: str, command: str, shell=False):
        """
        Adds a stage to the list of runners.

        :param prog_name: The name of the program to run.
        :param command: The command to execute.
        :param shell: Whether to execute the command through the shell.
        """
        stage = SubprocessRunner(prog_name, command, self.config, shell=shell)
        self.runners.append(stage)

    def add_stage_output_to_file(self, prog_name: str, stage: int, outfile_name: str, errfile_name: str):
        """
        Adds a stage that redirects output to a file.

        :param prog_name: The name of the program to run.
        :param stage: The index of the stage whose output to redirect.
        :param outfile_name: The name of the file to write the output to.
        """
        rd_stage = RedirectOutputRunner(prog_name, self.runners[stage], outfile_name, errfile_name, self.config)
        self.runners.append(rd_stage)

    def check_path(self, path: str):
        os.makedirs(self.save_dir, exist_ok=True)
        if not os.path.exists(path):
            raise FileNotFoundError(f"{path} not found")

    def summary(self) -> list[str]:
        """
        Provides a summary of the stages.

        :returns: A list of messages summarizing each stage.
        """
        messages = []
        for i, runner in enumerate(self.runners):
            messages.append(f"Step {i}: {runner.message}")
        return messages

    def run(self):
        """
        Executes all the stages.

        Returns:
            list: The output from each stage.
        """
        self.config.logger.info(f"Running: {self.heading}")
        self.output = []
        for i, runner in enumerate(self.runners):
            out = runner.run()
            self.output.append(out)
        # self.config.logger.flush()
        self.runners = []
        return all(self.output)
    
    @abstractmethod
    def setup(self):
        pass
