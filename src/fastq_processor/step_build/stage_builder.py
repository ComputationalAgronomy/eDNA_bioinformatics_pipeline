from abc import ABC, abstractmethod
import os

from fastq_processor.step_build import stage_config
from fastq_processor.step_build import subproces_runner


class StageBuilder(ABC):
    """
    Abstract base class for building and managing stages of runners.

    Attributes:
        heading (str): The heading for the stage builder.
        config (StageConfig): The configuration for the stage.
        runners (list of Runner): List of runners to execute.
        output (list): List of outputs from the executed runners.
    """

    def __init__(self, heading: str, config: stage_config.StageConfig):
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
        stage = subproces_runner.SubprocessRunner(prog_name, command, self.config, shell=shell)
        self.runners.append(stage)

    def add_stage_output_to_file(self, prog_name: str, stage: int, outfile_name: str, errfile_name: str):
        """
        Adds a stage that redirects output to a file.

        :param prog_name: The name of the program to run.
        :param stage: The index of the stage whose output to redirect.
        :param outfile_name: The name of the file to write the output to.
        """
        rd_stage = subproces_runner.RedirectOutputRunner(prog_name, self.runners[stage], outfile_name, errfile_name, self.config)
        self.runners.append(rd_stage)

    def check_infile(self):
        if not os.path.exists(self.infile):
            raise FileNotFoundError(f"{self.infile} not found")

    def check_savedir(self):
        os.makedirs(self.save_dir, exist_ok=True)

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
