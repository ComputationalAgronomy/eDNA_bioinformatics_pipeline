import shlex
import subprocess

from stage.runner import Runner


class SubprocessRunner(Runner):

    """
    Runner class for executing a command as a subprocess.

    Attributes:
        command (str): The command to execute.
        shell (bool): Whether to execute the command through the shell.
    """

    def __init__(self, prog_name: str, command: str, config, shell: bool = False):
        super().__init__(prog_name, config)

        self.command = command
        self.shell = shell

    def run(self) -> bool:
        """
        Executes the command as a subprocess.

        :returns: True if the execution is successful, False otherwise.
        """
        if self.verbose:
            print(self.message)

        if self.shell:
            args = self.command
        else:
            args = shlex.split(self.command)

        if not self.dry:
            self.logger.write(self.message)
            try:
                self.capture_output = subprocess.run(args,
                                                   capture_output=True,
                                                   shell=self.shell,
                                                   check=True,
                                                   text=True)
                self.logger.write(f"{Runner.MSG_LOG} COMPLETE: {self.prog_name}.\n")
                return True
            except subprocess.CalledProcessError as e:
                self.logger.write(f"{Runner.MSG_LOG} FAIL: {self.prog_name}. SubprocessError: {e}.\n")
                return False
            except FileNotFoundError as e:
                self.logger.write(f"{Runner.MSG_LOG} FAIL: {self.prog_name}. FileNotFoundError: {e}.\n")
                return False
            except Exception as e:
                self.logger.write(f"{Runner.MSG_LOG} FAIL: {self.prog_name}. Other Exception: {e}.\n")
                return False


class RedirectOutputRunner(Runner):
    """
    Runner class for redirecting the output of a subprocess to a file.

    Attributes:
        runner (SubprocessRunner): The subprocess runner whose output to redirect.
        outfile (str): The file to write the output to.
    """
    def __init__(self, prog_name: str, runner: SubprocessRunner, outfile: str, config):
        super().__init__(prog_name, config)

        if isinstance(runner, SubprocessRunner):
            self.runner = runner
            self.outfile = outfile
            self.message = f"{Runner.MSG_LOG} RedirectOutput: {self.prog_name}."
        else:
            info = type(runner)
            raise TypeError(f"Invalid instance type: {info}.")

    def run(self) -> bool:
        """
        Redirects the output of the runner to the specified file.

        :returns: True if the output redirection is successful, False otherwise.
        """
        if not self.dry:
            try:
                with open(self.outfile, "w") as f:
                    output = self.runner.capture_output.stdout
                    f.write(output)
                    return True
            except AttributeError as e:
                self.logger.write(f"{Runner.MSG_LOG} FAIL: {self.prog_name}. Error: {e}.\n")
            return False
