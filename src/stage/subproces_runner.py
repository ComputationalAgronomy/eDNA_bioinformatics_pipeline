import shlex
import subprocess

from stage.runner import Runner




class SubprocessRunner(Runner):
    def __init__(self, prog_name, command, config, shell=False):
        super().__init__(prog_name, config)

        self.command = command
        self.shell = shell
        self.comp_process = None


    def run(self):
        """Execute the command using subprocess

        `check_output and/or check_call` are for legacy use (python 2).
        Refactor it out to python 3 later.

        Return
        ------
        True or False:

        TODO(SW): Add features to fetch detail output.
        """
        if self.shell:
            args = self.command
        else:
            args = shlex.split(self.command)
        if self.verbose:
            print(self.message)
        if not self.dry:
            self.logger.write(self.message)
            try:
                # self.out = subprocess.check_output(args)
                self.comp_process = subprocess.run(args, capture_output=True,
                                                   shell=self.shell,
                                                   check=True, text=True)
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

    def __init__(self, prog_name, runner, outfile, config):
        super().__init__(prog_name, config)
        if isinstance(runner, SubprocessRunner):
            self.runner = runner
            self.outfile = outfile
            self.message = f"{Runner.MSG_LOG} RedirectOutput: {self.prog_name}."
        else:
            info = type(runner)
            raise TypeError(f"Invalid instance type: {info}.")

    def run(self):
        try:
            with open(self.outfile, "w") as f:
                output = self.runner.comp_process.stdout
                f.write(output)
                return True
        except AttributeError as e:
            self.logger.write(f"{Runner.MSG_LOG} FAIL: {self.prog_name}. Error: {e}.\n")
        return False
