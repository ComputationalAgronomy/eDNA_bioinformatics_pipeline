import shlex
import subprocess

from stage.runner import Runner


class SubprocessRunner(Runner):
    def __init__(self, prog_name, command, config, shell=False):
        super().__init__(prog_name, config)

        self.command = command
        self.shell = shell

    def run(self):

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
        if not self.dry:
            try:
                with open(self.outfile, "w") as f:
                    output = self.runner.capture_output.stdout
                    f.write(output)
                    return True
            except AttributeError as e:
                self.logger.write(f"{Runner.MSG_LOG} FAIL: {self.prog_name}. Error: {e}.\n")
            return False
