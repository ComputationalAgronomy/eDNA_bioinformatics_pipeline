import os
import subprocess

from analysis_toolkit.utils.base_logger import logger
#TODO(SW): If you only use these functions in ml_tree, put these there.
def check_mltree_overwrite(save_dir: str, prefix: str) -> str:
    # TODO(SW): Eventually, replace with argparse.
    """
    Check if an MLTree run with the given prefix already exists in the specified directory.
    If it does, prompt the user to decide whether to redo the run or not.

    :param save_dir: Directory where the output files are saved.
    :param prefix: Prefix for the output file names.
    :return: User input choice ('-redo', '--redo-tree', '--undo', or 'stop').
    """
    ckp_path = os.path.join(save_dir, prefix + '.ckp.gz') # e.g. save/dir/SpA.ckp.gz
    if os.path.exists(ckp_path):
        print(
            f"> MLTree checkpoint fileCheckpoint ({ckp_path}) indicates that a previous run successfully finished  already exists.\n"
            "Use `-redo` option if you really want to redo the analysis and overwrite all output files.\n"
            "Use `--redo-tree` option if you want to restore ModelFinder and only redo tree search.\n"
            "Use `--undo` option if you want to continue previous run when changing/adding options.\n"
        )
        user_input_choice = ['-redo', '--redo-tree', '--undo', 'stop']
        help_msg = f"({'/'.join(user_input_choice)}): "
        while True:
            user_input = input(help_msg).strip().lower()
            if user_input in user_input_choice:
                return user_input
            print("> Invalid input.")
    else:
        pass

def run_iqtree2(
        seq_path:str,
        save_dir: str,
        prefix: str,
        model:str = None,
        bootstrap: int = None,
        threads: int = None
    ) -> None:
    """
    Run IQTREE2 command.
    """
    checkpoint = check_mltree_overwrite(save_dir, prefix)
    if checkpoint == 'stop':
        logger.info("> Stopping the run.")
        return

    model = model or 'TEST'
    prefix_path = os.path.join(save_dir, prefix)
    threads = threads or 'AUTO'

    cmd = [
        'iqtree2', '-m', model, '-s', seq_path, '--prefix', prefix_path, '-nt', threads
    ]
    if bootstrap:
        cmd.extend(["-b", str(bootstrap)])
    if checkpoint:
        cmd.append(checkpoint)

    logger.info(f"Running IQTREE2 command: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        logger.info(f"IQTREE2 finished. Output files saved in: {save_dir}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error occurred during IQTREE2 run: {e}")