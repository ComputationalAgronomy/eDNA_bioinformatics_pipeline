import setuptools

import site
import sys


if __name__ == "__main__":
    # https://github.com/pypa/pip/issues/7953 .
    # site.ENABLE_USER_SITE = "--user" in sys.argv[1:]
    site.ENABLE_USER_SITE = 1
    # sys.argv[1:] = ["develop", "--user"]
    setuptools.setup()
