import os
import re
import sys

from argparse import ArgumentParser, RawTextHelpFormatter
from functools import partial
from textwrap import dedent, indent


def _progname(f):
    return os.path.basename(__file__)


def _dedent(msg):
    return re.sub(r'^\s*', '', msg, flags=re.M)


def usage(f, msg):
    p = _progname(f)

    offs = len("usage: ") + len(p) + 1

    return "{} [-h]\n".format(p) + indent(_dedent(msg), ' ' * offs)


def progress(msg, i, num):
    print((msg + ": {} / {}".format(i + 1, num)).ljust(80) + "\r",
          end='',
          file=sys.stderr,
          flush=True)


def progress_done():
    print(file=sys.stderr, flush=True)


def error(f, msg):
    p = _progname(f)

    print("{}: error: {}".format(p, msg), file=sys.stderr)

    sys.exit(1)


class NiceArgumentParser(ArgumentParser):
    def __init__(self, max_help_position=80, **kwargs):
        def formatter_class(prog):
            return RawTextHelpFormatter(
                prog, max_help_position=max_help_position)

        super().__init__(formatter_class=formatter_class, **kwargs)
