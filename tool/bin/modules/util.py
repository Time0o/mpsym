import os
import re
import sys

from functools import partial
from textwrap import dedent, indent


def _progname(f):
    return os.path.basename(f)


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
