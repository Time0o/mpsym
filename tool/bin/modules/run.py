import os

from subprocess import Popen, DEVNULL, PIPE
from tempfile import NamedTemporaryFile


def run_gap(script):
    with NamedTemporaryFile(mode='w+', delete=False) as f:
        f.write(script)

    p = Popen(['gap', '--nointeract', '-q', f.name], stdout=PIPE, stderr=DEVNULL)
    out, _ = p.communicate()
    out = out.decode('ascii')[:-1]

    os.remove(f.name)

    return out
