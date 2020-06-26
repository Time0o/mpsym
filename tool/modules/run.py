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


def normalize_gap_output(out):
    out = out.replace(' ', '')

    out_list = [c for c in out]

    i = 0
    while i < len(out_list):
        if out_list[i] == '\n' and out_list[i - 1] in ':),':
            out_list[i] = ''
        i += 1

    return ''.join(out_list)
