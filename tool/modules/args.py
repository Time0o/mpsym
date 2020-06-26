from argparse import ArgumentParser, RawTextHelpFormatter


class NiceArgumentParser(ArgumentParser):
    def __init__(self, max_help_position=80, **kwargs):
        def formatter_class(prog):
            return RawTextHelpFormatter(
                prog, max_help_position=max_help_position)

        super().__init__(formatter_class=formatter_class, **kwargs)
