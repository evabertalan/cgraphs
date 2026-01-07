import sys
from .view.view import start


def main(args):
    if len(args) < 1 or args[0] != "start":
        raise ValueError(
            "Invalid command. To start the tool use: python3 -m cgraphs start"
        )
    start()


if __name__ == "__main__":
    main(sys.argv[1:])
