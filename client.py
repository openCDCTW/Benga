import sys
from src.utils import transport


def main(folder, jobid):
    sender = transport.Sender("127.0.0.1", 5556)
    sender.send(folder, jobid)


if __name__ == "__main__":
    f, j = sys.argv[1:3]
    main(f, j)
