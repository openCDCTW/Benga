import logging
from utils import *




def main():
    root = "/media/pika/Workbench/workspace/pywgMLST"
    annotate_configs(joinpath(root, "upload_Assembly"), joinpath(root, "upload_Assembly"), joinpath(root, "output/db"))


if __name__ == "__main__":
    main()