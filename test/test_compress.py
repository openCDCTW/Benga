import tarfile

def to_gz(target, filename):
    tar = tarfile.open(filename + ".tar.gz", "w:gz")
    tar.add(target, arcname=filename)
    tar.close()

folder = "/media/pika/Workbench/workspace/pywgMLST/test_database"
compressed_file = "db"
to_gz(folder, compressed_file)

