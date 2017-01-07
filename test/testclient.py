import zmq
import time
import tarfile


def to_gz(target, filename):
    with tarfile.open(filename + ".tar.gz", "w:gz") as tar:
        tar.add(target, arcname=filename)

context = zmq.Context()
socket = context.socket(zmq.PUB)
socket.bind("tcp://127.0.0.1:5556")
time.sleep(1)
print("Connecting to server...")

# compress folder
folder = "/media/pika/Workbench/workspace/pywgMLST/test_database"
compressed_file = "db"  # use jobid
print("Compressing files...")
# to_gz(folder, compressed_file)

# transfer compressed file
with open(compressed_file + ".tar.gz", "rb") as f:
    socket.send_multipart(f.readlines())
    print("Sending file...")
    time.sleep(1)
socket.close()
context.term()

# delete compressed file