import zmq
import time
import tarfile
import os


def decompress(target, path):
    with tarfile.open(target, "r:gz") as tar:
        tar.extractall(path)


context = zmq.Context()
socket = context.socket(zmq.SUB)
socket.connect("tcp://127.0.0.1:5556")
print("connected!")
socket.setsockopt_string(zmq.SUBSCRIBE, "")

message = socket.recv_multipart()
time.sleep(10)
print("Received!")

path = "/media/pika/Workbench/workspace/pywgMLST"
jobid = "76824138-ccc9-11e6-9873-dc85de763f2b"
filename = os.path.join(path, jobid + ".tar.gz")
with open(filename, "wb") as f:
    for m in message:
        f.write(m)

# decompress file
decompress(filename, "/home/pika")
