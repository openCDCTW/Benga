import zmq
import time
import os

from src.utils import files


class Sender:
    def __init__(self, ip, port):
        self._ip = ip
        self._port = port

    """
    compress folder into a tar.gz file and send it to server
    """
    def send(self, filepath, compressed_name):
        print("Establishing socket...")
        self._context = zmq.Context()
        self._socket = self._context.socket(zmq.PUB)

        print("Connecting to server...")
        self._socket.bind("tcp://" + self._ip + ":" + str(self._port))
        time.sleep(1)

        print("Compressing files...")
        files.to_gz(filepath, compressed_name)

        print("Sending file...")
        compressed_file = compressed_name + ".tar.gz"
        with open(compressed_file, "rb") as f:
            self._socket.send_multipart(f.readlines())
            time.sleep(1)

        print("Deleting compressed file...")
        os.remove(compressed_file)

        print("Closing socket...")
        self._socket.close()
        self._context.term()


class Receiver:
    def __init__(self, ip, port):
        self._ip = ip
        self._port = port

    def receive(self):
        message = self._socket.recv_multipart()
        time.sleep(10)
        return message

    def run(self, filepath, jobtype, func):
        print("Establishing socket...")
        self._context = zmq.Context()
        self._socket = self._context.socket(zmq.SUB)
        self._socket.connect("tcp://" + self._ip + ":" + str(self._port))

        print("Connected!")
        self._socket.setsockopt_string(zmq.SUBSCRIBE, "")

        while True:
            print("Waiting for file...")
            message = self.receive()
            print("Received!")

            filename = files.joinpath(filepath, jobtype + ".tar.gz")
            with open(filename, "wb") as f:
                for m in message:
                    f.write(m)

            print("Decompress files...")
            jobid = files.decompress_gz(filename, filepath)

            # create folder
            job_dir = files.joinpath(filepath, jobtype)
            files.create_if_not_exist(job_dir)
            id_dir = files.joinpath(job_dir, jobid)
            files.create_if_not_exist(id_dir)

            func(files.joinpath(filepath, jobid), id_dir)

