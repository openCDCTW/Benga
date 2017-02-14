import sys
import time
import zmq

context = zmq.Context()

receiver = context.socket(zmq.PULL)
receiver.connect("tcp://localhost:5557")

sender = context.socket(zmq.PUSH)
sender.connect("tcp://localhost:5558")

while True:
    s = receiver.recv()

    sys.stdout.write("Work on " + s.decode("ascii") + "\n")
    sys.stdout.flush()
    time.sleep(int(s)*0.001)

    sender.send(s)
