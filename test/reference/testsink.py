import sys
import zmq

context = zmq.Context()

receiver = context.socket(zmq.PULL)
receiver.bind("tcp://*:5558")

while True:
    s = receiver.recv()
    sys.stdout.write("Result: " + s.decode("ascii") + "\n")
    sys.stdout.flush()

