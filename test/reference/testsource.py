import zmq
import random
import time
import sys

context = zmq.Context()
sender = context.socket(zmq.PUSH)
sender.bind("tcp://*:5557")
print("Sender ready!")

random.seed()
total_msec = 0
for task_nbr in range(100):

    workload = random.randint(1, 100)
    total_msec += workload
    msg = u'{}'.format(workload)

    sys.stdout.write(msg + " sended!\n")
    sys.stdout.flush()

    sender.send_string(msg)
    time.sleep(0.5)
