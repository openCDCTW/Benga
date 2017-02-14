import zmq
import sys
from threading import Thread

SERVERIP = "*"
SERVERPORT = 5570


def tprint(msg):
    sys.stdout.write(msg + '\n')
    sys.stdout.flush()


class Client(Thread):
    def __init__(self, id):
        super(Client, self).__init__()
        self.id = id

    def run(self):
        context = zmq.Context()
        receiver = context.socket(zmq.PULL)
        receiver.bind("tcp://*:5558")

        sender = context.socket(zmq.PUSH)
        sender.bind("tcp://*:5557")

        identity = u'worker-{}'.format(self.id)
        receiver.identity = identity.encode('ascii')
        print('Client {} started'.format(identity))

        poll = zmq.Poller()
        poll.register(receiver, zmq.POLLIN)
        # poll.register(sender, zmq.POLLIN)

        reqs = 0
        while True:
            reqs += 1
            print('Req #{} sent..'.format(reqs))
            sender.send_string(u'request #{}'.format(reqs))

            for i in range(5):
                sockets = dict(poll.poll(1000))
                if receiver in sockets:
                    msg = receiver.recv()
                    tprint('Client {} received: {}'.format(identity, msg))

            # msg = receiver.recv()
            # tprint('Client {} received: {}'.format(identity, msg))


def main():
    for i in range(3):
        client = Client(i)
        client.start()

if __name__ == "__main__":
    main()
