import zmq
import sys
from threading import Thread
import time
from random import randint

IP = "*"
PORT = 5570


def tprint(msg):
    """like print, but won't get newlines confused with multiple threads"""
    sys.stdout.write(msg + '\n')
    sys.stdout.flush()


class Server(Thread):
    def __init__(self):
        super(Server, self).__init__()

    def run(self):
        context = zmq.Context()
        frontend = context.socket(zmq.ROUTER)
        print("Server started")
        frontend.bind('tcp://{}:{}'.format(IP, PORT))

        backend = context.socket(zmq.DEALER)
        print("Backend started")
        backend.bind('inproc://backend')

        workers = []
        for i in range(5):
            worker = Worker(context)
            worker.start()
            workers.append(worker)

        zmq.proxy(frontend, backend)

        frontend.close()
        backend.close()
        context.term()


class Worker(Thread):
    def __init__(self, context):
        super(Worker, self).__init__()
        self.context = context

    def run(self):
        socket = self.context.socket(zmq.DEALER)
        socket.connect('inproc://backend')
        tprint('Worker started')
        while True:
            ident, msg = socket.recv_multipart()
            tprint('Worker received {} from {}'.format(msg, ident))

            replies = randint(0, 4)
            for i in range(replies):
                time.sleep(1. / (randint(1, 10)))
                time.sleep(5)
                print("send reply: {} {}".format(ident, msg))
                socket.send_multipart([ident, msg])


def main():
    server = Server()
    server.start()
    server.join()

if __name__ == "__main__":
    main()
