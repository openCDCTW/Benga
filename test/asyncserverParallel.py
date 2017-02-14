import zmq
import sys
from concurrent.futures import ThreadPoolExecutor, wait
import asyncio
import time

IP = "*"
PORT = 5570


def tprint(msg):
    sys.stdout.write(msg + '\n')
    sys.stdout.flush()


def g(i):
    time.sleep(2)
    tprint("task {} done!".format(i))
    return i


async def server():
    context = zmq.Context()

    receiver = context.socket(zmq.PULL)
    receiver.connect("tcp://localhost:5557")

    sender = context.socket(zmq.PUSH)
    sender.connect("tcp://localhost:5558")

    print("Server started")

    futures = []
    print("Backend started")
    with ThreadPoolExecutor() as executor:
        while True:
            msg = receiver.recv()
            msg = msg.decode("ascii")
            print("received message {", msg, "}")
            print("assign job")
            future = executor.submit(g, msg)
            futures.append(future)
            sender.send_string(msg + " submitted!")
            time.sleep(1)

            done, not_done = wait(futures)
            for f in done:
                print(f.result())
                futures.remove(f)
            print("NOT DONE: ", not_done)
            print()


def main():
    loop = asyncio.get_event_loop()
    loop.run_until_complete(server())

if __name__ == "__main__":
    main()
