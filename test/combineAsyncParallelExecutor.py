from concurrent.futures import ThreadPoolExecutor, wait
import asyncio
import time
import sys


def g(i):
    time.sleep(2)
    print("task {} done!".format(i))
    sys.stdout.flush()
    return i


async def foo():
    futures = []
    with ThreadPoolExecutor() as executor:
        for i in range(10):
            print("assign job")
            future = executor.submit(g, i)
            futures.append(future)
            time.sleep(1)

    done, not_done = wait(futures)
    print("DONE: ", [x.result() for x in done])
    print("NOT DONE: ", not_done)


loop = asyncio.get_event_loop()
loop.run_until_complete(foo())
