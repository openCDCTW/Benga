from concurrent.futures import ThreadPoolExecutor
import time


def foo(a):
    x = 0
    for i in range(0, 10000000):
        x *= i
    return a, x

a = [i for i in range(50)]

results = []

t = time.time()
with ThreadPoolExecutor(5) as executor:
    for x in a:
        r = executor.submit(foo, x)
        # r.add_done_callback(lambda y: print(y.result()))
        # r.add_done_callback(lambda y: print(y.result(), "callback 2"))
        results.append(r)
t = time.time() - t

print("results:")
for r in results:
    print(r.result())
print(t, " sec")

