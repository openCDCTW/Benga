import os
import time


def check_ready(async_results, outpath):
    deleted = []
    for r in async_results:
        if r.ready():
            for filename, file in r.get():
                print(filename)
                with open(os.path.join(outpath, filename), "w") as f:
                    f.write(file)
            deleted.append(r)
    for d in deleted:
        async_results.remove(d)


def prokka(inpath, outpath):
    print("sleep for 10 sec")
    time.sleep(10)
    results = []
    for name in os.listdir(inpath):
        filename = os.path.join(inpath, name)
        with open(filename, "r") as f:
            result = prokka.s(name, f.read()).apply_async()
        results.append(result)

    time.sleep(60)

    while results:
        time.sleep(10)
        check_ready(results, outpath)


def roary(inpath, outpath, ident_min, threads):
    print("sleep for 10 sec")
    time.sleep(10)
    # prepare uploads
    uploads = []
    for name in os.listdir(inpath):
        filename = os.path.join(inpath, name)
        with open(filename, "r") as f:
            uploads.append((name, f.read()))

    # run and wait for return
    result = roary.s(uploads, ident_min, threads).apply_async()
    while not result.ready():
        time.sleep(10)

    # receive and write result
    filename, file = result.get()
    print(filename)
    with open(os.path.join(outpath, filename), "w") as f:
        f.write(file)

