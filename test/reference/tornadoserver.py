from concurrent.futures import ThreadPoolExecutor
from tornado import gen, ioloop
from tornado.web import RequestHandler, Application


executor = ThreadPoolExecutor(8)


class ThreadPoolHandler(RequestHandler):

    def _do_loop(self, a):
        x = 0
        for i in range(0, 10000000):
            x *= i
        return a, x

    @gen.coroutine
    def get(self, chunksize=16384):
        self.request.files[""]

        with open(path, 'rb') as f:
            while True:
                data = f.read(chunksize)
                if not data:
                    break
                self.write(data)
        self.finish()


        n = int(self.get_argument("number", default=None, strip=False))
        a = [i for i in range(n)]
        for x in a:
            r = executor.submit(self._do_loop, x)
            # r.add_done_callback(lambda y: print(y.result()))
            yield r
            self.write("{}\n".format(x))
            self.flush()
        # self.render("index.html")
        return self.response.body


if __name__ == "__main__":
    application = Application([
        (r"/", ThreadPoolHandler)
    ])

    application.listen(8888)
    ioloop.IOLoop.current().start()
