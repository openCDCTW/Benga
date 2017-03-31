from celery import Celery

app = Celery('prokka')
app.config_from_object('prokka.celeryconfig')

if __name__ == '__main__':
    app.start()

