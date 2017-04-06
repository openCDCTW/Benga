from celery import Celery

app = Celery('roaryapp')
app.config_from_object('roaryapp.celeryconfig')

if __name__ == '__main__':
    app.start()

