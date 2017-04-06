from celery import Celery

app = Celery('prokkaapp')
app.config_from_object('prokkaapp.celeryconfig')

if __name__ == '__main__':
    app.start()

