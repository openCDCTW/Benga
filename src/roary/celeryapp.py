from celery import Celery

app = Celery('roary')
app.config_from_object('roary.celeryconfig')

if __name__ == '__main__':
    app.start()

