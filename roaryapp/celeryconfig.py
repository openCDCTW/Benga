BROKER_URL = 'pyamqp://admin:mypass@my-rabbit:5672'
CELERY_RESULT_BACKEND = 'redis://my-redis:6379/0'
CELERY_INCLUDE = ['roaryapp.tasks']
CELERY_ACCEPT_CONTENT = ['json']
CELERY_TASK_SERIALIZER = 'json'

